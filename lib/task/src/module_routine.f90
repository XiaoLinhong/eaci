module mod_routine

    use mod_hash, only : hash
    use mod_constant, only : FILLVALUE
    use mod_structure, only : pointInfo
    use mod_structure, only : optMeta

    use mod_neighbour, only : scanner

    use mod_tool, only : cal_distance

    use flogger, only : log_notice, log_warning

    implicit none

    contains

    subroutine get_this_city_date(obsData, obsErr, mdlData, opt, patch, innov, HP, R)

        implicit none
        ! Input Args
        real, dimension(:, :, :, :), intent(in) :: obsData ! nvar, nSite, 24, nDays
        real, dimension(:), intent(in) :: obsErr ! nvar
        real, dimension(:, :, :, :, :), intent(in) :: mdlData ! nvar, nSite, 24, nDays, mDim
        type(optMeta), intent(in) :: opt
        type(scanner), intent(in) :: patch

        ! Out Args
        real, dimension(:, :), allocatable, intent(out) :: innov ! oDim <= nvar, nSite, 24, nDays
        real, dimension(:, :), allocatable, intent(out) :: HP ! oDim, mDim <= nvar, nSite, 24, nDays
        real, dimension(:, :), allocatable, intent(out) :: R  ! oDim, oDim +

        ! local Vars
        integer :: i, j, k, ii, idx, nn
        integer :: iBeg, iEnd
        real :: length ! 特征长度，KM
        real :: decay ! 距离衰减系数
        integer :: oDim, mDim, nVar, nPoint, nSlice, nSite
        integer :: nVaildObs, nVaildMdl, nVaildMean
        real :: thisObs, thisMdl, thisMean, factor
        real, dimension(:), allocatable :: tmp1d
        real, dimension(:, :, :), allocatable :: obs3d
        real, dimension(:, :, :, :), allocatable :: mdl4d

        decay = 1.
        nSlice = size(obsData, 3)/opt%nTime ! 一天中时间分段
        nSite = size(patch%idx)
        mDim = size( mdlData, 5 )
        if (opt%city) then
            nPoint = size(patch%dcode)
            length = 50. ! 单位为km
        else
            nPoint = nSite
            length = 2. ! km
        end if
        oDim = nPoint*opt%nVar*nSlice

        allocate( tmp1d(nSite) )
        allocate( obs3d(nSite, opt%nTime, size(obsData, 4)) )
        allocate( mdl4d(nSite, opt%nTime, size(obsData, 4), mDim) )
        allocate( innov(oDim, 1) )
        allocate( HP(oDim, mDim) )
        allocate( R(oDim, oDim) )
        innov = 0.
        HP = 0.
        R = 0.
        idx = 0
        ! 各个时段*各个变量*各个站点
        do k = 1, nPoint
            do j = 1, opt%nVar
                do i = 1, nSlice ! 时间片段
                    iBeg = (i-1)*opt%nTime+1
                    iEnd = i*opt%nTime
                    idx = idx + 1
                    ! 处理观测空间
                    if (opt%city) then ! 需要求城市平均
                        nn = 0 ! 每个城市有多少个站点
                        do ii = 1, nSite
                            if ( patch%cityIds(ii) == patch%dcode(k) ) then ! 通过城市编号匹配
                                nn = nn + 1
                                tmp1d(nn) = patch%ratio(ii)
                                obs3d(nn, :, :) = obsData(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :)
                                mdl4d(nn, :, :, :) = mdlData(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :, :)
                            end if
                        end do
                        ! 城市中最近一个观测点的距离化系数作为代表系数
                        if (opt%localisation == 2) decay = minval(tmp1d(1:nn)) 
                    else ! 不需要求城市平均
                        if (opt%localisation == 2) decay = patch%ratio(k)
                        nn = 1
                        obs3d(nn, :, :) = obsData( opt%idxs(j), patch%idx(k), iBeg:iEnd, :)
                        mdl4d(nn, :, :, :) = mdlData( opt%idxs(j), patch%idx(k), iBeg:iEnd, :, :)
                    end if
                    factor = opt%ratio(j)*decay

                    ! 观测平均值
                    thisObs = 0.
                    nVaildObs = COUNT(obs3d(1:nn, :, :) /= FILLVALUE)
                    if (nVaildObs > 0) thisObs = sum(obs3d(1:nn, :, :) , obs3d(1:nn, :, :)  /= FILLVALUE)/nVaildObs

                    ! 观测误差矩阵
                    if (nVaildObs < size(obs3d(1:nn, :, :))/3. .or. thisObs <= 0. ) then ! 观测质量不行
                        R(idx, idx) = 9999.
                    else
                        R(idx, idx) = thisObs*obsErr( opt%idxs(j) ) ! 观测误差
                        R(idx, idx) = R(idx, idx) + (opt%delta/length)**0.5*R(idx, idx)/2. ! 代表性误差
                        ! sqr(delta/L)*err, L为代表性误差特征长度, err是估计出来的
                        R(idx, idx) = R(idx, idx) / factor ! 局地化
                    end if

                    ! 集合平均
                    thisMean = 0.
                    nVaildMean = COUNT(mdl4d(1:nn, :, :, :)  /= FILLVALUE)
                    if (nVaildMean > 0) thisMean = sum(mdl4d(1:nn, :, :, :) , mdl4d(1:nn, :, :, :)  /= FILLVALUE)/nVaildMean

                    ! innov = y_o - Hx_b
                    if (thisObs>0. .and. thisMean>0.) innov(idx, 1) = (thisObs - thisMean) !* factor

                    ! 集合成员
                    do ii = 1, mDim
                        thisMdl = 0.
                        nVaildMdl = COUNT(mdl4d(1:nn, :, :, ii) /= FILLVALUE)
                        if (nVaildMdl > 0) thisMdl = sum(mdl4d(1:nn, :, :, ii), mdl4d(1:nn, :, :, ii) /= FILLVALUE)/nVaildMdl
                        ! (X-mean(X))/sqt((m-1))
                        if (thisMdl>0.) HP(idx, ii) = ((thisMdl - thisMean)/(mDim-1.)**0.5)
                    end do
                    if ( sum(HP(idx, :)) > 0.1 ) call log_warning('HP is wrong')
                end do
            end do
        end do

        deallocate(obs3d, mdl4d, tmp1d)

        if (opt%inflation ) HP = get_lambda(innov, HP, R)*HP

    end subroutine get_this_city_date

    real function get_lambda(innov, HP, R) result(lambda)

        implicit none

        ! Input Args
        real, dimension(:, :), intent(in) :: innov ! oDim, 1
        real, dimension(:, :), intent(in) :: HP ! oDim, mDim
        real, dimension(:, :), intent(in) :: R ! oDim, oDim

        ! local varbile
        integer :: i
        real :: denominator
        real, dimension(1, 1) :: numerator
        real, dimension(size(innov), 1) :: MM !R^-1/2 * (y_o - Hx_b)
        real, dimension(size(innov), size(innov)) :: RR ! oDim, oDim
        real, dimension(size(HP, 2), size(HP, 2)) :: PP ! mDim, mDim

        lambda = 1.
        ! 分子
        RR = 0.
        do i = 1, size(innov)
            RR(i, i) = 1./R(i, i)**0.5
        end do
        MM = matmul(RR, innov)
        numerator = matmul(transpose(MM), MM) - size(innov)
        if (numerator(1, 1) <= 0.) return
        ! 分母
        do i = 1, size(innov)
            RR(i, i) = 1./R(i,i)
        end do
        PP = matmul(matmul(transpose(HP), RR), HP)
        ! trace
        denominator = 0.
        do i = 1, size(HP, 2)
            denominator = denominator + PP(i,i)
        end do
        if (denominator <= 0.) return
        ! lambda
        lambda = (numerator(1, 1)/denominator)**0.5
        if (lambda<1.) lambda = 1.
        if (lambda>10.) lambda = 10.

    end function get_lambda

end module mod_routine
