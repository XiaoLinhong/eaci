module mod_routine

    use mod_hash, only : hash
    use mod_constant, only : FILLVALUE
    use mod_structure, only : pointInfo
    use mod_structure, only : optMeta

    use mod_neighbour, only : scanner

    use mod_tool, only : cal_distance
    use mod_tool, only : normalize_positive_variable

    use flogger, only : log_notice, log_warning

    implicit none

    contains

    subroutine get_this_city_date(obsData, obsErr, rawData, mdlData, opt, patch, innov, HP, R, inflation)

        implicit none
        ! Input Args
        real, dimension(:, :, :, :), intent(in) :: obsData ! nvar, nSite, 24, nDays
        real, dimension(:), intent(in) :: obsErr ! nvar
        real, dimension(:, :, :, :), intent(in) :: rawData ! nvar, nSite, 24, nDays
        real, dimension(:, :, :, :, :), intent(in) :: mdlData ! nvar, nSite, 24, nDays, mDim
        type(optMeta), intent(in) :: opt
        type(scanner), intent(in) :: patch

        ! Out Args
        real, dimension(:, :), allocatable, intent(out) :: innov ! oDim <= nvar, nSite, 24, nDays
        real, dimension(:, :), allocatable, intent(out) :: HP ! oDim, mDim <= nvar, nSite, 24, nDays
        real, dimension(:, :), allocatable, intent(out) :: R  ! oDim, oDim
        real, intent(out) :: inflation

        ! local Vars
        integer :: i, j, k, ii, idx, nn
        integer :: iBeg, iEnd
        real, parameter :: length = 2. ! 特征长度，KM
        real :: decay ! 距离衰减系数
        integer :: oDim, mDim, nVar, nPoint, nSlice, nSite
        integer :: nVaildObs, nVaildMdl, nVaildMean
        real :: thisObs, thisMdl, thisMean, thisRaw, factor
        real, dimension(:), allocatable :: tmp1d
        real, dimension(:, :), allocatable :: HX
        real, dimension(:, :, :), allocatable :: obs3d
        real, dimension(:, :, :), allocatable :: raw3d
        real, dimension(:, :, :, :), allocatable :: mdl4d

        real, dimension(:, :), allocatable :: innov_A ! oDim <= nvar, nSite, 24, nDays
        real, dimension(:, :), allocatable :: HP_A ! oDim, mDim <= nvar, nSite, 24, nDays
        real, dimension(:, :), allocatable :: R_A  ! oDim, oDim +
    
        nSlice = size(obsData, 3)/opt%nTime ! 一天中时间分段
        nSite = size(patch%idx)
        mDim = size( mdlData, 5 )
        nPoint = nSite
        if (opt%city) nPoint = size(patch%dcode)
        oDim = nPoint*opt%nVar*nSlice

        allocate( tmp1d(nSite) )
        allocate( obs3d(nSite, opt%nTime, size(obsData, 4)) )
        allocate( raw3d(nSite, opt%nTime, size(obsData, 4)) )
        allocate( mdl4d(nSite, opt%nTime, size(obsData, 4), mDim) )
        allocate( HX(oDim, mDim) )
        ! 包含缺失值
        allocate( innov_A(oDim, 1) )
        allocate( HP_A(oDim, mDim) )
        allocate( R_A(oDim, oDim) )
        innov_A = 0.
        HP_A = 0.
        R_A = 0.
        idx = 0
        ! 各时段*各变量*各站点(或者各城市)
        do k = 1, nPoint
            decay = 1.
            do j = 1, opt%nVar
                do i = 1, nSlice ! 时间片段
                    iBeg = (i-1)*opt%nTime+1
                    iEnd = i*opt%nTime
                    ! 处理观测空间
                    if (opt%city) then ! 需要求城市平均
                        nn = 0 ! 每个城市有多少个站点
                        do ii = 1, nSite
                            if ( patch%cityIds(ii) == patch%dcode(k) ) then ! 通过城市编号匹配
                                nn = nn + 1
                                tmp1d(nn) = patch%ratio(ii)
                                obs3d(nn, :, :) = obsData(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :)
                                raw3d(nn, :, :) = rawData(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :)
                                mdl4d(nn, :, :, :) = mdlData(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :, :)
                            end if
                        end do
                        ! 城市中最近一个观测点的距离化系数作为代表系数
                        if (opt%localisation == 2) decay = maxval(tmp1d(1:nn)) 
                    else ! 不需要求城市平均
                        if (opt%localisation == 2) decay = patch%ratio(k)
                        nn = 1
                        obs3d(nn, :, :) = obsData( opt%idxs(j), patch%idx(k), iBeg:iEnd, :)
                        raw3d(nn, :, :) = rawData( opt%idxs(j), patch%idx(k), iBeg:iEnd, :)
                        mdl4d(nn, :, :, :) = mdlData( opt%idxs(j), patch%idx(k), iBeg:iEnd, :, :)
                    end if
                    factor = opt%ratio(j)*decay ! 变量局地化参数*空间距地化参数

                    ! 观测平均值
                    thisObs = 0.
                    nVaildObs = COUNT(obs3d(1:nn, :, :) /= FILLVALUE)
                    if (nVaildObs > 0) thisObs = sum(obs3d(1:nn, :, :) , obs3d(1:nn, :, :)  /= FILLVALUE)/nVaildObs
                    ! 观测质量不行
                    if (nVaildObs < size(obs3d(1:nn, :, :))/3. .or. thisObs <= 0. ) cycle

                    ! 模式预报
                    thisRaw = 0.
                    nVaildMean = COUNT(raw3d(1:nn, :, :) /= FILLVALUE)
                    if (nVaildMean > 0) thisRaw = sum(raw3d(1:nn, :, :) , raw3d(1:nn, :, :)  /= FILLVALUE)/nVaildMean

                    ! 观测误差矩阵
                    idx = idx + 1
                    R_A(idx, idx) = thisObs*obsErr( opt%idxs(j) ) ! 观测误差

                    ! 不采用膨胀的方式，感觉效果也没有很好！
                    ! if (thisRaw < thisObs*0.1 ) R_A(idx, idx) = thisRaw*obsErr( opt%idxs(j) ) ! 模式的值很小
                    R_A(idx, idx) = R_A(idx, idx) + (opt%delta/(length*nn))**0.2*R_A(idx, idx)/2. ! 代表性误差
                    ! sqr(delta/L)*err, L为代表性误差特征长度, err是估计出来的
                    R_A(idx, idx) = R_A(idx, idx) / factor ! 局地化，增大观测误差，减少矫正

                    ! innov = y_o - Hx_b: 应该要考虑一下极值的影响!
                    if (thisObs>0. .and. thisRaw>0.) innov_A(idx, 1) = (thisObs - thisRaw) !* factor

                    ! 集合平均
                    thisMean = 0. ! mdlData(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :, :)
                    nVaildMean = COUNT(mdl4d(1:nn, :, :, :)  /= FILLVALUE)
                    if (nVaildMean > 0) thisMean = sum(mdl4d(1:nn, :, :, :) , mdl4d(1:nn, :, :, :) /= FILLVALUE)/nVaildMean
                    ! 集合成员
                    do ii = 1, mDim
                        thisMdl = 0. ! nvar, nSite, 24, nDays, mDim
                        nVaildMdl = COUNT(mdl4d(1:nn, :, :, ii) /= FILLVALUE)
                        if (nVaildMdl > 0) thisMdl = sum(mdl4d(1:nn, :, :, ii), mdl4d(1:nn, :, :, ii) /= FILLVALUE)/nVaildMdl
                        HX(idx, ii) = thisMdl
                        ! (X-mean(X))/sqt((m-1))
                        if (thisMdl>0.) HP_A(idx, ii) = (thisMdl - thisMean)/((mDim-1.)**0.5)
                    end do
                    if ( sum(HP_A(idx, :))/mDim > 1. ) call log_warning('HP is wrong')
                end do
            end do
        end do
        ! 只保留有效部分
        allocate( innov(idx, 1) )
        allocate( HP(idx, mDim) )
        allocate( R(idx, idx) )
        innov = innov_A(1:idx, :)
        HP = HP_A(1:idx, :)
        R = R_A(1:idx, 1:idx)
        inflation = get_lambda(innov, HX(1:idx, :), R)
        deallocate(obs3d, raw3d, mdl4d, tmp1d, HX)
    end subroutine get_this_city_date

    real function get_lambda(innov, HX, R) result(lambda)

        implicit none

        ! Input Args
        real, dimension(:, :), intent(in) :: innov ! oDim, 1
        real, dimension(:, :), intent(in) :: HX ! oDim, mDim
        real, dimension(:, :), intent(in) :: R ! oDim, oDim

        ! local varbile
        integer :: i
        real :: denominator
        real, dimension(1, 1) :: numerator
        real, dimension(size(innov), 1) :: MM !R^-1/2 * (y_o - Hx_b)
        real, dimension(size(innov), size(innov)) :: RR ! oDim, oDim
        real, dimension(size(HX, 2), size(HX, 2)) :: PP ! mDim, mDim

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
        do i = 1, size(innov) ! 对角矩阵的逆
            RR(i, i) = 1./R(i,i)
        end do
        PP = matmul(matmul(transpose(HX), RR), HX)
        ! trace
        denominator = 0.
        do i = 1, size(HX, 2)
            denominator = denominator + PP(i,i)
        end do

        if (denominator <= 0.) return
        ! lambda
        lambda = (numerator(1, 1)/denominator)**0.5

        ! write(*, *) 'numerator:', numerator
        ! write(*, *) 'denominator:', denominator
        ! write(*, *) 'lambda: ', lambda
        if (lambda<1.) lambda = 1.
        if (lambda>50.) lambda = 50.

        ! write(*, *) 'lambda: ', lambda
    end function get_lambda

end module mod_routine
