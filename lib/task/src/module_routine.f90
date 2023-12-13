module mod_routine

    use mod_constant, only : FILLVALUE
    use mod_structure, only : pointInfo
    use mod_structure, only : optMeta

    use mod_neighbour, only : scanner

    use mod_tool, only : cal_distance

    use flogger, only : log_notice, log_warning

    implicit none

    contains

    subroutine get_this_city_date(obsData, obsErr, mdlMean, mdlData, PE, opt, patch, innov, HP, R, inflation)

        implicit none
        ! Input Args
        real, dimension(:, :, :, :), intent(in) :: obsData ! nvar, nSite, 24, nDays
        real, dimension(:), intent(in) :: obsErr ! nvar
        real, dimension(:, :, :, :), intent(in) :: mdlMean ! nvar, nSite, 24, nDays
        real, dimension(:, :, :, :, :), intent(in) :: mdlData ! nvar, nSite, 24, nDays, mDim
        real, dimension(:) :: PE ! mDim: emission coefficient

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
        real, parameter :: LENGHT = 2. ! 特征长度，KM
        real :: corr ! 相关系数
        real :: thisRate ! 排放能解释浓度多大的变化 < 1.0
        real :: decay ! 距离衰减系数
        integer :: oDim, mDim, nVar, nPoint, nSlice, nSite
        integer :: nVaildObs, nVaildMdl, nVaildMean
        real :: thisObs, thisMdl, thisMean, thisRaw, factor
        real, dimension(:), allocatable :: wgts
        real, dimension(:, :), allocatable :: HX
        real, dimension(:, :, :), allocatable :: obs3d
        real, dimension(:, :, :), allocatable :: raw3d
        real, dimension(:, :, :, :), allocatable :: mdl4d

        real, dimension(:, :), allocatable :: innov_A ! oDim <= nvar, nSite, 24, nDays
        real, dimension(:, :), allocatable :: HP_A ! oDim, mDim <= nvar, nSite, 24, nDays
        real, dimension(:, :), allocatable :: R_A  ! oDim, oDim
    
        nSlice = size(obsData, 3)/opt%nTime ! 一天中，时间分段
        nSite = size(patch%idx)
        mDim = size( mdlData, 5 )
        nPoint = nSite
        if (opt%city) nPoint = size(patch%dcode) ! 学要求城市平均
        oDim = nPoint*opt%nVar*nSlice ! 一个天中多个城市多个变量的偏差，同时订正一个城市

        ! write(*, *) nSite, opt%nTime, size(obsData, 4), mDim
        allocate( wgts(nSite) )
        allocate( obs3d(nSite, opt%nTime, size(obsData, 4)) )
        allocate( raw3d(nSite, opt%nTime, size(obsData, 4)) )
        allocate( mdl4d(nSite, opt%nTime, size(obsData, 4), mDim) )
        allocate( HX(oDim, mDim) )

        ! 包含缺失值
        allocate( innov_A(oDim, 1) )
        allocate( HP_A(oDim, mDim) )
        allocate( R_A(oDim, oDim) )
        R_A = 0.
        idx = 0 ! 有效数据

        ! 各时段*各变量*各站点(或者各城市)
        do k = 1, nPoint
            decay = 1. ! 初始化衰减系数
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
                                wgts(nn) = patch%ratio(ii)
                                obs3d(nn, :, :) = obsData(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :)
                                raw3d(nn, :, :) = mdlMean(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :)
                                mdl4d(nn, :, :, :) = mdlData(opt%idxs(j), patch%idx(ii), iBeg:iEnd, :, :)
                            end if
                        end do
                        ! 城市中最近一个观测点的距离化系数作为代表系数
                        if (opt%localisation == 2) decay = maxval(wgts(1:nn)) 
                        ! write(*, *) nn, patch%dcode(k), decay
                    else ! 不需要求城市平均
                        if (opt%localisation == 2) decay = patch%ratio(k)
                        nn = 1
                        obs3d(nn, :, :) = obsData( opt%idxs(j), patch%idx(k), iBeg:iEnd, :)
                        raw3d(nn, :, :) = mdlMean( opt%idxs(j), patch%idx(k), iBeg:iEnd, :)
                        mdl4d(nn, :, :, :) = mdlData( opt%idxs(j), patch%idx(k), iBeg:iEnd, :, :)
                    end if
                    factor = opt%ratio(j)*decay ! 变量局地化参数*空间距地化参数

                    ! 观测平均值: 时段整体平均
                    nVaildObs = COUNT(obs3d(1:nn, :, :) /= FILLVALUE)
                    if (nVaildObs > 0) thisObs = sum(obs3d(1:nn, :, :) , obs3d(1:nn, :, :)  /= FILLVALUE)/nVaildObs
                    ! 观测质量不行
                    if (nVaildObs < size(obs3d(1:nn, :, :))/3. .or. thisObs <= 0. ) cycle

                    ! 模式预报: 缺省值应该很少
                    nVaildMean = COUNT(raw3d(1:nn, :, :) /= FILLVALUE)
                    if (nVaildMean > 0) thisRaw = sum(raw3d(1:nn, :, :) , raw3d(1:nn, :, :)  /= FILLVALUE)/nVaildMean
                    ! 模式质量不行
                    if (nVaildObs < size(raw3d(1:nn, :, :))/3. .or. thisRaw <= 0. ) cycle

                    ! 集合平均：不应该有很多缺测
                    nVaildMean = COUNT(mdl4d(1:nn, :, :, :) /= FILLVALUE)
                    if ( nVaildMean < size(mdl4d(1:nn, :, :, :))/2. ) cycle
                    thisMean = sum(mdl4d(1:nn, :, :, :) , mdl4d(1:nn, :, :, :) /= FILLVALUE)/nVaildMean
                    
                    ! 判断相关系数: 去掉噪音
                    corr = get_coefficient(PE, mdl4d(1:nn, :, :, :))
                    if (abs(corr) < opt%validR) cycle ! 相关性弱的，就不考虑！
                    ! 有一些虚假负相关？剔除掉吧
                    if ( corr*decay <  opt%r1 .or. corr*decay > opt%r2 ) cycle
                    ! 模式的数值很小时, 偏差太大了
                    ! 排放能解释模式多大的变化: 模拟变化和误差之间的关系,
                    thisRate = get_change_rate(mdl4d(1:nn, :, :, :), thisObs)

                    if (opt%city) write(*, '(A10, 10F8.4, 10F8.4, I5)') patch%dcode(k), corr, thisRate, i
                    if (.not. opt%city) write(*, '(A10, 10F8.4, 10F8.4, I5)') patch%cityIds(k), corr, thisRate, i
                    ! 观测误差矩阵
                    idx = idx + 1
                    R_A(idx, idx) = thisObs*obsErr( opt%idxs(j) ) ! 观测误差
                    R_A(idx, idx) = R_A(idx, idx) + (opt%delta/(LENGHT*nn))**0.2*R_A(idx, idx)/2. ! 代表性误差
                    ! sqr(delta/L)*err, L为代表性误差特征长度, err是估计出来的
                    ! 局地化
                    R_A(idx, idx) = (R_A(idx, idx) / factor)**2 ! 局地化, 增大观测误差， 减少矫正
                    !R_A(idx, idx) = R_A(idx, idx) / factor ! 局地化, 增大观测误差， 减少矫正

                    !=========== 可以用鲁棒性更强的方式，求innov ！
                    innov_A(idx, 1) = (thisObs - thisRaw)*thisRate !* factor

                    ! 集合成员
                    do ii = 1, mDim
                        HX(idx, ii) = thisMean
                        nVaildMdl = COUNT(mdl4d(1:nn, :, :, ii) /= FILLVALUE)
                        if (nVaildMdl > 0) HX(idx, ii) = sum(mdl4d(1:nn, :, :, ii), mdl4d(1:nn, :, :, ii) /= FILLVALUE)/nVaildMdl
                        ! (X-mean(X))/sqt((m-1))
                         HP_A(idx, ii) = (HX(idx, ii) - thisMean)/((mDim-1.)**0.5)
                    end do
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
        deallocate(innov_A, HP_A, R_A)
        deallocate(obs3d, raw3d, mdl4d, wgts, HX)
    end subroutine get_this_city_date

    real function get_lambda(innov, HX, R) result(lambda)
        ! 计算膨胀系数
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
        if (lambda<1.) lambda = 1.
        if (lambda>10.) lambda = 10.
    end function get_lambda

    subroutine get_this_city_loc(cityID, siteInfo, siteLoc)
        ! ! 当前城市有多少个站点
        implicit none
        character(len=15), intent(in) :: cityID
        type(pointInfo),intent(in) :: siteInfo

        integer, dimension(:), allocatable, intent(inout) :: siteLoc

        integer :: i
        integer :: nsite
        integer, dimension(500) :: tmpCity

        nsite = 0
        do i = 1, siteInfo%n
            if ( trim(siteInfo%cityIds(i)) == trim(cityID) ) then
                    nsite = nsite + 1
                    tmpCity(nsite) = i
            end if
        end do
        allocate(siteLoc(nsite))
        siteLoc = tmpCity(1:nsite)

    end subroutine get_this_city_loc

    subroutine get_this_city_ratio(obsData, mdlMean, varIDX, siteLoc, ratio)
        ! 求矫正城市的浓度平均
        implicit none
        ! Input Args
        real, dimension(:, :, :, :), intent(in) :: obsData ! nvar, nSite, 24, nDays
        real, dimension(:, :, :, :), intent(in) :: mdlMean ! nvar, nSite, 24, nDays
        integer, intent(in) :: varIDX ! opt%idxs(1)
        integer, dimension(:), intent(in) :: siteLoc
    
        ! Out Args
        real, intent(out) :: ratio
    
        integer :: i, num
        integer :: nVaild
        real :: obsConc
        real :: mdlConc
        real :: thisObs
        real :: thisMdl
        real, dimension(size(obsData, 3), size(obsData, 4)) :: data2d ! 24, nDays
        num = 0
        thisObs = 0.
        thisMdl = 0.
        ratio = 1.0
        do i = 1, size(siteLoc)
            ! obs opt%idxs(1)
            data2d = obsData(varIDX, siteLoc(i), :, :)
            nVaild = COUNT(data2d /= FILLVALUE)
            obsConc = FILLVALUE
            if (nVaild > 0) obsConc = sum(data2d, data2d /= FILLVALUE)/nVaild
            if (nVaild < size(data2d)/3. .or. obsConc <= 0. ) cycle
    
            ! model
            data2d = mdlMean(varIDX, siteLoc(i), :, :)
            nVaild = COUNT(data2d /= FILLVALUE)
            mdlConc = FILLVALUE
            if (nVaild > 0) mdlConc = sum(data2d, data2d /= FILLVALUE)/nVaild
            if (nVaild < size(data2d)/3. .or. mdlConc <= 0. ) cycle
            thisObs = thisObs + obsConc
            thisMdl = thisMdl + mdlConc
            num = num + 1
        end do
        ! if (num > 0) write(*, *) thisObs/num, thisMdl/num
        if (num > 0) ratio = thisObs/thisMdl
    end subroutine get_this_city_ratio    

    real function get_coefficient(PE, mdl4d) result(r)
        ! 计算膨胀系数
        implicit none
        ! Input Args
        real, dimension(:), intent(in) :: PE ! mDim: emission coefficient
        real, dimension(:, :, :, :), intent(in) :: mdl4d ! nSite, nTime, nDay, mDim


        ! local varbile
        integer :: i, n
        integer :: nVaildMdl

        real :: denominator
        real :: mean1, mean2
        real, dimension(size(PE)) :: a1, a2

        r = 0.
        n = 0
        a1 = 0.
        a2 = 0.
        do i = 1, size(PE)
            nVaildMdl = COUNT(mdl4d(:, :, :, i) /= FILLVALUE)
            if (nVaildMdl > 0) then
                n = n + 1
                a1(n) = PE(i)
                a2(n) = sum(mdl4d(:, :, :, i), mdl4d(:, :, :, i) /= FILLVALUE)/nVaildMdl
            end if
        end do 

        mean1 = sum(a1)/n
        mean2 = sum(a2)/n
        denominator = sum((a1 - mean1)**2)**0.5 * sum((a2 - mean2)**2)**0.5

        if (denominator /= 0) then
            r = sum((a1 - mean1)*(a2 - mean2)) / denominator
        end if

    end function get_coefficient

    real function get_change_rate(mdl4d, thisObs) result(rate)
        ! 计算膨胀系数
        implicit none
        ! Input Args
        real, dimension(:, :, :, :), intent(in) :: mdl4d ! nSite, nTime, nDay, mDim
        real, intent(in), optional :: thisObs


        ! local varbile
        integer :: i, n
        integer :: nVaildMdl

        real :: thisMean, thisMax
        real, dimension(size(mdl4d, 4)) :: a1

        n = 0
        a1 = 0. ! 初始化
        rate = 0.5 ! 默认排放的变化只能解释一半的偏差问题！

        do i = 1, size(mdl4d, 4)
            nVaildMdl = COUNT(mdl4d(:, :, :, i) /= FILLVALUE)
            if (nVaildMdl > 0) then
                n = n + 1
                a1(n) = sum(mdl4d(:, :, :, i), mdl4d(:, :, :, i) /= FILLVALUE)/nVaildMdl
            end if
        end do

        thisMax = maxval(a1)
        if (thisMax /= 0 .and. n>0) then
            thisMean = sum(a1)/n
            if (present(thisObs)) then
                rate = (sum((a1 - thisMean)**2)**0.5/n)/abs(thisMean - thisObs)
            else
                rate = (sum((a1 - thisMean)**2)**0.5/n)/thisMean
            end if
        end if
        ! 排放对观测的解释程度
        if (rate > 0.95) rate = 0.95
        if (rate < 0.01) rate = 0.01

    end function get_change_rate

end module mod_routine

