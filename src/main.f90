program main
    ! Emission Adjustment Coefficient Inversion

    ! namelist
    use mod_constant, only : FILLVALUE

    use mod_structure, only : param
    use mod_structure, only : pointInfo

    use mod_cfg, only : get_cfg

    use mod_csvio, only : read_info
    use mod_csvio, only : read_raw_obs
    use mod_csvio, only : read_raw_mdl
    use mod_csvio, only : read_raw_adj
    use mod_csvio, only : write_data_csv

    use mod_neighbour, only : scanner

    use mod_routine, only : get_this_city_date

    use mod_tool, only : get_filename
    use mod_tool, only : does_file_exist
    use mod_tool, only : normalize_positive_variable
    use mod_tool, only : inverse_normalize_positive_variable

    use mod_enkf, only : enkf

    use mod_datetime, only : datetime, create_datetime
    use mod_timedelta, only : timedelta, create_timedelta

    use string, only : to_str
    use flogger, only : log_warning, log_notice

    implicit none

    logical, parameter :: norm = .false.
    type(param) :: cfg ! 配置信息
    type(pointInfo) :: siteInfo
    type(pointInfo) :: cityInfo

    type(scanner) :: thisPatch

    integer :: i, j, k, ii, jj
    real, dimension(:, :), allocatable :: adjMean ! nvar, nCity
    real, dimension(:, :), allocatable :: coff ! nvar, nCity ! 转换系数
    real, dimension(:, :, :), allocatable :: adjData ! nvar, nCity, mDim
    real, dimension(:, :, :, :), allocatable :: obsData ! nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :, :, :), allocatable :: rawData ! nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :, :, :, :), allocatable :: mdlData ! nvar, nSite, cfg%nHour, nDay, mDim
    integer, dimension(:, :, :, :), allocatable :: nn ! nvar, nSite, cfg%nHour, nDay, 有效数据个数

    type(datetime) :: thisTime
    character(len=256) :: thisFile

    ! ENKF
    real, dimension(:, :), allocatable :: x_a ! old adj
    real, dimension(:, :), allocatable :: x_b ! CITY_NVAR, nDim
    real, dimension(:, :), allocatable :: P  ! 1, mDim
    real, dimension(:, :), allocatable :: innov ! oDim, 1 <= nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :), allocatable :: HP ! oDim, mDim <= nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :), allocatable :: R  ! oDim, oDim

    cfg = get_cfg()

    ! 读取矫正系数
    thisFile = get_filename(cfg%adjInfo%fileName, flag=trim(to_str(1, 2)))
    cityInfo = read_info(thisFile, flag='city')
    allocate(adjData(cfg%adjInfo%nVar, cityInfo%n, cfg%mDim))
    do i = 1, cfg%mDim
        thisFile = get_filename(cfg%adjInfo%fileName, flag=trim(to_str(i, 2)))
        if (cfg%debug .and. i==1) call log_notice(thisFile)
        call read_raw_adj(thisFile, adjData(:, :, i) )
    end do

    ! 正态化
    if (norm) then 
        allocate(coff(cfg%adjInfo%nVar, cityInfo%n))
        do j = 1, cityInfo%n
            do k = 1, cfg%adjInfo%nVar
                call normalize_positive_variable(adjData(k, j, :), coff(k, j))
                do ii = 1, cfg%nHour
                    do jj = 1, cfg%nday ! 模式平均值可能会为负数！
                        call normalize_positive_variable(mdlData(k, j, ii, jj, :))
                    end do
                end do
            end do
        end do
    end if

    ! 扰动平均，通常清楚下, 这个都应该为1
    allocate(adjMean(cfg%adjInfo%nVar, cityInfo%n))
    adjMean = sum(adjData, dim=3)/cfg%mDim ! cfg%nVar, nCity

    ! 站点信息
    siteInfo = read_info(cfg%siteFileName)

    ! 读取站点数据
    allocate(obsData(cfg%obsInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday))
    thisTime = create_datetime(cfg%begTime, '%Y%m%d%H')
    do i = 1, cfg%nday
        do j = 1, cfg%nHour
            thisFile = thisTime%format(cfg%obsInfo%fileName)
            if (cfg%debug .and. i==1 .and. j==1) call log_notice(thisFile)
            call read_raw_obs(thisFile, siteInfo%locs, obsData(:, :, j, i))
            thisTime = thisTime + create_timedelta(hours=1)
        end do
    end do

    allocate(mdlData(cfg%obsInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday, cfg%mDim))
    ! 读取集合预报数据
    do i = 1, cfg%mDim
        do j = 1, cfg%nday
            thisTime = create_datetime(cfg%begTime, '%Y%m%d%H') + create_timedelta(days=j-1)
            thisFile = thisTime%format(cfg%mdlInfo%fileName)
            thisFile = get_filename(thisFile, flag=trim(to_str(i, 2)))
            if (cfg%debug .and. i==1 .and. j==1) call log_notice(thisFile)
            call read_raw_mdl(thisFile, thisTime%hour, siteInfo%locs, mdlData(:, :, :, j, i))
        end do
    end do

    ! 读取模式预报数据
    allocate(nn(cfg%obsInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday))
    allocate(rawData(cfg%obsInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday))
    if (cfg%rawfileName/='-') then
        do j = 1, cfg%nday
            thisTime = create_datetime(cfg%begTime, '%Y%m%d%H') + create_timedelta(days=j-1)
            thisFile = thisTime%format(cfg%rawfileName)
            thisFile = get_filename(thisFile, flag=trim(to_str(i, 2)))
            if (cfg%debug .and. j==1) call log_notice(thisFile)
            call read_raw_mdl(thisFile, thisTime%hour, siteInfo%locs, rawData(:, :, :, j))
        end do
    else
        if (cfg%debug) call log_warning('mean of ensemble will be used to calculate invo!')
        rawData = FILLVALUE
        nn = COUNT(mdlData /= FILLVALUE, DIM=5)
        where(nn>0) rawData = sum(mdlData, DIM=5)/nn
    end if

    allocate( x_b(cfg%nVar, cityInfo%n) )
    x_b = 1.
    allocate( P(1, cfg%mDim) ) ! 排放扰动
    ! do i = 1, cityInfo%n
    do i = 295, 295 !拉萨
    ! do i = 234, 234 !海口
    ! do i = 1, 1 ! 北京
    ! do i = 123, 123 ! 南昌
    ! do i = 2, 3 ! 天津、唐山
        if (cfg%debug) call log_notice(cityInfo%ids(i))
        ! !$OMP PARALLEL DO PRIVATE(thisPatch, P, innov, HP, R)
        ! do j = 1, cfg%nVar
        do j = 1, 1
            ! 扫描目标位置周围的观测点位，不同变量的检索范围可以不一样
            call thisPatch%scan(cfg%opts(j)%radius, cfg%opts(j)%length, cityInfo%lons(i), cityInfo%lats(i), siteInfo)
            ! write(*, *) thisPatch%dcode
            ! 计算排放扰动项
            P = ( adjData(cfg%opts(j)%idx, i:i, :) - adjMean(cfg%opts(j)%idx, i) )/(cfg%mDim-1.)**0.5 !
            ! 处理目标位置的数据: 局地化，膨胀，过滤缺省值
            call get_this_city_date(obsData, cfg%obsInfo%error(1:cfg%obsInfo%nVar), rawData, mdlData, cfg%opts(j), thisPatch, innov, HP, R)
            if (size(innov) == 0 ) then
                call log_warning(cityInfo%ids(i)// 'has no data!')
                cycle
            end if
            ! EnKF
            if (norm) x_b(j, i) = coff(j, i)
            call enkf(x_b(j:j, i:i), P, innov, HP, R, cfg%opts(j)%lowRank)
            ! x_b(j, i) = inverse_normalize_positive_variable(x_b(j, i), coff(j, i)) ! 逆正态化
            if (norm) call inverse_normalize_positive_variable(x_b(j, i)) ! 逆正态化
            ! 诊断信息
            ! if (cfg%debug .and. j==1) write(*, '(20A10)') (trim(thisPatch%dcode(k)), k=1, size(thisPatch%dcode))
            if (cfg%debug) write(*, *) cfg%opts(j)%name,'oDim: ', to_str(size(innov)), ' ncity: ', to_str(size(thisPatch%dcode))
            deallocate(innov, HP, R)
            if(x_b(j, i) < cfg%opts(j)%vmin) x_b(j, i) = cfg%opts(j)%vmin ! 处理极小值
            if(x_b(j, i) > cfg%opts(j)%vmax) x_b(j, i) = cfg%opts(j)%vmax + (x_b(j, i) - cfg%opts(j)%vmax)**0.5 ! 处理极大值
        end do
        ! !$OMP END PARALLEL DO
    end do

    ! 写出去
    call write_data_csv(cfg%outFileName, x_b, cityInfo, cfg%opts(1:cfg%nVar)%name)
    ! 累乘权重
    if (does_file_exist(cfg%dftFileName)) then
        allocate( x_a(cfg%nVar, cityInfo%n) )
        call read_raw_adj(cfg%dftFileName, x_a)
        x_b = x_a * x_b
        if(x_b(j, i) < cfg%opts(j)%vmin/2.) x_b(j, i) = cfg%opts(j)%vmin/2. ! 处理极小值
        if(x_b(j, i) > cfg%opts(j)%vmax*2.) x_b(j, i) = cfg%opts(j)%vmax*2. + (x_b(j, i) - cfg%opts(j)%vmax*2.)**0.5 ! 处理极大值
    end if
    if (trim(cfg%dftFileName) /= '-') call write_data_csv(cfg%dftFileName, x_b, cityInfo, cfg%opts(1:cfg%nVar)%name)

end program main
