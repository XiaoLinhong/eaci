program main
    ! Emission Adjustment Coefficient Inversion
    use mod_constant, only : FILLVALUE

    ! namelist
    use mod_structure, only : param
    use mod_cfg, only : get_cfg

    ! 点位信息
    use mod_structure, only : pointInfo
    use mod_csvio, only : read_info
    ! I/O
    use mod_csvio, only : read_raw_obs
    use mod_csvio, only : read_raw_mdl
    use mod_csvio, only : read_raw_adj
    use mod_csvio, only : write_data_csv

    ! 临近点位
    use mod_neighbour, only : scanner

    use mod_routine, only : get_this_city_loc
    use mod_routine, only : get_this_city_date
    use mod_routine, only : get_this_city_ratio
    ! 核心算法
    use mod_enkf, only : enkf

    ! 处理文件名
    use mod_tool, only : get_filename
    use mod_tool, only : does_file_exist
    use mod_tool, only : get_idx_in_str_list
    ! 时间模块
    use mod_datetime, only : datetime, create_datetime
    use mod_timedelta, only : timedelta, create_timedelta
    ! 处理字符串
    use string, only : to_str
    ! 提示信息
    use flogger, only : log_warning, log_notice, log_print

    implicit none

    type(param) :: cfg ! 配置信息
    type(pointInfo) :: siteInfo ! 站点信息
    type(pointInfo) :: cityInfo ! 城市信息

    type(scanner) :: thisPatch ! 领域

    integer :: idx ! 当前城市在观测中的位置
    integer :: i, j, k, ii, jj
    ! 待订正变量
    real, dimension(:, :, :), allocatable :: adjMean ! nvar, nCity, nSector
    real, dimension(:, :, :, :), allocatable :: adjData ! nvar, nCity, nSector, mDim
    ! 观测数据
    real, dimension(:, :, :, :), allocatable :: obsData ! nvar, nSite, cfg%nHour, nDay
    ! 模式数据
    real, dimension(:, :, :, :), allocatable :: mdlMean ! nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :, :, :, :), allocatable :: mdlData ! nvar, nSite, cfg%nHour, nDay, mDim
    integer, dimension(:, :, :, :), allocatable :: nValid ! nvar, nSite, cfg%nHour, nDay, 有效数据个数

    type(datetime) :: thisTime
    character(len=256) :: thisFile

    real, dimension(:, :, :), allocatable :: x_a ! old adj

    ! ENKF
    real, dimension(:, :, :), allocatable :: x_b ! nvar, nCity, nSector
    real, dimension(:, :), allocatable :: P  ! 1, mDim
    real, dimension(:, :), allocatable :: innov ! oDim, 1 <= nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :), allocatable :: HP ! oDim, mDim <= nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :), allocatable :: R  ! oDim, oDim
    real :: inflation

    ! obs/model
    real :: ratio
    integer, dimension(:), allocatable :: siteLoc ! 当前城市的站点位置

    cfg = get_cfg()

    ! 读取城市信息
    cityInfo = read_info(cfg%cityFileName, flag='city')

    ! 读取矫正系数: 分行业扰动
    allocate(adjData(cfg%adjInfo%nVar, cityInfo%n, cfg%nSector, cfg%mDim))
    do i = 1, cfg%mDim
        thisFile = get_filename(cfg%adjInfo%fileName, flag=trim(to_str(i, 2)))
        if (cfg%debug .and. i==1) call log_notice(thisFile)
        call read_raw_adj(thisFile, cityInfo%locs, cfg%sectorNames(1:cfg%nSector), adjData(:, :, :, i))
    end do

    ! 扰动平均，通常情况下, 这个都应该为1, 不考虑缺省问题
    allocate(adjMean(cfg%adjInfo%nVar, cityInfo%n, cfg%nSector))
    adjMean = sum(adjData, dim=4)/cfg%mDim ! cfg%nVar, nCity, cfg%nSector

    ! 站点信息： 模式和观测数据站点对齐
    siteInfo = read_info(cfg%siteFileName)

    ! 读取站点数据
    allocate(obsData(cfg%obsInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday))
    obsData = FILLVALUE
    thisTime = create_datetime(cfg%begTime, '%Y%m%d%H')
    do i = 1, cfg%nday
        do j = 1, cfg%nHour
            thisFile = thisTime%format(cfg%obsInfo%fileName)
            if (cfg%debug .and. i==1 .and. j==1) call log_notice(thisFile)
            call read_raw_obs(thisFile, siteInfo%locs, obsData(:, :, j, i))
            thisTime = thisTime + create_timedelta(hours=1)
        end do
    end do

    allocate(mdlData(cfg%mdlInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday, cfg%mDim))
    ! 读取集合预报数据
    do i = 1, cfg%mDim
        thisTime = create_datetime(cfg%begTime, '%Y%m%d%H')
        do j = 1, cfg%nday
            thisFile = get_filename(thisTime%format(cfg%mdlInfo%fileName), flag=trim(to_str(i, 2)))
            if (cfg%debug .and. i==1 .and. j==1) call log_notice(thisFile)
            call read_raw_mdl(thisFile, thisTime%hour, siteInfo%locs, mdlData(:, :, :, j, i))
        end do
    end do

    ! 读取模式预报数据
    allocate(mdlMean(cfg%mdlInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday))
    if (cfg%rawfileName /= '-') then
        thisTime = create_datetime(cfg%begTime, '%Y%m%d%H')
        do j = 1, cfg%nday
            thisFile = get_filename(thisTime%format(cfg%rawfileName), flag=trim(to_str(i, 2)))
            if (cfg%debug .and. j==1) call log_notice(thisFile)
            call read_raw_mdl(thisFile, thisTime%hour, siteInfo%locs, mdlMean(:, :, :, j))
            thisTime = thisTime + create_timedelta(days=j)
        end do
    else
        allocate(nValid(cfg%mdlInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday))
        if (cfg%debug) call log_warning('mean of ensemble will be used to calculate invo!')
        mdlMean = FILLVALUE
        nValid = COUNT(mdlData /= FILLVALUE, DIM=5)
        where( nValid>0 ) mdlMean = sum(mdlData, DIM=5, MASK=mdlData/=FILLVALUE)/nValid
    end if

    allocate( x_b(cfg%nVar, cityInfo%n, cfg%nSector) )
    allocate( P(1, cfg%mDim) ) ! 排放扰动

    ! 累乘权重
    !allocate( x_a(cfg%nVar, cityInfo%n, cfg%nSector) )
    !x_a = 1.
    !if (does_file_exist(cfg%dftFileName)) then
    !    call read_raw_adj(cfg%dftFileName, cityInfo%locs, cfg%sectorNames(1:cfg%nSector), x_a)
    !end if

    x_b = 1.
    do i = 1, cityInfo%n
    !do i = 74, 74 ! 南京市 xiaolh
    !do i = 1, 1 ! 北京市 xiaolh
        if (cfg%debug) call log_notice(cityInfo%ids(i))
        ! 当前城市的站点位置
        call get_this_city_loc(cityInfo%ids(i), siteInfo, siteLoc)
        ! 直辖县级市、澳门、香港、台湾等, 行政区里面可能是没有观测站点的
        if (size(siteLoc) == 0) call log_print('    '//cityInfo%ids(i)// ' has no obs site')
        do j = 1, cfg%nVar ! xiaolh
            ! 扫描目标位置周围的观测点位， 不同变量的检索范围可以不一样
            call thisPatch%scan(cityInfo%lons(i), cityInfo%lats(i), cfg%opts(j)%radius, cfg%opts(j)%length, siteInfo)

            ! 诊断信息
            !if (cfg%debug .and. j==3) write(*, '(15A10)') (trim(thisPatch%dcode(k)), k=1, size(thisPatch%dcode))

            ! 求obs城市浓度/model城市浓度
            call get_this_city_ratio(obsData, mdlMean, cfg%opts(j)%idxs(1), siteLoc, ratio)

            ! 当前城市有的观测站点，但是数据无效
            if ( size(siteLoc) > 0 .and. ratio == 1.0 ) then
                call log_warning('    '//trim(cfg%opts(j)%name)// 'has no valid site!')
                cycle
            end if

            do k = 1, cfg%nSector
                ! 这个污染物可能不扰动
                if (maxval( adjData(cfg%opts(j)%idx, i, k, :) ) > 1.0 ) then 
                    ! 处理目标位置的数据: 局地化，膨胀，过滤缺省值
                    write(*, *) '+++++++++++++++++++++++++++++++++++++++++++++++'
                    call log_print('    '//trim(cfg%opts(j)%name))

                    call get_this_city_date(obsData, cfg%obsInfo%error(1:cfg%obsInfo%nVar), mdlMean, &
                        mdlData, adjData(cfg%opts(j)%idx, i, k, :), cfg%opts(j), thisPatch, innov, HP, R, inflation)

                    if (cfg%opts(j)%inflation) HP = HP*inflation
            
                    ! EnKF
                    if (size(innov) == 0 ) then ! 有数据，但是相关性不够
                        call log_warning('    '//trim(cfg%opts(j)%name)// ' obs is invalid!')
                        cycle
                    end if
            

                    ! 计算排放扰动项 nvar, nCity, nSector, mDim:
                    P = ( adjData(cfg%opts(j)%idx, i:i, k, :) - adjMean(cfg%opts(j)%idx, i, k) )/(cfg%mDim-1.)**0.5 !
                    if (cfg%opts(j)%inflation) P = P*inflation
                    call enkf(x_b(j:j, i:i, k), P, innov, HP, R, cfg%opts(j)%lowRank)

                    ! 限制: 增加鲁棒性, 保证偏差的正态分布，一次调整不能太多！
                    if(x_b(j, i, k) < cfg%opts(j)%vmin) x_b(j, i, k) = cfg%opts(j)%vmin ! 处理极小值
                    if(x_b(j, i, k) > cfg%opts(j)%vmax) x_b(j, i, k) = cfg%opts(j)%vmax ! 处理极大值

                    deallocate(innov, HP, R)
                end if
            end do
        end do
        deallocate(siteLoc)
    end do

    ! 写出去
    call write_data_csv(cfg%outFileName, x_b, cityInfo, cfg%opts(1:cfg%nVar)%name, cfg%sectorNames(1:cfg%nSector))

    ! 累乘权重
    if (does_file_exist(cfg%dftFileName)) then
        allocate( x_a(cfg%nVar, cityInfo%n, cfg%nSector) )
        call read_raw_adj(cfg%dftFileName, cityInfo%locs, cfg%sectorNames(1:cfg%nSector), x_a)
        x_b = x_a * x_b
    end if

    if (trim(cfg%dftFileName) /= '-') call write_data_csv(cfg%dftFileName, x_b, cityInfo, cfg%opts(1:cfg%nVar)%name, cfg%sectorNames(1:cfg%nSector))

end program main
