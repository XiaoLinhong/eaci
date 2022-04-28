program main
    ! Emission Adjustment Coefficient Inversion

    ! namelist
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

    use mod_enkf, only : enkf

    use mod_datetime, only : datetime, create_datetime
    use mod_timedelta, only : timedelta, create_timedelta

    use string, only : to_str
    use flogger, only : log_warning, log_notice

    implicit none

    type(param) :: cfg ! 配置信息
    type(pointInfo) :: siteInfo
    type(pointInfo) :: cityInfo

    type(scanner) :: thisPatch

    integer :: i, j, k
    real, dimension(:, :), allocatable :: adjMean ! ratio, nCity
    real, dimension(:, :, :), allocatable :: adjData ! ratio, nCity, mDim
    real, dimension(:, :, :, :), allocatable :: obsData ! nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :, :, :, :), allocatable :: mdlData ! nvar, nSite, cfg%nHour, nDay, mDim

    type(datetime) :: thisTime
    character(len=256) :: thisFile

    ! ENKF
    real, dimension(:, :), allocatable :: x_b ! CITY_NVAR, nDim
    real, dimension(:, :), allocatable :: P  ! 1, mDim
    real, dimension(:, :),    allocatable :: innov ! oDim, 1 <= nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :), allocatable :: HP ! oDim, mDim <= nvar, nSite, cfg%nHour, nDay
    real, dimension(:, :), allocatable :: R  ! oDim, oDim

    cfg = get_cfg()

    ! 站点信息
    siteInfo = read_info(cfg%siteFileName)

    ! 读取站点数据
    allocate(obsData(cfg%obsInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday))
    allocate(mdlData(cfg%obsInfo%nVar, siteInfo%n, cfg%nHour, cfg%nday, cfg%mDim))

    thisTime = create_datetime(cfg%begTime, '%Y%m%d%H')
    do i = 1, cfg%nday
        do j = 1, cfg%nHour
            thisFile = thisTime%format(cfg%obsInfo%fileName)
            if (cfg%debug .and. i==1) call log_notice(thisFile)
            call read_raw_obs(thisFile, siteInfo%locs, obsData(:, :, j, i))
            thisTime = thisTime + create_timedelta(hours=1)
        end do
    end do

    ! 读取集合预报数据
    do i = 1, cfg%mDim
        do j = 1, cfg%nday
            thisTime = create_datetime(cfg%begTime, '%Y%m%d%H') + create_timedelta(days=j-1)
            thisFile = thisTime%format(cfg%mdlInfo%fileName)
            thisFile = get_filename(thisFile, flag=trim(to_str(i, 2)))
            if (cfg%debug .and. i==1) call log_notice(thisFile)
            call read_raw_mdl(thisFile, thisTime%hour, siteInfo%locs, mdlData(:, :, :, j, i))
        end do
    end do

    ! 读取矫正系数
    thisFile = get_filename(cfg%adjInfo%fileName, flag=trim(to_str(1, 2)))
    cityInfo = read_info(thisFile, flag='city')
    allocate(adjData(cfg%adjInfo%nVar, cityInfo%n, cfg%mDim))

    do i = 1, cfg%mDim
        thisFile = get_filename(cfg%adjInfo%fileName, flag=trim(to_str(i, 2)))
        if (cfg%debug .and. i==1 ) call log_notice(thisFile)
        call read_raw_adj(thisFile, adjData(:, :, i) )
    end do

    ! 扰动平均
    adjMean = sum(adjData, dim=3)/cfg%mDim ! cfg%nVar, nCity

    allocate( x_b(cfg%nVar, cityInfo%n) )
    x_b = 1.
    allocate( P(1, cfg%mDim) )
    do i = 1, cityInfo%n
        call log_notice(cityInfo%ids(i))
        do j = 1, cfg%nVar 
            call thisPatch%scan(cfg%opts(j)%radius, cityInfo%lons(i), cityInfo%lats(i), siteInfo)
            P = ( adjData(cfg%opts(j)%idx, i:i, :) - adjMean(cfg%opts(j)%idx, i) )/(cfg%mDim-1.)**0.5 ! pe
            call get_this_city_date(obsData, cfg%obsInfo%error(1:cfg%obsInfo%nVar), mdlData, cfg%opts(j), thisPatch, innov, HP, R)
            if (size(innov) > 0) call enkf(x_b(j:j, i:i), P, innov, HP, R, cfg%opts(j)%lowRank)
            deallocate(innov, HP, R)
            write(*, *) cfg%opts(j)%name, shape(innov)
            if (cfg%debug .and. j==1) write(*, *) thisPatch%dcode
            if (cfg%debug .and. j==1) write(*, *) thisPatch%dcode
            if (cfg%debug .and. j==1) write(*, *) 'innov:', shape(innov)
            if (size(innov) == 0 .and. j==1 ) call log_warning(cityInfo%ids(i)// 'has no data!')
        end do
    end do
    where(x_b<0.1) x_b = 0.05
    call write_data_csv(cfg%outFileName, x_b, cityInfo, cfg%opts(1:cfg%nVar)%name)

end program main