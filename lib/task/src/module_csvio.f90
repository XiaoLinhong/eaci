module mod_csvio

    use mod_constant, only : FILLVALUE

    use mod_hash, only : hash
    use mod_structure, only : pointInfo

    use mod_tool, only : counte_row
    use mod_tool, only : get_filename
    use mod_tool, only : does_file_exist
    use mod_tool, only : get_idx_in_str_list

    use mod_datetime, only : datetime, create_datetime
    use mod_timedelta, only : timedelta, create_timedelta

    use flogger, only : log_warning, log_error

    character(len=8), parameter :: YMDH = '%Y%m%d%H' ! 时间格式

    contains

    type(pointInfo) function read_info(vFile, flag) result(info)
        ! 获取点位、城市信息（经纬度）
        implicit none
        ! Input Args
        character(*), intent(in)  :: vFile
        character(*), intent(in), optional  :: flag
        ! local
        integer           :: i
        integer           :: ierror ! Flag for open file error
        character(len=15) :: buff

        if (.not. does_file_exist(vfile)) call log_error(trim(vfile) // ' not exist ...')
        ! 分配空间
        info%n = counte_row(vFile) - 1
        allocate(info%ids(info%n))
        allocate(info%lats(info%n))
        allocate(info%lons(info%n))

        call info%locs%reserve(info%n*2) ! 输出化哈希表空间

        ! 读取文件
        open(57, file=vFile, status='old', iostat=ierror, form='formatted')
        read(57, *, iostat=ierror) ! 头信息
        do i = 1, info%n
            if( present(flag) ) then ! city
                read(57, *, iostat=ierror) info%ids(i), buff, info%lons(i), info%lats(i)
            else ! site
                read(57, *, iostat=ierror) info%ids(i), info%lons(i), info%lats(i)
            end if
            call info%locs%set(info%ids(i), i)
        end do
        close(57)
    end function read_info

    subroutine read_raw_adj(vFile, locs, sectorNames, data3d)
        ! 读取调优系数
        implicit none
        ! Input Args
        character(*), intent(in) :: vFile
        type(hash), intent(inout) :: locs
        character(*), dimension(:), intent(in) :: sectorNames

        ! Out Args
        real, dimension(:, :, :), intent(out) :: data3d ! nvar, nSite, nSector

        ! local
        integer :: i, ierror
        integer :: iSector, iCity
        character(len=16), dimension(7) :: buff
        character(len=16) :: sectorName, thisId
        real, dimension(size(data3d, 1) ) :: data1d ! 一行记录

        data3d = FILLVALUE
        if (.not. does_file_exist(vfile) ) then 
            call log_warning(trim(vfile) // ' not exist ...')
            return
        end if
        open(50, file=vFile, form='formatted', iostat=ierror, status='old')
        read(50, *, iostat=ierror)
        do while (.true.) ! 一直读到最后一行
            ! FLAG, ID,   sector, name,    lon,      lat, MM, HH,  PM25,   PMC,    CO,   NOx,   SO2, NMVOC,   NH3
            read(50, *, iostat=ierror) buff(1), thisId,  sectorName, buff(3:7), data1d
            if(ierror /= 0) exit
            iCity = locs%get(thisId)
            iSector = get_idx_in_str_list(sectorName, sectorNames)
            if ( iCity>0 .and. iSector>0) data3d(:, iCity, iSector) = data1d
        end do
        close(50)

    end subroutine read_raw_adj

    subroutine read_raw_obs(vFile, locs, data2d)
        ! 读取观测数据
        implicit none
        ! Input Args
        character(*), intent(in)  :: vFile
        type(hash), intent(inout) :: locs

        ! Out Args
        real, dimension(:, :), intent(out) :: data2d ! nvar, nSite

        integer :: idx, ierror
        character(len=15) :: thisId
        real, dimension(size(data2d, 1) ) :: data1d ! 一行记录

        data2d = FILLVALUE

        if (.not. does_file_exist(vfile) ) then 
            call log_warning(trim(vfile) // ' not exist ...')
            return
        end if

        open(50, file=vFile, form='formatted', iostat=ierror, status='old')

        do while (.true.) ! 一直读到最后一行
            read(50, *, iostat=ierror) thisId, data1d
            if(ierror /= 0) exit
            call control_quality(data1d)
            idx = locs%get(thisId)
            if ( idx>0 ) data2d(:, idx) = data1d
        end do
        close(50)
    end subroutine read_raw_obs

    subroutine read_raw_mdl(vFile, nhour, locs, data3d)
        ! 读取模式数据
        implicit none
        ! Input Args
        character(*), intent(in)  :: vFile
        integer, intent(in) :: nhour ! 锁定到某个小时, 一般为0，代表读取24个小时
        type(hash), intent(inout) :: locs

        ! Out Args
        real, dimension(:, :, :), intent(out) :: data3d ! nvar, nSite, 24

        ! local Vars
        integer :: idx
        integer :: ierror
        integer :: iHour ! 第几个小时
        type(datetime) :: bTime ! 文件里面的数据时间
        type(timedelta) :: delta
        character(len=15) :: thisId ! 站点号
        character(len=10) :: buff, iTime ! 文件里面的时间字符串
        real, dimension(size(data3d, 1)) :: data1d ! 每个站点， 每个时次的所有变量的数据

        data3d = FILLVALUE

        if (.not. does_file_exist(vfile)) then 
            call log_warning(trim(vfile) // ' not exist ...')
            return
        end if

        open(50, file=vFile, form='formatted', iostat=ierror, status='old')
        write(*, *) trim(vfile)
        do while (.true.) ! 一直读到最后一行
            read(50, *, iostat=ierror) buff, iTime, thisId, data1d
            if(ierror /= 0) exit
            ! 时间控制
            if (bTime%year == 1 ) bTime = create_datetime(iTime, YMDH)
            delta = create_datetime(iTime, YMDH) - bTime! 计算时效
            iHour = int(delta%total_hours()) + 1 - nhour
            if (iHour > size(data3d, 3) ) exit
            call control_quality(data1d)
            idx = locs%get(thisId)
            if (idx > 0 .and. iHour > 0 ) data3d(:, idx, iHour) = data1d
        end do
        close(50)

    end subroutine read_raw_mdl

    subroutine control_quality(data1d)
        ! 控制数据质量, 变量顺序固定
        implicit none
        real, dimension(:), intent(inout) :: data1d

        ! 沙尘
        if (data1d(2) > 200 .and. data1d(1) < 100) data1d(1:2) = FILLVALUE
        if (data1d(2) > 400 .and. data1d(1) < 150) data1d(1:2) = FILLVALUE
        ! 浓度倒挂
        if (data1d(2) < data1d(1) .and. data1d(1) > 0. ) data1d(2) = data1d(1) + 1.0
        ! PM10 => PMC
        if ( data1d(2) > 0. .and. data1d(1) > 0. ) data1d(2) = data1d(2) - data1d(1)
        ! 质控
        where(data1d > 1000.) data1d = FILLVALUE
        ! CO: mg => ug
        data1d(3) = data1d(3)*1000.
        where(data1d <  0.) data1d = FILLVALUE
        where(data1d == 0.) data1d = 0.01 ! 尽量不要出现为0的情况

    end subroutine control_quality

    subroutine write_data_csv(vFile, data, info, varNames, sectorNames)
        ! 写出调优系数表
        implicit none
        ! Input Args
        character(*), intent(in)   :: vFile
        real, dimension(:, :, :), intent(in) :: data
        type(pointInfo), intent(in)       :: info ! 城市信息
        character(*), dimension(:), intent(in) :: varNames
        character(*), dimension(:), intent(in) :: sectorNames

        ! local
        integer           :: i, j
        integer           :: unit ! Open unit 
        integer           :: ierror    ! Flag for open file error 
  
        open(99, file=trim(vFile), form='formatted', iostat=ierror)
        write(99, 103) 'FLAG,       ID,   sector,   name,       lon,       lat, MM,  HH', (ADJUSTR(trim(varNames(i))), i=1, size(varNames))
        do j = 1, size(sectorNames)
            do i = 1, info%n
                write(99, 104) '2,', trim(info%ids(i)),',', trim(sectorNames(j)),',', 'name', info%lons(i), info%lats(i), '01,', '00',  data(:, i, j)
            end do
        end do
        103 format( A63, 10(:,',',A10) )
        104 format( A4, A10, A1, A9, A1, A7, 2(:,',',F10.2), ',', A4, A4, 10(:,',',F10.3) )
        close(99)
    end subroutine write_data_csv

end module mod_csvio
