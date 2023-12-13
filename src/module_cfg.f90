module mod_cfg
    ! 读取配置文件
    use mod_constant
    use mod_structure, only : param
    use mod_structure, only : csvMeta, optMeta

    use mod_tool, only : get_idx_in_str_list

    use string, only : to_str
    use flogger, only : log_error, log_notice

    implicit none

    private

    public get_cfg

    contains

    function get_input_name() result(input)
        ! 获取配置文件名称
        implicit none
        character(len=LENFILENAME)   :: input     ! 配置文件名称
        integer                      :: fnlen     ! The len of the name of input file. 
        integer                      :: ierror    ! Flag for open file error 
        !  read the input file name
        if (command_argument_count() == 0) then
            input = 'namelist.input'
        else
            call get_command_argument(1, input, fnlen, ierror)
        end if
    end function get_input_name

    type(param) function get_cfg result(p)
        ! 获取配置参数
        implicit none

        ! local 
        integer :: i, j, idx
        integer  :: iHandle ! 句柄 
        integer  :: ierror ! Flag for open file error 
        character(len=LENFILENAME) :: fileName ! 配置文件名称

        ! share
        logical :: debug = .false.
        integer :: mDim
        integer :: nHour
        character(len=LENVAR) :: begTime ! 开始时间
        character(len=LENFILENAME) :: siteFileName
        character(len=LENFILENAME) :: cityFileName
        character(len=LENFILENAME) :: outFileName
        character(len=LENFILENAME) :: dftFileName = '-'
        character(len=LENVAR), dimension(MAXVAR) :: sectorNames = '-' ! 行业名称

        ! source
        type(csvMeta) :: obsInfo
        type(csvMeta) :: mdlInfo
        type(csvMeta) :: adjInfo
        character(len=LENFILENAME) :: rawFileName = '-'

        ! algorithm
        integer :: nTime = 24
        integer :: localisation = 1
        real :: delta = 27.
        real :: radius = 500.
        real :: length = 100.
        logical :: city = .false.
        logical :: lowRank = .false.
        logical :: inflation = .false.
        real :: vmin = 0.1
        real :: vmax = 10.
        real :: r1 = 0.2
        real :: r2 = 1.0
        real :: validR = 0.25
        ! enkf
        type(optMeta), dimension(MAXVAR) :: opts

        NAMELIST /share/ debug, mDim, nHour, begTime, cityFileName, siteFileName, outFileName, dftFileName, sectorNames
        NAMELIST /source/ obsInfo, mdlInfo, adjInfo, rawFileName
        NAMELIST /default/ nTime, localisation, delta, radius, length, city, lowRank, inflation, vmin, vmax, r1, r2, validR

        NAMELIST /custom/ opts

        iHandle = 55
        fileName = get_input_name()
        open(iHandle, file=fileName, status='old', iostat=ierror, form='formatted')
        if(ierror /= 0) call log_error(trim(fileName) // ' not exist ...')

        ! share
        read(iHandle, nml=share)
        p%debug = debug
        p%mDim = mDim
        p%nHour = 24
        p%nDay = nHour/24
        if ( p%nDay == 0 ) then ! 不满一天
            p%nDay = 1
            p%nHour = nHour
        end if
        p%begTime = begTime
        p%siteFileName = siteFileName
        p%cityFileName = cityFileName
        p%outFileName = outFileName
        p%dftFileName = dftFileName
        p%nSector = count(sectorNames /= "-")
        p%sectorNames = sectorNames(1:p%nSector)

        ! source
        read(iHandle, nml=source)
        obsInfo%nVar = count(obsInfo%varNames /= "-")
        mdlInfo%nVar = count(mdlInfo%varNames /= "-")
        adjInfo%nVar = count(adjInfo%varNames /= "-")
        p%obsInfo = obsInfo
        p%mdlInfo = mdlInfo
        p%adjInfo = adjInfo
        p%rawFileName = rawFileName

        ! default: algorithm
        read(iHandle, nml=default)

        ! custom: algorithm
        read(iHandle, nml=custom)
        p%nVar = count(opts%name /= "-")
        do i = 1, p%adjInfo%nVar
            if (opts(i)%nTime == 0) opts(i)%nTime = nTime
            if (p%nHour /= 24) opts(i)%nTime = p%nHour ! 不满一天
            if (opts(i)%delta == 0.)  opts(i)%delta = delta
            if (opts(i)%radius == 0.)  opts(i)%radius = radius
            if (opts(i)%length == 0.)  opts(i)%length = length
            if (opts(i)%vmin == 0.)  opts(i)%vmin = vmin
            if (opts(i)%vmax == 0.)  opts(i)%vmax = vmax
            if (opts(i)%r1 == 0.)  opts(i)%r1 = r1
            if (opts(i)%r2 == 0.)  opts(i)%r2 = r2
            if (opts(i)%validR == 0.)  opts(i)%validR = validR
            if (.not. opts(i)%city)  opts(i)%city = city
            if (.not. opts(i)%lowRank)  opts(i)%lowRank = lowRank
            if (.not. opts(i)%inflation) opts(i)%inflation = inflation
            if (opts(i)%localisation == 0) opts(i)%localisation = localisation
        end do

        ! 变量之间的映射关系
        idx = 0
        do i = 1, p%adjInfo%nVar
            ! 输出调优系数变量和排放变量的对应关系
            if (opts(i)%name /= "-") then
                idx = idx + 1
                p%opts(idx) = opts(i)
                p%opts(idx)%idx = get_idx_in_str_list(opts(i)%name, p%adjInfo%varNames(1:p%adjInfo%nVar))
                p%opts(idx)%nVar = count(opts(i)%varNames /= "-") ! 浓度变量
                do j = 1, p%opts(idx)%nVar ! 每个调优系数变量对应的浓度变量，在obs中的位置
                    p%opts(idx)%idxs(j) = get_idx_in_str_list(opts(i)%varNames(j), p%obsInfo%varNames(1:p%obsInfo%nVar) )
                end do
            end if
        end do

        ! 诊断
        if (p%debug) then
            call log_notice('===================================================================')
            write(*, '(A, 10A10)') 'conc: ', (trim(p%obsInfo%varNames(j)), j = 1, p%obsInfo%nVar)
            write(*, '(A, 10A10)') 'emis: ', (trim(p%adjInfo%varNames(j)), j = 1, p%adjInfo%nVar)
            write(*, '(15A10)') 'emis', 'emis-idx', 'conc', 'conc-idx', 'ratio', 'radius', 'length', 'nTime',  'vmin', 'vmax' , 'r1', 'r2', 'validR'
            do i = 1, p%nVar
                do j = 1, p%opts(i)%nVar
                    write(*, '(15A10)') trim(p%opts(i)%name), to_str(p%opts(i)%idx), &
                    trim(p%opts(i)%varNames(j)), to_str(p%opts(i)%idxs(j)), to_str(p%opts(i)%ratio(j), 2, 3),  &
                    to_str(p%opts(i)%radius, 2, 9), to_str(p%opts(i)%length, 2, 9), to_str(p%opts(i)%nTime), &
                    to_str(p%opts(i)%vmin, 2, 7), to_str(p%opts(i)%vmax, 2, 9),&
                    to_str(p%opts(i)%r1, 2, 7), to_str(p%opts(i)%r2, 2, 9), to_str(p%opts(i)%validR, 2, 9)
                end do
            end do
            call log_notice('===================================================================')
        end if

    end function get_cfg

end module mod_cfg
