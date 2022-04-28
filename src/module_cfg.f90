module mod_cfg
    ! 读取配置文件
    use mod_constant
    use mod_structure, only : param
    use mod_structure, only : csvMeta, optMeta

    use mod_tool, only : get_idx_in_str_list

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
        character(len=LENFILENAME) :: outFileName

        ! source
        type(csvMeta) :: obsInfo
        type(csvMeta) :: mdlInfo
        type(csvMeta) :: adjInfo

        ! algorithm
        integer :: nTime
        integer :: localisation
        real :: delta
        real :: radius
        logical :: city = .false.
        logical :: lowRank = .false.
        logical :: inflation = .false.
    
        ! enkf
        type(optMeta), dimension(MAXVAR) :: opts

        NAMELIST /share/ debug, mDim, nHour, begTime, siteFileName, outFileName
        NAMELIST /source/ obsInfo, mdlInfo, adjInfo
        NAMELIST /default/ nTime, localisation, delta, radius, city, lowRank, inflation

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
        p%outFileName = outFileName

        ! source
        read(iHandle, nml=source)
        obsInfo%nVar = count(obsInfo%varNames /= "-")
        mdlInfo%nVar = count(mdlInfo%varNames /= "-")
        adjInfo%nVar = count(adjInfo%varNames /= "-")
        p%obsInfo = obsInfo
        p%mdlInfo = mdlInfo
        p%adjInfo = adjInfo

        ! default: algorithm
        read(iHandle, nml=default)

        ! custom: algorithm
        read(iHandle, nml=custom)
        p%nVar = count(opts%name /= "-")
        do i = 1, p%nVar
            if (opts(i)%nTime == 0) opts(i)%nTime = nTime
            if (opts(i)%delta == 0.)  opts(i)%delta = delta
            if (opts(i)%radius == 0.)  opts(i)%radius = radius
            if (.not. opts(i)%city)  opts(i)%city = city
            if (.not. opts(i)%lowRank)  opts(i)%lowRank = lowRank
            if (.not. opts(i)%inflation) opts(i)%inflation = inflation
            if (opts(i)%localisation == 0) opts(i)%localisation = localisation
        end do
        p%opts = opts

        ! 数据之间的变量映射关系
        do i = 1, p%adjInfo%nVar
            ! 输出调优系数变量和排放变量的对应关系
            p%opts(i)%idx = get_idx_in_str_list(p%opts(i)%name, p%obsInfo%varNames(1:p%obsInfo%nVar) )
            p%opts(i)%nVar = count(p%opts(i)%varNames /= "-") ! 浓度变量
            do j = 1, p%opts(i)%nVar ! 每个调优系数变量对应的浓度变量，在obs中的位置
                p%opts(i)%idxs(j) = get_idx_in_str_list(p%opts(i)%varNames(j), p%obsInfo%varNames(1:p%obsInfo%nVar) )
            end do
        end do

        ! 诊断
        if (p%debug) then
            do i = 1, p%adjInfo%nVar
                call log_notice(p%adjInfo%varNames(i))
                write(*, '(10A10)') (trim(p%opts(i)%varNames(j)), j = 1, p%opts(i)%nVar)
                write(*, '(10I10)') p%opts(i)%idxs(1:p%opts(i)%nVar)
                write(*, '(10F10.2)') p%opts(i)%ratio(1:p%opts(i)%nVar)
            end do
        end if

    end function get_cfg

end module mod_cfg
