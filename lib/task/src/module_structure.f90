module mod_structure
    ! 一些简单的数据结构

    use mod_constant
    use mod_hash, only : hash

    implicit none

    private

    public pointInfo
    public csvMeta
    public optMeta ! 算法参数
    public param ! namelist

    type csvMeta ! 站点数据
        integer :: nVar
        character(len=LENFILENAME) :: fileName
        character(len=LENVAR), dimension(MAXVAR) :: varNames = '-'
        real, dimension(MAXVAR) :: error = 1.0
    end type csvMeta

    type optMeta ! 排放数据源描述信息
        real :: delta = 0.
        real :: radius = 0.
        real :: length = 0.
        integer :: nTime = 0
        integer :: localisation = 0
        logical :: city = .false.
        logical :: lowRank = .false.
        logical :: inflation = .false.

        integer :: idx
        character(len=LENVAR) :: name = '-'

        integer :: nVar
        integer, dimension(MAXVAR) :: idxs
        real, dimension(MAXVAR) :: ratio = 1.0
        character(len=LENVAR), dimension(MAXVAR) :: varNames = '-'
    end type optMeta

    type param ! 配置参数

        ! share
        logical :: debug = .false.
        integer :: mDim
        integer :: nDay
        integer :: nHour
        character(len=LENVAR) :: begTime ! 开始时间
        character(len=LENFILENAME) :: siteFileName
        character(len=LENFILENAME) :: outFileName
        character(len=LENFILENAME) :: dftFileName = '-'
        
        ! file
        type(csvMeta) :: obsInfo
        type(csvMeta) :: mdlInfo
        type(csvMeta) :: adjInfo
        character(len=LENFILENAME) :: rawFileName = '-'

        ! enkf
        integer :: nVar
        type(optMeta), dimension(MAXVAR) :: opts

    end type param

    type pointInfo ! 配置参数
        integer :: n
        type(hash) :: locs
        character(len=15), dimension(:), allocatable :: ids ! 站点号
        real, dimension(:), allocatable :: lats
        real, dimension(:), allocatable :: lons
    end type pointInfo

end module mod_structure
