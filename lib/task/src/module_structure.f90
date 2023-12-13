module mod_structure
    ! 一些简单的数据结构
    use mod_constant
    use mod_hash, only : hash

    implicit none

    private

    public pointInfo ! 点位信息
    public csvMeta ! 文件信息
    public optMeta ! 算法参数
    public param ! 配置参数

    type csvMeta ! 文件信息
        integer :: nVar ! 变量个数
        character(len=LENFILENAME) :: fileName ! 文件名
        character(len=LENVAR), dimension(MAXVAR) :: varNames = '-' ! 变量名
        real, dimension(MAXVAR) :: error = 1.0 ! 观测误差, 单位 %
    end type csvMeta

    type optMeta ! 算法参数
        real :: delta = 0. ! 模式的网格分辨率
        real :: radius = 0. ! 局地化检索半径
        real :: length = 0. ! 去相关化长度
        real :: vmin = 0. ! 调整的最小值
        real :: vmax = 0. ! 调整的最大值
        real :: r1 = 0. ! 相关系数剔除下界
        real :: r2 = 0. ! 相关系数剔除上界
        real :: validR = 0. ! 有效相关系数 abs(r) > validR
        integer :: nTime = 0 ! 多长时间平均
        integer :: localisation = 0 ! 距地化方案
        logical :: city = .false. ! 是否计算城市平均
        logical :: lowRank = .false. ! 是否进行降维
        logical :: inflation = .false. ! 是否进行膨胀

        integer :: idx ! 排放变量索引
        character(len=LENVAR) :: name = '-' ! 排放变量名

        integer :: nVar ! 浓度变量个数
        integer, dimension(MAXVAR) :: idxs ! 变量索引
        real, dimension(MAXVAR) :: ratio = 1.0 ! 浓度变量相关性：0 - 1
        character(len=LENVAR), dimension(MAXVAR) :: varNames = '-' ! 浓度变量名
    end type optMeta

    type param ! 配置参数

        ! share
        logical :: debug = .false. ! 是否打印诊断信息
        integer :: mDim ! 集合成员个数
        integer :: nDay ! 检索多少天的模拟数据
        integer :: nHour ! 检索多小时的模拟数据 nHour = nDay*24
        character(len=LENVAR) :: begTime ! 开始时间
        ! sectors
        integer :: nSector
        character(len=LENVAR), dimension(MAXVAR) :: sectorNames = '-' ! 行业名称
        ! meta
        character(len=LENFILENAME) :: siteFileName ! 站点
        character(len=LENFILENAME) :: cityFileName ! 站点
        character(len=LENFILENAME) :: outFileName ! 输出文件名
        character(len=LENFILENAME) :: dftFileName = '-' ! 上一次调优的文件, 可以有，可以没有
        ! file
        type(csvMeta) :: obsInfo ! 观测文件信息
        type(csvMeta) :: mdlInfo ! 集合预报文件信息
        type(csvMeta) :: adjInfo ! 扰动系数文件信息
        character(len=LENFILENAME) :: rawFileName = '-' ! 模式预报文件信息，可以有，可以没有
        ! enkf: 需要调整的排放变量
        integer :: nVar
        type(optMeta), dimension(MAXVAR) :: opts
    end type param

    type pointInfo ! 点位参数
        integer :: n ! 站点个数
        type(hash) :: locs ! ids => idx, 站点号，索引位置
        character(len=15), dimension(:), allocatable :: ids ! 站点号
        character(len=15), dimension(:), allocatable :: cityIds ! 城市编号
        real, dimension(:), allocatable :: lats ! 维度
        real, dimension(:), allocatable :: lons ! 经度
    end type pointInfo

end module mod_structure
