module mod_constant
    ! 公用参数
    implicit none

    integer, parameter :: MAXVAR = 16 ! 变量个数

    integer, parameter :: LENVAR = 32 ! 普通字符串宽度
    integer, parameter :: LENFILENAME= 256 ! 文件名宽度

    real, parameter :: FILLVALUE = -999.

end module mod_constant
