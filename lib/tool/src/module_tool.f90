module mod_tool
    ! 一些公共小程序

    implicit none
    real, parameter :: FILLVALUE = -999.

    contains

    recursive function replace_string(str, pattern, replace) result(res)
        ! 字符串替换, 检索所有的可能
        implicit none
        character(*), intent(in) :: str
        character(*), intent(in) :: pattern
        character(*), intent(in) :: replace
        character(:), allocatable :: res

        integer i

        i = index(str, pattern)
        if (i == 0) then
            allocate(character(len_trim(str))::res)
            res = str
            return
        end if
        allocate(character((len_trim(str)-len(pattern)+len_trim(AdjustL(replace))))::res)

        res = str(1:i-1) // trim(AdjustL(replace)) // str(i+len(pattern):len_trim(str))

        res = replace_string(res, pattern, replace)

    end function replace_string

    function get_filename(str, site, mm, hh, domain, flag, source, model, mech) result(filename)
        ! 字符串替换, 检索所有的可能
        implicit none
        character(*), intent(in) :: str
        character(*), optional, intent(in) :: site
        character(*), optional, intent(in) :: mm
        character(*), optional, intent(in) :: hh
        character(*), optional, intent(in) :: domain
        character(*), optional, intent(in) :: flag
        character(*), optional, intent(in) :: source
        character(*), optional, intent(in) :: model
        character(*), optional, intent(in) :: mech

        character(:), allocatable :: filename

        allocate(character(len_trim(str))::filename)
        filename = str
        if (present(site)) filename = replace_string(filename, '[SITE]', site)
        if (present(domain)) filename = replace_string(filename, '[DOMAIN]', domain)
        if (present(mm)) filename = replace_string(filename, '[MM]', mm)
        if (present(hh)) filename = replace_string(filename, '[HH]', hh)
        if (present(flag)) filename = replace_string(filename, '[FLAG]', flag)
        if (present(source)) filename = replace_string(filename, '[SOURCE]', source)
        if (present(model)) filename = replace_string(filename, '[MDL]', model)
        if (present(mech)) filename = replace_string(filename, '[MECH]', mech)

    end function get_filename

    logical function does_file_exist(fileName) result(lstat)
        ! 在哪个高度层？
        implicit none
        character(len=*), intent(in) :: filename

        inquire(file=fileName, exist=lstat)
    end function does_file_exist

    integer function counte_cols(fileName) result(num)
        ! 读取一个文件的列数
        implicit none
        ! Input Args
        character(*), intent(in) :: fileName
        ! local var
        integer :: ierror
        character(32) :: head
        character(32), dimension(256) :: buff
        num = 0
        ierror = 0
        open(67, file=fileName, status='old', action='read')
        do while (ierror == 0)
            read(67, *, iostat=ierror) buff(1:num+1)
            if (num == 1) head = buff(1)
            if (num > 1 .and. buff(1) /= head) then ! 已经读到第二行了
                num = num -1
                exit
            end if
            if (ierror == 0) then
                num = num + 1
                ! rewind(67)
                backspace(67) ! 通过跳行来判断
            end if
        end do
        close(67)
    end function counte_cols

    integer function counte_row(fileName) result(num)
        ! 读取一个文件的行数
        implicit none
        ! Input Args
        character(*), intent(in) :: fileName
        ! local var
        integer :: ierror
        character(32) :: buff
        num = 0
        ierror = 0
        open(67, file=fileName, status='old', action='read')
        do while (ierror == 0)
            read(67, *, iostat=ierror) buff
            num = num + 1
        end do
        num = num - 1
        close(67)
    end function counte_row

    real function get_bes(array, missValue) result(bes)
        !------------------------------------------------
        ! 
        ! Purpose :
        ! To realize the best easy systematic estimator 
        ! (BES; Wonnacott and Wonnacott 1972, section 7.3) 
        ! a robust replace for mean 
        !------------------------------------------------
        !
        ! Profile
        ! 1. sort(arr)
        ! 2. caluete the quantile
        ! 3. bes =  (quantile1 + quantile2*2 +quantile3)/4.
        ! 
        !------------------------------------------------
        implicit none
        real, dimension(:), intent(in) :: array ! 传的都是地址？
        real, intent(in) :: missValue ! 缺省值

        ! local
        real, dimension(:), allocatable :: sorted
        integer :: i, j, nvalid

        ! 
        integer   :: location1                  ! the locatin of 1/4 quantile
        integer   :: location2                  ! the locatin of 1/2 quantile
        integer   :: location3                  ! the locatin of 3/4 quantile

        bes = missValue
        nvalid = count(array /= missValue)
        if (nvalid < 2 ) return

        allocate(sorted(nvalid))
        j = 0
        do i = 1, size(array)
            if (array(i) == missValue) cycle
            j = j + 1
            sorted(j) = array(i)
        end do

        if(nvalid <= 4) then ! 太少没有意义
            bes = sum(sorted)/real(nvalid)
        else
            call qSort(sorted, nvalid)
            if (nvalid < 7) then ! 剔除极端值
                bes = sum(sorted(1:nvalid-1))/real(nvalid-2)
            else
                location1      = int(real(nvalid+1)/4.) ! 靠近左边
                location2      = int(real(nvalid+1)/2.) + 1 ! 靠近右边
                location3      = int(real(nvalid+1)*3./4.) ! 靠近左边
                bes =  (sorted(location1) + sorted(location2)*2 + sorted(location3))/4.
            end if
        end if
        if (bes == 0.) bes = 0.00001 ! 会出现0的情况
        deallocate(sorted)
    end function get_bes

    recursive subroutine qSort(a, na)
    ! DUMMY ARGUMENTS
    implicit none
    integer, intent(in)     :: na
    real, intent(inout) :: a(na)

    ! LOCAL VARIABLES
    integer  :: left, right
    real     :: random
    real     :: pivot
    real     :: temp
    integer  :: marker

    if (na > 1) then

       call random_number(random)
       pivot = a(int(random * real(na - 1)) + 1)
       left = 0
       right = na + 1

       do while (left < right)

          right = right - 1
          do while (a(right) > pivot)
             right = right - 1
          end do

          left = left + 1
          do while (a(left) < pivot)
             left = left + 1
          end do

          if (left < right) then
             temp = a(left)
             a(left) = a(right)
             a(right) = temp
          end if
       end do

       if (left == right) then
          marker = left + 1
       else
          marker = left
       end if

       call qSort(a(:marker - 1), marker - 1)
       call qSort(a(marker:), na - marker + 1)
    end if

    end subroutine qSort

    integer function get_idx_in_str_list(name, names) result(idx)
        ! 找位置！
        implicit none
        character(*), intent(in) :: name
        character(*), dimension(:), intent(in) :: names

        ! local var
        integer :: i

        idx = 0
        if (trim(name) == 'xxx') then
            idx = size(names) + 1
            return
        end if
        do i = 1, size(names)
            if (trim(name) == trim(names(i))) then
                idx = i
                exit
            end if
        end do
    end function

    real function cal_distance(lon1, lat1, lon2, lat2) result(d)
        ! haversine’ formula to calculate the great-circle distance between two points
        implicit none
        real, intent(in) :: lon1
        real, intent(in) :: lat1
        real, intent(in) :: lon2
        real, intent(in) :: lat2

        ! local
        real, parameter :: EARTH_RADIUS_M = 6370000.   ! same as MM5 system
        real, parameter :: PI = 3.1416926
        real, parameter :: RAD_PER_DEG = PI/180.

        real :: deltaLamda, deltaPhi
        real :: phi1, phi2
        real :: a

        phi1 = lat1 * RAD_PER_DEG ! φ, λ in radians
        phi2 = lat2 * RAD_PER_DEG

        deltaLamda = (lon2-lon1) * RAD_PER_DEG

        a = acos( sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2) * cos(deltaLamda) ) * EARTH_RADIUS_M
        d = a/1000. ! in KM
    end function cal_distance

    subroutine normalize_positive_variable(array, coff)
        ! normal distribution
        implicit none
        real, dimension(:), intent(inout) :: array
        real, optional, intent(out) :: coff

        ! local
        integer :: n
        real :: u, u_

        n = count(array/=FILLVALUE)
        if (n<2) return
        u = sum(array, array/=FILLVALUE)/real(n)

        where (array/=FILLVALUE .and. array<u ) array = u + 1 - u/array
        if (present(coff)) then
            u_ = sum(array, array/=FILLVALUE)/real(n)
            ! write(*, *) 'mean: ', u_
            coff = u_ ! 求平均值 
            ! 感觉这个没有必要
            ! where (array/=FILLVALUE) array = array*u/u_ ! 保持均值为1，对扰动系数有用
        end if

    end subroutine normalize_positive_variable

    subroutine inverse_normalize_positive_variable(raw, coff)
        ! inverse normal distribution, 只对均值为1的变量进行逆变换
        implicit none
        real, intent(inout) :: raw
        real, optional, intent(in) :: coff

        ! local
        if (present(coff)) raw = raw*coff ! 感觉这个没有必要
        if (raw<1) raw = 1/(2-raw) ! 均值为1 !!!!
    end subroutine inverse_normalize_positive_variable

end module mod_tool
