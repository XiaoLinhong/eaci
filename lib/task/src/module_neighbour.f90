module mod_neighbour

    use mod_hash, only : hash
    use mod_constant, only : FILLVALUE
    use mod_structure, only : pointInfo

    use mod_tool, only :  cal_distance

    use flogger, only : log_notice, log_warning, log_print

    implicit none
    type scanner
        integer, dimension(:), pointer :: idx
        real, dimension(:), pointer :: ratio
        character(len=16), dimension(:), pointer :: dcode
        character(len=16), dimension(:), pointer :: cityIds
        contains ! 绑定过程
          procedure :: scan
          procedure :: get_dcode
    end type scanner

    contains

    subroutine scan(self, radius, cityLon, cityLat, siteInfo)
        ! 需要哪些点位？
        implicit none
        class(scanner), intent(inout) :: self
        real, intent(in) :: radius
        real, intent(in) :: cityLon
        real, intent(in) :: cityLat
        type(pointInfo),intent(in) :: siteInfo

        ! local
        integer :: d
        integer :: i, j, idx
        type(hash) :: siteLocs
        type(hash) :: distance
        character(len=16) :: thisId

        if (associated(self%idx)) deallocate(self%idx)
        if (associated(self%ratio)) deallocate(self%ratio)
        if (associated(self%dcode)) deallocate(self%dcode)
        if (associated(self%cityIds)) deallocate(self%cityIds)

        call siteLocs%reserve(siteInfo%n)
        call distance%reserve(siteInfo%n)

        do i = 1, siteInfo%n
            d = cal_distance(cityLon, cityLat, siteInfo%lons(i), siteInfo%lats(i))
            if (d < radius) then
                call siteLocs%set(siteInfo%ids(i), i)
                call distance%set(siteInfo%ids(i), d)
            end if
        end do

        allocate( self%idx(siteLocs%n_keys))
        allocate( self%ratio(siteLocs%n_keys))
        allocate( self%cityIds(siteLocs%n_keys))
        j = 0
        do i = 1, siteInfo%n
            thisId = siteInfo%ids(i)
            idx = siteLocs%get(thisId)
            if ( idx > 0 ) then
                j = j + 1
                self%idx(j) = idx
                self%ratio(j) = 1./distance%get(thisId)**0.5
                self%cityIds(j) = thisId(1:6) ! 用来城市平均, 六位的城市编码
            end if
        end do
        call self%get_dcode
    end subroutine scan

    subroutine get_dcode(self)
        ! 获取行政区域编码
        class(scanner), intent(inout) :: self

        ! local var
        integer :: i, j, idx
        character(16), dimension(size(self%cityIds))  :: buffer

        buffer = '-'
        idx = 1
        do i = 1, size(self%cityIds)
            if ( .not. any (buffer(1:idx) == self%cityIds(i)) ) then
                buffer(idx) = self%cityIds(i)
                idx = idx + 1
            end if
        end do
        allocate(self%dcode(idx-1))
        self%dcode(:) = buffer(1:idx-1)
    end subroutine get_dcode

end module mod_neighbour
