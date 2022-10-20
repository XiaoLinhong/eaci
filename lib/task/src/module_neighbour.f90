module mod_neighbour

    use mod_structure, only : pointInfo

    use mod_tool, only :  cal_distance

    use flogger, only : log_notice, log_warning, log_print

    implicit none
    type scanner
        integer, dimension(:), pointer :: idx ! 检索范围中点 在 站点列表中的位置
        real, dimension(:), pointer :: ratio ! 衰减率
        character(len=15), dimension(:), pointer :: dcode ! 一共包含多少城市
        character(len=15), dimension(:), pointer :: cityIds ! 每个站点对应的城市ID
        contains ! 绑定过程
          procedure :: scan
          procedure :: get_dcode
    end type scanner

    contains

    subroutine scan(self, cityLon, cityLat, radius, length, siteInfo)
        ! 需要哪些点位？
        implicit none
        class(scanner), intent(inout) :: self
        real, intent(in) :: radius
        real, intent(in) :: length
        real, intent(in) :: cityLon
        real, intent(in) :: cityLat
        type(pointInfo),intent(in) :: siteInfo

        ! local
        real :: d
        integer :: i
        integer :: nSite

        real, dimension(:), allocatable :: distance
        integer, dimension(:), allocatable :: siteLocs
        character(15), dimension(:), allocatable :: cityIDs

        ! 资源复用
        if (associated(self%idx)) deallocate(self%idx)
        if (associated(self%ratio)) deallocate(self%ratio)
        if (associated(self%dcode)) deallocate(self%dcode)
        if (associated(self%cityIds)) deallocate(self%cityIds)

        allocate( cityIDs(siteInfo%n) )
        allocate( distance(siteInfo%n) )
        allocate( siteLocs(siteInfo%n) )

        nSite = 0
        do i = 1, siteInfo%n
            d = cal_distance(cityLon, cityLat, siteInfo%lons(i), siteInfo%lats(i))
            if (d < radius) then
                nSite = nSite + 1
                cityIDs(nSite) = siteInfo%cityIds(i)
                distance(nSite) = d
                siteLocs(nSite) = i
            end if
        end do

        allocate( self%idx(nSite))
        allocate( self%ratio(nSite))
        allocate( self%cityIds(nSite))

        do i = 1, nSite
            self%idx(i) = siteLocs(i)
            self%ratio(i) = exp( - distance(i)**2/(2.*length**2) )
            self%cityIds(i) = cityIds(i) ! 用来求城市平均
        end do
        ! 城市编号
        call self%get_dcode
        ! 释放资源
        deallocate(cityIDs)
        deallocate(distance)
        deallocate(siteLocs)

    end subroutine scan

    subroutine get_dcode(self)
        ! 获取行政区域编码
        class(scanner), intent(inout) :: self

        ! local var
        integer :: i, j, idx
        character(15), dimension(size(self%cityIds))  :: buffer

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
