#define TYPEPARAM character(15)

module mod_hash

    ! use intrinsic, only : ichar

    implicit none
    ! key and value, 键值对
    type pair
        TYPEPARAM :: key
        integer   :: value = 0
    end type pair

    ! 链表法解决冲突
    type node
        type(pair), private, allocatable :: kv
        type(node), private, pointer     :: next => null()

        contains ! 绑定过程
        procedure :: set => set_node ! 设置键值对
        procedure :: get => get_node ! 获取键值对
        procedure :: clear => clear_node ! 清理节点单元

    end type node

    type hash

        integer :: n_keys = 0
        integer :: n_buckets = 0
        type(node), private, dimension(:), pointer :: buckets => null()

        contains ! 绑定过程
        procedure :: set => set_value ! 设置键值对
        procedure :: get => get_value ! 获取键值对
        procedure :: clear => clear_hash ! 清理节点单元
        procedure :: reserve ! 分配hash table的空间

    end type hash

    contains

    subroutine reserve(self, n_buckets)
        ! 分配hash table的空间
        implicit none
        class(hash), intent(inout) :: self
        integer,     intent(in)    :: n_buckets
        ! local
        integer,     dimension(29) :: sizes
        integer :: i
    
        if (self%n_keys > 0) stop 'Cannot reserve when fhash is not empty.'
    
        sizes = (/5, 11, 23, 47, 97, 199, 409, 823, 1741, 3469, 6949, 14033, &
          & 28411, 57557, 116731, 236897, 480881, 976369,1982627, 4026031, &
          & 8175383, 16601593, 33712729, 68460391, 139022417, 282312799, &
          & 573292817, 1164186217, 2147483647/)
        do i = 1, size(sizes)
          if (sizes(i) >= n_buckets) then
            self%n_buckets = sizes(i)
            allocate(self%buckets(self%n_buckets))
            return
          end if
        end do
    end subroutine

    subroutine set_value(self, key, value)
        ! 插入键值对
        implicit none
        class(hash), intent(inout) :: self

        TYPEPARAM, intent(in) :: key
        integer :: value

        integer :: address
        logical :: new

        address = modulo(hash_value(key), self%n_buckets) + 1

        ! 没有去检查空间够不够
        call self%buckets(address)%set(key, value, new)
    
        if (new) self%n_keys = self%n_keys + 1
      end subroutine

      function hash_value(key) result(hash)
        ! 哈希值， 特别关键
        implicit none
        TYPEPARAM, intent(in) :: key
        integer :: hash
        ! local
        integer :: i
        hash = 0
        do i = 1, len(key)
            hash = ICHAR(key(i:i))*i + hash
        end do
    end function

    recursive subroutine set_node(self, key, value, new)
        ! 设置一个节点
        implicit none
        class(node), intent(inout) :: self
        TYPEPARAM, intent(in) :: key
        integer :: value
        logical, optional, intent(out) :: new

        if ( .not. allocated(self%kv) ) then
            allocate(self%kv)
            self%kv%key = key
            self%kv%value = value
            if (present(new)) new = .true.
        else if (self%kv%key == key) then
            self%kv%value = value
            if (present(new)) new = .false.
        else
            if ( .not. associated(self%next) ) allocate(self%next)
            call self%next%set(key, value, new)
        end if
    end subroutine

    function get_value(self, key) result(value)
        ! 获取value
        implicit none
        class(hash), intent(inout) :: self
        TYPEPARAM, intent(in)      :: key
        integer :: value

        ! local
        integer :: address

        address = modulo(hash_value(key), self%n_buckets) + 1
        value = self%buckets(address)%get(key)
      end function

      recursive function get_node(self, key) result(value)
        ! 获取一个节点
        implicit none
        class(node), intent(inout) :: self
        TYPEPARAM, intent(in) :: key
        integer :: value

        if (.not. allocated(self%kv)) then
          ! Not found. (Initial node in the bucket not set)
            value = 0
        else if (self%kv%key == key) then
            value = self%kv%value
        else if (associated(self%next)) then
            value = self%next%get(key)
        end if
      end function

    recursive subroutine clear_node(self)
        implicit none
        class(node), intent(inout) :: self

        if (associated(self%next)) then
            call self%next%clear()
            deallocate(self%next)
            nullify(self%next)
        end if
    end subroutine

    subroutine clear_hash(self)
        implicit none
        class(hash), intent(inout) :: self
        ! local
        integer :: i
        if ( .not. associated(self%buckets)) return
        do i = 1, size(self%buckets)
          if (associated(self%buckets(i)%next)) then
            call self%buckets(i)%next%clear()
            deallocate(self%buckets(i)%next)
          end if
        end do
        deallocate(self%buckets)
        self%n_keys = 0
        self%n_buckets = 0
    end subroutine

end module mod_hash
