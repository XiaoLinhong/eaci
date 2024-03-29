module mod_enkf
    ! x_a = x_b + BH'(R + HBH')^-1(y_o-Hx_b)

    implicit none

    contains

    subroutine enkf(x_b, P, innov, HP, R, lowRank)
        ! 每个网格点的分析值仅利用周边一定距离内的观测数据进行计算
        implicit none
        ! Input Args
        real, dimension(:, :), intent(inout) :: x_b ! 1, 1 : 排放
        real, dimension(:, :), intent(in) :: innov ! oDim, 1
        real, dimension(:, :), intent(in) :: P ! 1, mDim: 排放
        real, dimension(:, :), intent(in) :: HP ! oDim, mDim: 浓度
        real, dimension(:, :), intent(in) :: R ! oDim, oDim
        logical, intent(in) :: lowRank

        ! local varbile
        ! R + HP(HP)'
        integer :: oDim, mDim
        real, dimension(size(x_b, 1), size(innov)) :: PEHP ! 1, oDim
        real, dimension(size(innov),  size(innov)) :: RR ! oDim, oDim
        real, dimension(size(innov),  size(innov)) :: RR_I ! oDim, oDim
        real, dimension(size(x_b, 1), size(innov)) :: gainMatrix ! 1, oDim

        RR = R
        oDim = size(HP, 1)
        mDim = size(HP, 2)
        ! (R + HP(HP): 模式和观测误差和
        ! RR = R + matmul(HP, transpose(HP))
        call sgemm('n','t', oDim, oDim, mDim, &
                               1.0, HP, oDim, &
                                    HP, oDim, &
                               1.0, RR, oDim)
        ! 特征值分解: R -> Z*eig*Z` 
        ! 可以对矩阵进行降维度: (R + HBH')^-1 = (R + HP(HP)')^-1 = VE^-1V'
        if (lowRank)  RR_I = low_rank(RR) ! HP(HP)': 正定矩阵
        if (.not. lowRank) RR_I = get_penrose_inv( RR ) ! 这是伪逆
        
        ! gainMatrix = P(HP)'RR_I
        PEHP = matmul(P, transpose(HP)) ! 排放和浓度的关系
        gainMatrix = matmul(PEHP, RR_I)
        x_b = x_b + matmul(gainMatrix, innov)
        ! where( x_b > 10.) x_b = 10.
        ! where( x_b < 0.1) x_b = 0.1

        if (.true.) then ! 诊断输出
            write(*, *) '===================================== R: ', shape(R)
            call print_matrix(R)
            write(*, *) '===================================== RR: ', shape(RR)
            call print_matrix(RR)
            write(*, *) '===================================== RR_I: ', shape(RR_I)
            call print_matrix(RR_I)
            write(*, *) '===================================== check inv: '
            call print_matrix( matmul(RR, RR_I) )
            write(*, *) '===================================== PE: ', shape(P)
            call print_matrix(P)
            write(*, *) '===================================== HP: ', shape(HP)
            call print_matrix(HP)
            write(*, *) '===================================== innov: ', shape(innov)
            call print_matrix(innov)
            write(*, *) '===================================== PE*HP: ', shape(PEHP)
            call print_matrix(PEHP)
            write(*, *) '===================================== gainMatrix: ', shape(gainMatrix)
            call print_matrix(gainMatrix)
            write(*, *) '====================================='
            write(*, '(A10, 10F10.2)'), 'x_b:', x_b
        end if 

    end subroutine enkf

    function get_penrose_inv(A) result(A_I)
        ! 
        ! Abtract: To calculate the Penrose inverse of A
        implicit none
        ! declare the dummy varibles
        real, dimension(:, :), intent(in)  :: A
        real, dimension(size(A,1), size(A,2)) :: A_I
    
        ! declare the local varibles    
        integer :: n
        integer,dimension(size(A,1))        :: IPIV ! The pivot indices from DGETRF; for1<=i<=N
        real(kind=8),dimension(size(A,1))   :: WORK
        integer                             :: INFO ! = 0:  successful exit
        real(kind=8),dimension(size(A,1), size(A,2)) :: A_tmp

        integer :: i
        A_tmp = dble(A)
        n = size(A, 1)
        call dgetrf( n, n, A_tmp, n, IPIV, INFO ) ! LU分解(Sivan Toledo递归算法)
        ! computes the inverse of a matrix using the LU factorization computed
        ! by DGETRF.
        call dgetri( n, A_tmp, n, IPIV, WORK, n, INFO ) 

        if (INFO /= 0) return
        A_I = sngl(A_tmp)
       end function get_penrose_inv

    function get_penrose_inv1(A) result(A_I)
        ! Abtract: To calculate the Penrose inverse of A
        implicit none
        ! declare the dummy varibles
        real, dimension(:, :), intent(in)  :: A
        real, dimension(size(A,1), size(A,2)) :: A_I

        ! declare the local varibles
        integer :: n
        integer :: INFO ! = 0:  successful exit
        integer, dimension(size(A,1)) :: IPIV ! The pivot indices from DGETRF; for 1<=i<=N
        real, dimension(size(A,1))  :: WORK

        A_I = 0.
        n = size(A, 1)
        call sgetrf( n, n, A, n, IPIV, INFO ) ! LU分解(Sivan Toledo递归算法)
        ! if (INFO /= 0) stop 'Matrix is numerically singular!'
        if (INFO /= 0) write(*, *) 'Matrix is numerically singular!'
        if (INFO /= 0) return 
        ! computes the inverse of a matrix using the LU factorization computed by DGETRF.
        call sgetri( n, A_I, n, IPIV, WORK, n, INFO )
        if (INFO /= 0) return

    end function get_penrose_inv1


    subroutine get_eigenvalues(oDim, R, Z, eigenvalues)
        ! Compute eigenvalue decomposition of R -> Z*eig*Z` 
        ! Returns eigenvectors and eigenvalues in ascending order.
        implicit none
        ! Input Args
        integer, intent(in) :: oDim
        real, dimension(:, :), intent(in) :: R ! oDim, oDim
        ! Out Args
        real, dimension(:, :), intent(out) :: Z ! oDim, oDim
        real, dimension(:), intent(out) :: eigenvalues ! oDim

        ! Local Var
        real, dimension(oDim, oDim) ::  RR
        real, dimension(oDim) ::  EE

        real, dimension(oDim*8) :: FWORK
        integer, dimension(oDim*5) :: IWORK
        integer, dimension(oDim) :: IFAIL
        real abstol, VL, VU
        integer neig, INFO
        real, external :: SLAMCH

        RR = R
        abstol = 2.0*SLAMCH('S')
        call ssyevx('V', 'A', 'U', oDim, RR, oDim, VL, VU, 1, 1, abstol, &
                    neig, EE, Z, oDim, FWORK, 8*oDim, IWORK, IFAIL, INFO )
        if (INFO /= 0) then
            call print_matrix(R)
            stop 'ssyevx ierr'
        end if
        eigenvalues = EE
    end subroutine get_eigenvalues

    function low_rank(RR) result(RR_I)
        ! RR = Z*EE*Z'
        implicit none
        real, dimension(:, :), intent(in) :: RR ! oDim, oDim
        real, dimension(size(RR, 1), size(RR, 1)) :: RR_I ! oDim, oDim

        integer :: oDim
        real, dimension(size(RR, 1)) :: eigenvalues ! oDim, oDim
        real, dimension(size(RR, 1), size(RR, 1)) :: Z ! oDim, oDim
        real, dimension(size(RR, 1), size(RR, 1)) :: EE ! diag(1./eigenvalues)
        real :: sum1, sum2
        integer :: i, idx

        oDim = size(RR, 1)
        call get_eigenvalues(oDim, RR, Z, eigenvalues)

        ! Significant eigenvalues
        EE = 0.
        idx = 1
        sum1 = 0.
        sum2 = sum( eigenvalues )
        do i = oDim, 1, -1
            if (sum1/sum2 < 0.99) then
            ! if (eigenvalues(i)/sum2 > 0.01 ) then
                ! write(*, *) 'eigenvalues: ', eigenvalues(i)
                sum1 = sum1 + eigenvalues(i)
                EE(i, i) = 1./eigenvalues(i)
            else
                idx = i
                ! write(*, *) 'eigenvalues region: ', eigenvalues(idx), eigenvalues(oDim)
                ! write(*, *) 'Significant eigenvalues: ', oDim-idx, 'of ', oDim
                exit
            end if
        end do
        ! Z*EE*Z'
        RR_I = matmul(matmul(Z(:, idx:oDim), EE(idx:oDim, idx:oDim)), transpose(Z(:, idx:oDim)) )
    end function low_rank

    subroutine print_matrix(data)
        ! 获取站点列表
        implicit none
        ! Input Args
        real, dimension(:, :), intent(in) :: data
  
        ! local
        integer :: i
        integer :: nx, ny

        nx = min(19, size(data, 1))

        do i = 1, nx
            if (size(data, 2) > 19) then
                write(*, '(20F9.3)') data(i, 1:19)
            else
                write(*, '(20F9.3)') data(i, :)
            end if
        end do

    end subroutine print_matrix

end module mod_enkf
