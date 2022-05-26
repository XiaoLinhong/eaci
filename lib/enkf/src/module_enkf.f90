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
        real, dimension(size(innov), size(innov)) :: RR ! oDim, oDim
        real, dimension(size(x_b, 1), size(innov)) :: gainMatrix ! 1, oDim

        ! Evaluate RR = (R + HBH')
        ! RR = R + matmul(HP, transpose(HP))
        ! write(*, *) '===================================== HP(HP)^T:'
        ! call print_matrix(matmul(HP, transpose(HP)))

        RR = R
        oDim = size(HP, 1)
        mDim = size(HP, 2)
        call sgemm('n','t', oDim, oDim, mDim, &
                               1.0, HP, oDim, &
                                    HP, oDim, &
                               1.0, RR, oDim)

        if (.not. lowRank) RR = get_penrose_inv( RR )
        ! 特征值分解: R -> Z*eig*Z` 
        ! 可以对矩阵进行降维度: (R + HBH')^-1 = (R + HP(HP)')^-1 = VE^-1V'
        if (lowRank) call low_rank(RR) ! HP(HP)': 正定矩阵

        ! gainMatrix = P(HP)'RR
        gainMatrix = matmul(matmul(P, transpose(HP)), RR)
        x_b = x_b + matmul(gainMatrix, innov)

        write(*, *) '===================================== EP: '
        call print_matrix(P)
        write(*, *) '===================================== HP: '
        call print_matrix(HP)
        write(*, *) '===================================== innov: '
        call print_matrix(innov)
        write(*, *) '===================================== R: '
        call print_matrix(R)
        write(*, *) '===================================== RR: '
        call print_matrix(RR)
        write(*, *) '===================================== PE*HP: '
        call print_matrix(matmul(P, transpose(HP)))
        write(*, *) '===================================== x_b: '
        write(*, '(10F10.2)') x_b

    end subroutine enkf

    function get_penrose_inv(A) result(A_I)
        ! Abtract: To calculate the Penrose inverse of A
        implicit none
        ! declare the dummy varibles
        real, dimension(:, :), intent(in)  :: A
        real, dimension(size(A,1), size(A,2)) :: A_I

        ! declare the local varibles
        integer :: n
        integer :: INFO ! = 0:  successful exit
        integer, dimension(size(A,1))       :: IPIV ! The pivot indices from DGETRF; for 1<=i<=N
        real(kind=8), dimension(size(A,1))  :: WORK
        real(kind=8), dimension(size(A,1), size(A,2)) :: A_I_d

        A_I = 0.
        A_I_d = dble(A)
        n = size(A, 1)
        call dgetrf( n, n, A_I_d, n, IPIV, INFO ) ! LU分解(Sivan Toledo递归算法)
        ! if (INFO /= 0) stop 'Matrix is numerically singular!'
        if (INFO /= 0) write(*, *) 'Matrix is numerically singular!'
        if (INFO /= 0) return 
        ! computes the inverse of a matrix using the LU factorization computed by DGETRF.
        call dgetri( n, A_I_d, n, IPIV, WORK, n, INFO )
        if(INFO /= 0) return
        A_I = A_I_d

    end function get_penrose_inv

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
        if (INFO /= 0) stop 'ssyevx ierr'
        eigenvalues = EE
    end subroutine get_eigenvalues

    subroutine low_rank(RR)
        ! RR = Z*EE*Z'
        implicit none
        real, dimension(:, :), intent(inout) :: RR ! oDim, oDim

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
        RR = matmul(matmul(Z(:, idx:oDim), EE(idx:oDim, idx:oDim)), transpose(Z(:, idx:oDim)) )
    end subroutine low_rank

    subroutine print_matrix(data)
        ! 获取站点列表
        implicit none
        ! Input Args
        real, dimension(:, :), intent(in) :: data
  
        ! local
        integer :: i
        integer :: nx, ny

        nx = min(20, size(data, 1))

        do i = 1, nx
            if (size(data, 2) > 20) then
                write(*, '(20F8.3)') data(i, 1:20)
            else
                write(*, '(20F8.3)') data(i, :)
            end if
        end do

    end subroutine print_matrix

end module mod_enkf
