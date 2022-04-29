module mod_enkf

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
        ! eigenvalue decomposition
        if (lowRank) call low_rank(RR)

        ! gainMatrix = P(HP)'RR
        gainMatrix = matmul(matmul(P, transpose(HP)), RR)
        x_b = x_b + matmul(gainMatrix, innov)
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
        real(kind=8), dimension(oDim, oDim) ::  RR
        real(kind=8), dimension(oDim, oDim) ::  ZZ
        real(kind=8), dimension(oDim) ::  EE

        real(kind=8), dimension(oDim*8) :: FWORK
        integer, dimension(oDim*5) :: IWORK
        integer, dimension(oDim) :: IFAIL
        real(kind=8) abstol, ddum
        integer neig, INFO
        real, external :: DLAMCH

        RR = dble(R)
        abstol = 2.0*DLAMCH('S')
        call dsyevx('V', 'A', 'U', oDim, RR, oDim, ddum, ddum, 1, 1, abstol, &
                    neig, EE, ZZ, oDim, FWORK, 8*oDim, IWORK, IFAIL, INFO )
        if (INFO /= 0) stop 'dsyevx ierr'
        Z = ZZ
        eigenvalues = EE
    end subroutine get_eigenvalues

    subroutine low_rank(RR)
        ! RR = Z*EE*Z'
        implicit none
        real, dimension(:, :), intent(inout) :: RR ! oDim, oDim

        integer :: oDim
        real, dimension(size(RR, 1)) :: eigenvalues ! oDim, oDim
        real, dimension(size(RR, 1), size(RR, 1)) :: Z ! oDim, oDim
        real, dimension(size(RR, 1), size(RR, 1)) :: EE ! diag(eigenvalues)
        real :: sum1, sum2
        integer :: i, idx

        oDim = size(RR, 1)
        call get_eigenvalues(oDim, RR, Z, eigenvalues)

        ! Significant eigenvalues
        EE = 0.
        idx = 0
        sum1 = 0.
        sum2 = sum( eigenvalues )
        do i = oDim, 1, -1
            if (sum1/sum2 < 0.99 ) then
                sum1 = sum1 + eigenvalues(i)
                EE(i, i) = 1./eigenvalues(i)
            else
                idx = i
                exit
            end if
        end do
        RR = matmul(matmul(Z(:, idx+1:oDim), EE(idx+1:oDim, idx+1:oDim)), transpose(Z(:, idx+1:oDim)) )
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
