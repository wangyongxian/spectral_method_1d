#define rkind 4

module cheb_spect_method
    implicit none
!    integer, parameter :: rkind = kind(1.0)
    real(rkind), parameter :: PI = 4.0 * atan(1.0)

contains

    subroutine get_deri1_matrix(n, deri)
        implicit none
        integer, intent(in) :: n
        real(rkind), dimension(0:n, 0:n), intent(out) :: deri

        integer :: k, p

        deri(:, :) = 0.0

        do k = 0, n
            do p = k+1, n, 2
                deri(k, p) = 2.0 * p       ! eq (6.3.32) of "A first course in CFD"
            end do
        end do
        deri(0, :) = deri(0, :) * 0.5
    end subroutine get_deri1_matrix

    subroutine get_deri2_matrix(n, deri)
        implicit none
        integer, intent(in) :: n
        real(rkind), dimension(0:n, 0:n), intent(out) :: deri
        integer :: k, p

        deri(:, :) = 0.0

        do k = 0, n
            do p = k+2, n, 2
                deri(k, p) = p * (p * p - k * k)
            end do
        end do
        deri(0, :) = deri(0, :) * 0.5
    end subroutine get_deri2_matrix

    subroutine get_gauss_lobatto_node(n, x)
        implicit none
        integer, intent(in) :: n
        real(rkind), dimension(0:n), intent(out) :: x

        integer :: i

        do i = 0, n
            x(i) = cos(PI * i / n)      ! eq (6.3.22) of "A first course in CFD"
        end do
    end subroutine get_gauss_lobatto_node

    subroutine eval_cheb_upto_nth_order(n, x, fval)
        implicit none
        integer, intent(in) :: n
        real(rkind), intent(in) :: x
        real(rkind), dimension(0:n), intent(out) :: fval
        integer :: i

        fval(0) = 1.0
        fval(1) = x
        do i = 2, n
            fval(i) = 2.0 * x * fval(i-1) - fval(i-2)
        end do
    end subroutine eval_cheb_upto_nth_order

    subroutine get_spectrum_of_f(n, fun, source)
        implicit none
        integer, intent(in) :: n
        real(rkind), dimension(0:n), intent(out) :: source
        interface
            pure elemental function fun(x)
                real(rkind), intent(in) :: x
                real(rkind) :: fun
            end function fun
        end interface

        real(rkind), dimension(0:n, 0:n) :: mat
        integer :: i, j

        ! eq (2.4.15) of "Spectral Methods Fundamentals in Single Domains"
        do i = 0, n
        do j = 0, n
            mat(i, j) = 2.0 / n * cos(PI * i * j / n)
        end do
        end do
        mat(0, :) = mat(0, :) * 0.5
        mat(n, :) = mat(n, :) * 0.5
        mat(:, 0) = mat(:, 0) * 0.5
        mat(:, n) = mat(:, n) * 0.5

        ! eq (6.3.23) of "A first course in CFD"
        do i = 0, n
            source(i) = fun( cos( PI * i / n ) )
        end do

        source(:) = matmul(mat(:,:), source(:))
    end subroutine get_spectrum_of_f

    subroutine gauss_elimination(n, a, b)
        implicit none
        integer, intent(in) :: n
        real(rkind), dimension(0:n, 0:n), intent(inout) :: a
        real(rkind), dimension(0:n), intent(inout) :: b

        integer :: i, j, k
        real(rkind) :: d

        do k = 0, n-1
            do i = k+1, n
                d = a(i, k) / a(k, k)
                do j = k+1, n
                    a(i, j) = a(i, j) - d * a(k, j)
                end do
                b(i) = b(i) - d * b(k)
            end do
        end do

        call solve_utriangle_system(n, a, b)
    end subroutine gauss_elimination

    subroutine solve_utriangle_system(n, a, b)
        implicit none
        integer, intent(in) :: n
        real(rkind), dimension(0:n, 0:n), intent(in) :: a
        real(rkind), dimension(0:n), intent(inout) :: b

        integer :: i, j

        do j = n, 0, -1
            b(j) = b(j) / a(j, j)
            do i = 0, j-1
                b(i) = b(i) - a(i, j) * b(j)
            end do
        end do
    end subroutine solve_utriangle_system
end module cheb_spect_method

module case01
    use cheb_spect_method
    implicit none

    !---------------------------------------------------------------------------
    !
    !   Let us consider the 1-D second-order linear (P)DE proglem
    !
    !       u"(x) - 4 u'(x) + 4 u(x) = exp(x) + C,  for x in [-1, 1] ........(1)
    !
    !   with the Dirichlet boundary conditions
    !
    !       u(-1) = 0  and u(1) = 0 .........................................(2)
    !
    !   and where C is a constant : C = -4 e / (1 + e^2)
    !
    !   The exact solution of the system (1)-(2) is
    !
    !       u(x) = exp(x) - sinh(1) / sinh(2) * exp(2x) + C / 4
    !
    !   where sinh(x) = (exp(x) - exp(-x)) / 2
    !
    !---------------------------------------------------------------------------

    integer, parameter              :: n = 8
    real(rkind), dimension(0:n, 0:n):: deri1, deri2, coefMat
    real(rkind), dimension(0:n)     :: source

    integer, parameter              :: m = 20
    real(rkind), dimension(0:m)     :: xval, u_calc, u_true
    real(rkind), dimension(0:m, 0:n):: cheb_nodes_val

contains

    subroutine chebyshev_tau_method()
        implicit none
        integer     :: i, j

        coefMat(:, :) = 0.0
        do i = 0, n
            coefMat(i, i) = 1.0
        end do

        ! 一阶导数项
        call get_deri1_matrix(n, deri1)
        ! 二阶导数项
        call get_deri2_matrix(n, deri2)

        coefMat = deri2 - 4.0 * deri1 + 4.0 * coefMat

        ! 右端源项的谱
        call get_spectrum_of_f(n, f, source)

        ! 两个边界条件
        call eval_cheb_upto_nth_order( n, -1.0, coefMat(n-1, :) )
        call eval_cheb_upto_nth_order( n, +1.0, coefMat(n, :) )
        source(n-1) = 0.0
        source(n) = 0.0

        write(*,*) "L = "
        write(*,*) (coefMat(i, :), i = 0, n)
        write(*,*) "rhs = "
        write(*,*) source

        ! 求解线性方程组
        call gauss_elimination( n, coefMat, source )

        write(*,*) "solution = "
        write(*,*) source

        ! 与真解比对误差
        do i = 0, m
            xval(i) = 2.0 / m * i + (-1.0)
            !xval(i) = cos( PI * i / m )
            call eval_cheb_upto_nth_order(n, xval(i), cheb_nodes_val(i, :))
        end do
        u_true(:) = u( xval(:) )
        u_calc(:) = matmul( cheb_nodes_val, source )

        write(*,*) "#x, u_calc, u_true, u_error"
        do i = 0, m
            write(*,*) xval(i), u_calc(i), u_true(i), u_true(i) - u_calc(i)
        end do
    end subroutine chebyshev_tau_method

    pure elemental function f(x) result(y)
        implicit none
        real(rkind), intent(in) :: x
        real(rkind)             :: y, c, e

        e = exp(1.0)
        c = - 4.0 * e / (1 + e * e)
        y = exp(x) + c
    end function f

    pure elemental function u(x) result(y)
        implicit none
        real(rkind), intent(in) :: x
        real(rkind)             :: y, c, e, c0

        e = exp(1.0)
        c = - 4.0 * e / (1 + e * e)
        c0 = (exp(1.0) - exp(-1.0)) / (exp(2.0) - exp(-2.0))
        y = exp(x) - c0 * exp(x + x) + 0.25 * c
    end function u
end module case01

module case02
    use cheb_spect_method
    implicit none

    !---------------------------------------------------------------------------
    !
    !   Let us consider the 1-D second-order linear (P)DE proglem
    !
    !       u"(x) = 2,  for x in [-1, 1] ........(1)
    !
    !   with the Dirichlet boundary conditions
    !
    !       u(-1) = 0  and u(1) = 0 .............(2)
    !
    !   The exact solution of the system (1)-(2) is
    !
    !       u(x) = x * x - 1
    !
    !---------------------------------------------------------------------------

    integer, parameter              :: n = 4
    real(rkind), dimension(0:n, 0:n):: deri2, coefMat
    real(rkind), dimension(0:n)     :: source

    integer, parameter              :: m = 10
    real(rkind), dimension(0:m)     :: xval, u_calc, u_true
    real(rkind), dimension(0:m, 0:n):: cheb_nodes_val

contains

    subroutine chebyshev_tau_method()
        implicit none
        integer     :: i, j

        ! 二阶导数项
        call get_deri2_matrix(n, deri2)

        ! 右端源项的谱
        call get_spectrum_of_f(n, f, source)

        coefMat(2:n, :) = deri2(0:n-2, :)
        source(2:n) = source(0:n-2)

        ! 两个边界条件
        call eval_cheb_upto_nth_order( n, -1.0, coefMat(0, :) )
        call eval_cheb_upto_nth_order( n, +1.0, coefMat(1, :) )
        source(0) = 0.0
        source(1) = 0.0

        write(*,*) "L = "
        write(*,*) (coefMat(i, :), i = 0, n)
        write(*,*) "rhs = "
        write(*,*) source

        ! 求解线性方程组
        call gauss_elimination( n, coefMat, source )

        write(*,*) "solution = "
        write(*,*) source

        ! 与真解比对误差
        do i = 0, m
            xval(i) = 2.0 / m * i + (-1.0)
            !xval(i) = cos( PI * i / m )
            call eval_cheb_upto_nth_order(n, xval(i), cheb_nodes_val(i, :))
        end do
        u_true(:) = u( xval(:) )
        u_calc(:) = matmul( cheb_nodes_val, source )

        write(*,*) "#x, u_calc, u_true, u_error"
        do i = 0, m
            write(*,*) xval(i), u_calc(i), u_true(i), u_true(i) - u_calc(i)
        end do
    end subroutine chebyshev_tau_method

    pure elemental function f(x) result(y)
        implicit none
        real(rkind), intent(in) :: x
        real(rkind)             :: y

        y = 2.0
    end function f

    pure elemental function u(x) result(y)
        implicit none
        real(rkind), intent(in) :: x
        real(rkind)             :: y

        y = x * x - 1.0
    end function u
end module case02

program cheb_spect_method_driver
    use case02
    implicit none

    call chebyshev_tau_method()

end program cheb_spect_method_driver
