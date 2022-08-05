module ode_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  public :: ode_rk4

contains

  function ode_rk4(f, x, dt, params) result(xo)
    real(kind=dp), dimension(:, :), intent(in) :: x
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:), intent(in) :: params

    interface

      function f(x, params) result(res)
        use, intrinsic :: iso_fortran_env, only: dp => real64

        real(kind=dp), dimension(:, :), intent(in) :: x
        real(kind=dp), dimension(:), intent(in) :: params
        real(kind=dp), dimension(:, :), allocatable :: res
      end function f

    end interface

    real(kind=dp), dimension(:, :), allocatable :: &
      k1, k2, k3, k4, xo
    integer :: n
 
    n = size(x, 1)
    allocate(k1(n, n), k2(n, n), k3(n, n), k4(n, n), xo(n, n))
    k1 = dt * f(x, params)
    k2 = dt * f(x + 0.5d0 * k1, params)
    k3 = dt * f(x + 0.5d0 * k2, params)
    k4 = dt * f(x + k3, params)
    xo = (k1 + 2.0d0 * (k2 + k3) + k4) / 6.0d0 

  end function ode_rk4

end module ode_module
