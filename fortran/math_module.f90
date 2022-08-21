module math_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  real(kind=dp), parameter, public :: math_tau = 6.283185307179586d0

contains

  function math_random_integer(a, b) result(i)
    integer, intent(in) :: a, b

    integer :: i
    real(kind=dp) :: u

    call random_number(u)
    i = floor((b - a + 1) * u) + a

  end function math_random_integer

  function math_random_uniform(a, b) result(u)
    real(kind=dp), optional, intent(in) :: a, b

    real(kind=dp) :: u

    call random_number(u)
    u = 1.0d0 - u ! (0, 1]
    if (present(a) .and. present(b)) then
      u = (b - a) * u + a
    end if

  end function math_random_uniform

  function math_random_normal() result(x)
    real(kind=dp) :: x

    real(kind=dp), dimension(2) :: u
    call random_number(u)
    u(:) = 1.0d0 - u(:)
    x = sqrt(-2.0d0 * log(u(1))) * cos(math_tau * u(2))

  end function math_random_normal

end module math_module
