module fd_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  public :: fd_laplacian, fd_jacobian

contains

  pure function fd_laplacian(p, i, j) result(lp)
    real(kind=dp), dimension(:, :), intent(in) :: p
    integer, intent(in) :: i, j
    real(kind=dp) :: lp

    lp = p(i+1, j) + p(i-1, j) + p(i, j+1) + p(i, j-1) - 4.0d0 * p(i, j)

  end function fd_laplacian

  pure function fd_jacobian(p, q, i, j) result(ja)
    real(kind=dp), dimension(:, :), intent(in) :: p, q
    integer, intent(in) :: i, j
    real(kind=dp) :: ja, j1, j2, j3

    j1 = (q(i, j+1) - q(i, j-1)) * (p(i+1, j) - p(i-1, j)) &
       - (q(i+1, j) - q(i-1, j)) * (p(i, j+1) - p(i, j-1))
    j2 = (q(i+1, j+1) - q(i+1, j-1)) * p(i+1, j) &
       - (q(i-1, j+1) - q(i-1, j-1)) * p(i-1, j) &
       - (q(i+1, j+1) - q(i-1, j+1)) * p(i, j+1) &
       + (q(i+1, j-1) - q(i-1, j-1)) * p(i, j-1)
    j3 = q(i, j+1) * (p(i+1, j+1) - p(i-1, j+1)) &
       - q(i, j-1) * (p(i+1, j-1) - p(i-1, j-1)) &
       - q(i+1, j) * (p(i+1, j+1) - p(i+1, j-1)) &
       + q(i-1, j) * (p(i-1, j+1) - p(i-1, j-1))
    ja = (j1 + j2 + j3) / 12.0d0

  end function fd_jacobian

end module fd_module
