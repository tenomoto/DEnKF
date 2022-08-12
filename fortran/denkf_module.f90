module denkf_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  real(kind=dp), dimension(:, :), allocatable, private :: &
    hphr, kmat, pf

  public :: denkf_init, denkf_clean, denkf_analysis

contains

  subroutine denkf_init(n, nobs)
    integer, intent(in) :: n, nobs

    allocate(hphr(nobs, nobs), kmat(n, nobs), pf(n, n))

  end subroutine denkf_init

  subroutine denkf_clean()

    deallocate(hphr, kmat, pf)

  end subroutine denkf_clean

  subroutine denkf_analysis(x, p12, pf, y, hmat, rmat)
    use la_module, only: la_inv
    real(kind=dp), dimension(:), intent(inout) :: x
    real(kind=dp), dimension(:, :), intent(inout) :: p12
    real(kind=dp), dimension(:), intent(in) :: y
    real(kind=dp), dimension(:, :), intent(in) :: pf, hmat, rmat

    integer :: nobs

    nobs = size(y)
    hphr(:, :) = matmul(matmul(hmat, pf), transpose(hmat)) + rmat
    call la_inv(hphr)
    kmat(:, :) = matmul(matmul(pf, transpose(hmat)), hphr)
    x = x + matmul(kmat, y - matmul(hmat, x))
    p12 = p12 - 0.5d0 * matmul(kmat, matmul(hmat, p12))

  end subroutine denkf_analysis

end module denkf_module
