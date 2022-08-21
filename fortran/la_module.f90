module la_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  public :: la_inv

contains

  subroutine la_inv(amat)
    real(kind=dp), dimension(:, :), intent(inout) :: amat

    integer :: n, lda, lwork = -1, info, i 
    integer, dimension(:), allocatable :: ipiv
    real(kind=dp), dimension(1) :: dummy
    real(kind=dp), dimension(:), allocatable :: work

    lda = size(amat, 1)
    n = size(amat, 2)
    allocate(ipiv(n))
    call dgetrf(n, n, amat, lda, ipiv, info)
    call dgetri(n, amat, lda, ipiv, dummy, lwork, info)
    if (info /= 0) then
      print *, "error initializing dgetri: ", info
      stop
    end if
    lwork = nint(dummy(1))
    allocate(work(lwork))
    call dgetri(n, amat, lda, ipiv, work, lwork, info)

    deallocate(work)

  end subroutine la_inv

end module la_module
