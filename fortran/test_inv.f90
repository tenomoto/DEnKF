program test_inv
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use la_module, only: la_inv
  implicit none

! adopted from https://www.nag-j.co.jp/lapack/dgetri.htm
  integer, parameter :: n = 4
  real(kind=dp), dimension(n, n) :: amat, ainv
  integer :: i
  amat = reshape( &
    [ [ 1.80d0,  5.25d0,  1.58d0, -1.11d0], &
      [ 2.88d0, -2.95d0, -2.69d0, -0.66d0], &
      [ 2.05d0, -0.95d0, -2.90d0, -0.59d0], &
      [-0.89d0, -3.80d0, -1.04d0,  0.80d0] ], shape(amat))
  ainv(:, :) = amat(:, :)

  print *, "amat="
  do i=1, n
    print "(4f7.3)", amat(i, :)
  end do
  call la_inv(ainv)
  print *, "ainv="
  do i=1, n
    print "(4f7.3)", ainv(i, :)
  end do
  amat = matmul(amat, ainv)
  print *, "amat @ ainv="
  do i=1, n
    print "(4f7.3)", amat(i, :)
  end do

end program test_inv

