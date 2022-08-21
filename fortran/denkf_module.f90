module denkf_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  public :: denkf_analysis

contains

  subroutine denkf_analysis(x, p12, pf, y, obsloc, rmat)
    use la_module, only: la_inv
    real(kind=dp), dimension(:, :), intent(inout) :: x
    real(kind=dp), dimension(:, :, :), intent(inout) :: p12
    real(kind=dp), dimension(:), intent(in) :: y 
    integer, dimension(:), intent(in) :: obsloc
    real(kind=dp), dimension(:, :, :), intent(in) :: pf
    real(kind=dp), dimension(:, :), intent(in) :: rmat

    integer :: n, nobs, nens, nvar, i, j, k
    real(kind=dp), dimension(:), allocatable :: innov
    real(kind=dp), dimension(:, :), allocatable :: incr, hp, hphr
    real(kind=dp), dimension(:, :, :), allocatable :: ph, kmat

    n = size(x, 1)
    nvar = size(x, 2)
    nobs = size(y)
    nens = size(p12, 2)
    allocate(hphr(nobs, nobs), kmat(n, nobs, nvar), &
      innov(nobs), ph(n, nobs, nvar), incr(n, nvar), hp(nobs, nens))

    print *, "hphr"
!    ph(:, :) = matmul(pf, transpose(hmat))
    do k=1, nvar
      do j=1, nobs
        ph(:, j, k) = pf(:, obsloc(j), k)
      end do
    end do
    print *, "ph(q): min=", minval(ph(:, :, 1)), " max=", maxval(ph(:, :, 1))
    print *, "ph(p): min=", minval(ph(:, :, 2)), " max=", maxval(ph(:, :, 2))
!    hphr(:, :) = matmul(hmat, ph) + rmat
    do i=1, nobs
      hphr(i, :) = ph(obsloc(i), :, 2) + rmat(i, :)
    end do
    print *, "hphr: min=", minval(hphr), " max=", maxval(hphr)
    print *, "inverse"
    call la_inv(hphr)
    print *, "inv(hphr): min=", minval(hphr), " max=", maxval(hphr)
    print *, "kmat"
    do k=1, nvar
      kmat(:, :, k) = matmul(ph(:, :, k), hphr)
    end do
    print *, "kmat(q): min=", minval(kmat(:, :, 1)), maxval(kmat(:, :, 1))
    print *, "kmat(p): min=", minval(kmat(:, :, 2)), maxval(kmat(:, :, 2))
    print *, "analysis"
!    innov(:) = y - matmul(hmat, x)
    do j=1, nobs
      innov(j) = y(j) - x(obsloc(j), 2)
    end do
    print *, "innov: min=", minval(innov), maxval(innov)
    do k=1, nvar
      incr(:, k) = matmul(kmat(:, :, k), innov)
      x(:, k) = x(:, k) + incr(:, k)
    end do
    print *, "incr(q): min=", minval(incr(:, 1)), maxval(incr(:, 1))
    print *, "incr(p): min=", minval(incr(:, 2)), maxval(incr(:, 2))
!    p12 = p12 - 0.5d0 * matmul(kmat, matmul(hmat, p12))
    do i=1, nobs
      hp(i, :) = p12(obsloc(i), :, 2)
    end do
    print *, "pf12(q): min=", minval(p12(:, :, 1)), " max=", maxval(p12(:, :, 1))
    print *, "pf12(p): min=", minval(p12(:, :, 2)), " max=", maxval(p12(:, :, 2))
    do k=1, nvar
      p12(:, :, k) = p12(:, :, k) - 0.5d0 * matmul(kmat(:, :, k), hp)
    end do
    print *, "pa12(q): min=", minval(p12(:, :, 1)), " max=", maxval(p12(:, :, 1))
    print *, "pa12(p): min=", minval(p12(:, :, 2)), " max=", maxval(p12(:, :, 2))

  end subroutine denkf_analysis

end module denkf_module
