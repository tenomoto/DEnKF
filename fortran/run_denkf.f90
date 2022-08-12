program run_denkf
  use, intrinsic :: iso_fortran_env, only: &
    sp => real32, dp => real64
  use denkf_module, only: denkf_init, denkf_clean, denkf_analysis
  implicit none
  
  integer, parameter :: &
    imax = 129, jmax = 129, n = imax * jmax, &
    nens = 25, nobs = 300, dobs = n / nobs, &
    un = 51
  real(kind=dp), parameter :: tau = 6.283185307179586d0
  character(len=*), parameter :: &
    tfile = "true.dat"!, afile = "analysis.dat"
  real(kind=dp), parameter :: &
    rstd=2.0d0, rsig = 4.0d0, &
    inflation_factor = 1.06d0, &
    r0 = 15.d0, r0r0 = r0 * r0

  real(kind=dp), dimension(2) :: u
  real(kind=dp), dimension(:), allocatable :: x, y
  real(kind=dp), dimension(:, :), allocatable :: p12, rho, pf, rmat, hmat
  real(kind=sp), dimension(:), allocatable :: sbuf

  integer :: i, i1, i2, k1, j1, j2, k2, obsloc
  character(len=8) :: efile

  allocate(x(n), p12(n, nens), pf(n, n), sbuf(n), hmat(nobs, n), rmat(nobs, nobs))
  x(:) = 0.0d0
  hmat(:, :) = 0.0d0
  rmat(:, :) = 0.0d0
  do i=1, nobs
    rmat(i, i) = rsig
  end do

  open(unit=un, file=tfile, access="direct", &
    status="old", action="read", recl=n*4)
  read(unit=un, rec=1) sbuf
  close(unit=un)
  call random_number(u(1))
  obsloc = floor(imax * u(1)) + 1
  do i=1, nobs
    print *, obsloc, dobs, sqrt(-2.0d0 * log(1.0d0 - u(1))) * cos(tau * (1.0d0 - u(2)))
    call random_number(u)
    y(i) = sbuf(obsloc) + rstd * &
      sqrt(-2.0d0 * log(1.0d0 - u(1))) * cos(tau * (1.0d0 - u(2)))
    obsloc = obsloc + dobs
    hmat(i, obsloc) = 1.0d0
  end do
  stop

  do i=1, nens
    write(unit=efile, fmt="(i4.4'.dat')") i
    print *, efile
    open(unit=un, file=efile, access="direct", &
      status="old", action="read", recl=n*4)
    read(unit=un, rec=1) sbuf
    close(unit=un)
    p12(:, i) = sbuf(:)
    x(:) = x(:) + sum(p12(:, i))
  end do
  x(:) = x(:) / n
  do i=1, nens
    p12(:, i) = (p12(:, i) - x(:)) / sqrt(nens - 1.0)
  end do

  pf(:, :) = matmul(p12, transpose(p12))
  if (r0 > 0 .and. r0 < 0) then
    k1 = 1
    k2 = 1
    do j2=1, jmax
      do i2=1, imax
        do j1=1, jmax
          do i1=1, imax
            pf(k1, k2) = inflation_factor * exp(-((i1 - i2)**2 + (j1 - j2)**2) / r0r0) * pf(k1, k2)
            k1 = k1 + 1
          end do
        end do
        k2 = k2 + 1
      end do
    end do
  end if
  call denkf_init(n, nobs)
  call denkf_analysis(x, p12, pf, y, hmat, rmat)
  call denkf_clean()

end program run_denkf
