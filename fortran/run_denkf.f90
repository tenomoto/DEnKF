program run_denkf
  use, intrinsic :: iso_fortran_env, only: &
    sp => real32, dp => real64
  use denkf_module, only: denkf_analysis
  implicit none
  
  integer, parameter :: un = 51, nvar = 2
  character(len=*), parameter :: &
    lfile = "loc.dat", oofile = "obsoff.dat", var = "qp"

  integer :: i, j, k, imax, jmax, nobs, obsoff, dobs, n, nens, step, nstep
  real(kind=sp), dimension(:), allocatable :: ybuf
  real(kind=dp) :: rstd, rsig, inflation_factor
  real(kind=sp), dimension(:, :), allocatable :: sbuf, rho
  real(kind=dp), dimension(:), allocatable :: y
  integer, dimension(:), allocatable :: obsloc
  real(kind=dp), dimension(:, :), allocatable :: x, rmat
  real(kind=dp), dimension(:, :, :), allocatable :: p12, pf
  character(len=6) :: stepstr
  character(len=8) :: odir, fdir, adir
  character(len=11) :: ofile, ffile, afile
  logical :: localize

  namelist /qg_grid/ imax, jmax
  namelist /denkf_cycle/ inflation_factor, localize, step, nstep, nens, obsoff, fdir, adir
  namelist /denkf_obs/ odir, nobs, rsig

  read(unit=5, nml=qg_grid)
  read(unit=5, nml=denkf_cycle)
  read(unit=5, nml=denkf_obs)
  rstd = sqrt(rsig)

  write(unit=stepstr, fmt="(i0.6)") step
  n = imax * jmax
  print *, "n=", n, " nens=", nens, " nobs=", nobs
  allocate(sbuf(n, nens), x(n, nvar), p12(n, nens, nvar), pf(n, n, nvar), rho(n, n), &
    ybuf(nobs), y(nobs), obsloc(nobs), rmat(nobs, nobs))
  x(:, :) = 0.0d0
  rmat(:, :) = 0.0d0
  do i=1, nobs
    rmat(i, i) = rsig
  end do

  print *, "read obs"
  dobs = n / nobs
  print *, "obsoff=", obsoff
  obsloc(1) = obsoff
  do i=2, nobs
    obsloc(i) = obsloc(i-1) + dobs
  end do
  ofile = "o"//stepstr//".dat"
  open(unit=un, file=trim(odir)//"/"//ofile, access="direct", &
    status="old", action="read", recl=size(ybuf)*4)
  read(unit=un, rec=1) ybuf
  close(unit=un)
  y(:) = ybuf(:)

  print *, "read forecast"
  do k=1, nvar
    ffile = var(k:k)//stepstr//".dat"
    open(unit=un, file=fdir//"/"//ffile, access="direct", &
      status="old", action="read", recl=size(sbuf)*4)
      read(unit=un, rec=1) sbuf
    close(unit=un)
    p12(:, :, k) = sbuf(:, :)
    do i=1, nens
      x(:, k) = x(:, k) + p12(:, i, k)
    end do
    x(:, k) = x(:, k) / nens
    do i=1, nens
      p12(:, i, k) = (p12(:, i, k) - x(:, k)) / sqrt(nens - 1.0)
    end do
  end do

  print *, "calculate prior covariance"
  do j=1, n
    do i=1, n
      pf(i, j, 1) = sum(p12(i, :, 1) * p12(j, :, 2)) ! Pf(q, p)
      pf(i, j, 2) = sum(p12(i, :, 2) * p12(j, :, 2)) ! Pf(p, p)
    end do
  end do
  if (localize) then
    print *, "localize prior covariance"
    open(unit=un, file=lfile, access="direct", &
      status="old", action="read", recl=size(rho)*4)
      read(unit=un, rec=1) rho
    close(unit=un)
    print *, "rho: min=", minval(rho), " max=", maxval(rho)
    do k=1, nvar
      pf(:, :, k) = rho(:, :) * pf(:, :, k)
    end do
  end if
  print *, "inflate prior covariance by", inflation_factor
  pf(:, :, :) = inflation_factor * pf(:, :, :)

  print *, "calculate analysis with DEnKF"
  call denkf_analysis(x, p12, pf, y, obsloc, rmat)
  do i=1, nens
    p12(:, i, :) = x(:, :) + p12(:, i, :) * sqrt(nens - 1.0)
  end do

  print *, "save analysis"
  do k=1, nvar
    afile = var(k:k)//stepstr//".dat"
    open(unit=un, file=adir//"/"//afile, access="direct", &
      status="replace", action="write", recl=n*nens*4)
      write(unit=un, rec=1) real(p12(:, :, k), kind=sp)
    close(unit=un)
  end do

end program run_denkf
