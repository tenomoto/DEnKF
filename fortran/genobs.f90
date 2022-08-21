program genobs
  use, intrinsic :: iso_fortran_env, only: &
    sp => real32, dp => real64
  use math_module, only: math_random_integer, math_random_normal

  implicit none

  integer, parameter :: un = 51

  integer :: imax, jmax
  character(len=4) :: init
  integer :: nobs, nstep, tsave, nsave
  real(kind=dp) :: dt
  character(len=8) tdir, odir
  character(len=11) pfile, ofile

  integer, dimension(:), allocatable :: obsoff
  real(kind=sp), dimension(:), allocatable :: sbuf, y
  real(kind=dp) :: rsig, rstd
  integer :: i, j, k, n, iobs, dobs

  namelist /qg_grid/ imax, jmax
  namelist /qg_time/ init, nstep, tsave, nsave, dt
  namelist /denkf_true/ tdir
  namelist /denkf_obs/ odir, nobs, rsig

  read(unit=5, nml=qg_grid)
  read(unit=5, nml=qg_time)
  read(unit=5, nml=denkf_true)
  read(unit=5, nml=denkf_obs)

  n = imax * jmax
  allocate(sbuf(n), y(nobs), obsoff(nstep/nsave+1))
  rstd = sqrt(rsig)
  dobs = n / nobs

  k = 1
  do i=tsave, nstep, nsave
    write(unit=pfile, fmt="('p'i0.6'.dat')") i
    open(unit=un, file=trim(tdir)//"/"//pfile, access="direct", &
      status="old", action="read", recl=size(sbuf)*4)
    read(unit=un, rec=1) sbuf
    close(unit=un)
    obsoff(k) = math_random_integer(1, imax)
    iobs = obsoff(k)
    do j=1, nobs
      y(j) = sbuf(iobs) + real(rstd * math_random_normal(), kind=sp)
      iobs = iobs + dobs
    end do
    write(unit=ofile, fmt="('o'i0.6'.dat')") i
    open(unit=un, file=trim(odir)//"/"//ofile, access="direct", &
      status="replace", action="write", recl=size(y)*4)
    write(unit=un, rec=1) y
    close(unit=un)
    k = k + 1
  end do
  open(unit=un, file=trim(odir)//"/"//"obsoff.txt", &
    status="replace", action="write")
  write(unit=un, fmt="(i4)") obsoff
  close(unit=un)

end program genobs
