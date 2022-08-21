program genloc
  use, intrinsic :: iso_fortran_env, only: &
    sp => real32, dp => real64
  implicit none
  
  integer, parameter :: un = 51

  integer :: i1, i2, j1, j2, imax, jmax
  real(kind=sp) :: r0, r0r0
  character(len=*), parameter :: rfile = "loc.dat"

  real(kind=sp), dimension(:, :, :, :), allocatable :: rho

  namelist /qg_grid/ imax, jmax
  namelist /loc_param/ r0

  read(unit=5, nml=qg_grid)
  read(unit=5, nml=loc_param)
  print *, "r0=", r0
  allocate(rho(imax, jmax, imax, jmax))

  if (r0 > 0) then
    r0r0 = r0 * r0
    do j2=1, jmax
      do i2=1, imax
        do j1=1, jmax
          do i1=1, imax
            rho(i1, j1, i2, j2) = exp(-0.5e0 * ((i1 - i2)**2 + (j1 - j2)**2) / r0r0)
          end do
        end do
      end do
    end do
  end if

  print *, "rho: min=", minval(rho), " max=", maxval(rho)
  open(unit=un, file=rfile, access="direct", &
    status="replace", action="write", recl=size(rho)*4)
    write(unit=un, rec=1) rho
  close(unit=un)

end program genloc
