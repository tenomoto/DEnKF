program run_denkf
  use, intrinsic :: iso_fortran_env, only: &
    sp => real32, dp => real64
  use denkf_module, only: denkf_init, denkf_clean, denkf_analysis
  implicit none
  
  integer, parameter :: &
    imax = 129, jmax = 129, n = imax * jmax, &
    nens = 25, nobs = 300, un = 51
  char(len=*), parameter :: &
    ffile = "forecast.dat", afile = "analysis.dat"

  real(kind=sp), dimension(:), allocatable :: x
  real(kind=sp), dimension(:, :), allocatable :: p12

  integer :: i
  char(len=8) :: efile

  allocate(x(n), p12(n, nens), )
  open(unit=un, file=ffile, access="direct", &
    status="old", action="read", recl=n*4)
  read(unit=un, rec=1) x
  close(unit=un)
  do i=1, nens
    open(unit=un, file=efile, access="direct", &
      status="old", action="read", recl=n*4)
    read(unit=un, rec=1) p12(:, i)
    close(unit=un)
    p12(:, i) = (p12(:, i) - x(:)) / sqrt(nens - 1)
  end do
  call denkf_init(n, nobs)
  call denkf_init(n, nobs)
  call denkf_clean()

end program run_denkf
