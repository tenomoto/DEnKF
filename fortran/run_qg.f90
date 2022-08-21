program run_qg
  use, intrinsic :: iso_fortran_env, only: &
    sp => real32, dp => real64
  use math_module
  use mg_module, only: mg_itermax, mg_tol
  use qg_module, only: qg_init, qg_step, qg_psi2q, qg_save, qg_load, &
    qg_beta, qg_f, qg_eps, qg_a, qg_tau0, qg_psi
  use ode_module, only: ode_rk4
  implicit none

  integer, parameter :: un = 51

  integer :: imax, jmax, nstep, tsave, nsave
  real(kind=dp) ::  dt
  real(kind=dp), dimension(:), allocatable :: y
  real(kind=dp), dimension(:, :), allocatable :: q, dq
  real(kind=dp), dimension(1) :: dummy_params = [0.0d0]
  character(len=10) :: fname
  character(len=4) :: init
  integer :: i, j, r

  namelist /mg/ mg_itermax, mg_tol
  namelist /qg_grid/ imax, jmax
  namelist /qg_time/ init, nstep, tsave, nsave, dt
  namelist /qg_param/ qg_beta, qg_f, qg_eps, qg_a, qg_tau0

  read(unit=5, nml=mg)
  read(unit=5, nml=qg_grid)
  read(unit=5, nml=qg_time)
  read(unit=5, nml=qg_param)
  
  call qg_init(imax, jmax)
  allocate(q(imax, jmax), dq(imax, jmax))
  q(:, :) = 0.0d0 
  if (init == "rand") then
    call random_number(q(2:imax-1, 2:jmax-1))
    q(2:imax-1, 2:jmax-1) = 2.0d0 * q(2:imax-1, 2:jmax-1) - 1.0d0
  else if (init == "file") then
    call qg_load(un, "q", 0, q)
    call qg_load(un, "p", 0, qg_psi)
!    call qg_psi2q(qg_psi, q)
  end if

  r = 1
  do i=0, nstep
    print "(a, i6, 4(a, e12.4))", "step ", i, " p: min=", minval(qg_psi), " max=", maxval(qg_psi), &
      " q: min=", minval(q), " max=", maxval(q)
    if (i >= tsave .and. mod(i, nsave) == 0) then
      call qg_save(un, "q", i, q)
      call qg_save(un, "p", i, qg_psi)
      r = r + 1
    end if
    dq(:, :) = ode_rk4(qg_step, q, dt, dummy_params)
    q(2:imax-1, 2:jmax-1) = q(2:imax-1, 2:jmax-1) + dq(2:imax-1, 2:jmax-1)
  end do

  deallocate(q, dq)

end program run_qg
