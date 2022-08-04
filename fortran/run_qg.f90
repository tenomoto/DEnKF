program run_qg
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use mg_module, only: mg_itermax, mg_tol
  use qg_module, only: qg_init, qg_step, qg_save, &
    qg_c, qg_f, qg_eps, qg_a, qg_tau0, qg_psi
  use ode_module, only: ode_rk4
  implicit none

  integer, parameter :: &
    n = 129, nstep = 1000, nsave = 100, un = 51
  real(kind=dp), parameter ::  dt = 1.5d0
  integer, dimension(3), parameter :: itermax = [1, 1, 5000]

  real(kind=dp), dimension(:), allocatable :: y
  real(kind=dp), dimension(:, :), allocatable :: q, dq
  real(kind=dp), dimension(1) :: dummy_params = [0.0d0]
  character(len=10) :: fname
  integer :: i, j, r

!  mg_itermax = [1, 1, 1000]
!  mg_tol = 1.0d-4
!  qg_c = 0.0d0
!  qg_f = 0.0d0
!  qg_eps = 1.0d0
!  qg_a = 2.0d-12
!  qg_tau0 = 0.0d0
  mg_itermax = [1, 1, 100]
  call qg_init(n)
  allocate(q(n, n), dq(n, n))
  q(:, :) = 0.0d0 
!  call random_number(q(2:n-1, 2:n-1))
!  q(2:n-1, 2:n-1) = 2.0d0 * q(2:n-1, 2:n-1) - 1.0d0
  call qg_save(un, "q", 0, q)
  call qg_save(un, "p", 0, qg_psi)

  r = 1
  do i=0, nstep
    print "(a, i6, 4(a, e12.4))", "step ", i, " p: min=", minval(qg_psi), " max=", maxval(qg_psi), &
      " q: min=", minval(q), " max=", maxval(q)
    if (mod(i, nsave) == 0) then
      call qg_save(un, "q", i, q)
      call qg_save(un, "p", i, qg_psi)
      r = r + 1
    end if
    dq(:, :) = ode_rk4(qg_step, q, dt, dummy_params)
    q(2:n-1, 2:n-1) = q(2:n-1, 2:n-1) + dq(2:n-1, 2:n-1)
  end do

  deallocate(q, dq)

end program run_qg
