module qg_module
  use, intrinsic :: iso_fortran_env, only: &
    sp => real32, dp => real64
  use :: fd_module, only: fd_laplacian, fd_jacobian
  use :: mg_module, only: mg_vcycle
  implicit none

  real(kind=dp), private, parameter :: tau = 6.283185307179586d0
  real(kind=dp), public  :: &
    qg_c = 1.0d0, qg_f = 1.6d3, qg_eps = 1.0d-5, &
    qg_a = 2.0d-12, qg_tau0 = -tau
  real(kind=dp), dimension(:, :), allocatable, public :: qg_psi

  real(kind=dp), dimension(:), allocatable, private :: y
  real(kind=dp), private :: d, ddr, d2r

  public :: qg_init, qg_clean, qg_step, qg_save

contains

  subroutine qg_init(n)
    integer, intent(in) :: n

    integer :: j

    allocate(qg_psi(n, n), y(n))
    do j=1, n
      y(j) = 1.0d0 * (j - 1) / (n - 1)
    end do
    d = y(2) - y(1)
    d2r = 1.0d0 / (d * 2)
    ddr = 1.0d0 / (d * d)
    qg_psi(:, :) = 0.0d0

  end subroutine qg_init

  subroutine qg_clean()

    deallocate(qg_psi, y)

  end subroutine qg_clean

  function qg_step(q, params) result(dqdt)
    real(kind=dp), dimension(:, :), intent(in) :: q
    real(kind=dp), dimension(:), intent(in) :: params

    real(kind=dp) :: psix, lap3psi, jac
    real(kind=dp), dimension(:, :), allocatable :: &
      dqdt, lap1psi, lap2psi
    integer :: n, i, j

    n = size(q, 1)
    allocate(dqdt(n, n), lap1psi(n, n), lap2psi(n, n))
    dqdt(:, :) = 0.0d0
    lap1psi(:, :) = 0.0d0
    lap2psi(:, :) = 0.0d0

    call mg_vcycle(qg_psi, q, d, qg_f)
    do j=2, n-1
      do i=2, n-1
        lap1psi(i, j) = q(i, j) + qg_f * qg_psi(i, j)
      end do
    end do
    do j=2, n-1
      do i=2, n-1
        lap2psi(i, j) = fd_laplacian(lap1psi, i, j) * ddr
      end do
    end do
    do j=2, n-1
      do i=2, n-1
        psix = (qg_psi(i+1, j) - qg_psi(i-1, j)) * d2r
        lap3psi = fd_laplacian(lap2psi, i, j) * ddr
        jac = fd_jacobian(qg_psi, q, i, j) * ddr
        dqdt(i, j) = -qg_c * psix - qg_eps * jac - qg_a * lap3psi &
          + qg_tau0 * sin(tau * y(j))
      end do
    end do

  end function qg_step

  subroutine qg_save(un, vname, i, array)
    integer, intent(in) :: un, i
    character(len=1), intent(in) :: vname
    real(kind=dp), dimension(:, :), intent(in) :: array

    character(len=11) :: fname

    write(unit=fname, fmt="(a1,i0.6,'.dat')") vname, i
    open(unit=un, file=fname, status="replace", &
      action="write", access="direct", recl=size(array)*4)
    write(unit=un, rec=1) real(array, kind=sp)
    close(unit=un)

  end subroutine qg_save

end module qg_module
