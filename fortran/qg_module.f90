module qg_module
  use, intrinsic :: iso_fortran_env, only: &
    sp => real32, dp => real64
  use math_module, only: math_tau
  use fd_module, only: fd_laplacian, fd_jacobian
  use mg_module, only: mg_vcycle
  implicit none

  real(kind=dp), public  :: &
    qg_beta = 1.0d0, qg_f = 1.6d3, qg_eps = 1.0d-5, &
    qg_a = 2.0d-12, qg_tau0 = -math_tau, qg_d
  real(kind=dp), dimension(:, :), allocatable, public :: qg_psi

  real(kind=dp), dimension(:), allocatable, private :: y
  real(kind=dp), private :: ddr, d2r

  public :: qg_init, qg_clean, qg_psi2q, qg_step, qg_save, qg_load

contains

  subroutine qg_init(imax, jmax)
    integer, intent(in) :: imax, jmax

    integer :: j

    allocate(qg_psi(imax, jmax), y(jmax))
    do j=1, jmax
      y(j) = 1.0d0 * (j - 1) / (jmax - 1)
    end do
    qg_d = y(2) - y(1)
    d2r = 1.0d0 / (qg_d * 2)
    ddr = 1.0d0 / (qg_d * qg_d)
    qg_psi(:, :) = 0.0d0

  end subroutine qg_init

  subroutine qg_clean()

    deallocate(qg_psi, y)

  end subroutine qg_clean

  subroutine qg_psi2q(psi, q)
    real(kind=dp), dimension(:, :), intent(in) :: psi
    real(kind=dp), dimension(:, :), intent(inout) :: q

    integer :: i, j, imax, jmax

    imax = size(q, 1)
    jmax = size(q, 2)
    do j=2, jmax-1
      do i=2, imax-1
        q(i, j) = fd_laplacian(psi, i, j) - qg_f * psi(i, j)
      end do
    end do

  end subroutine qg_psi2q

  function qg_step(q, params) result(dqdt)
    real(kind=dp), dimension(:, :), intent(in) :: q
    real(kind=dp), dimension(:), intent(in) :: params

    real(kind=dp) :: psix, lap3psi, jac
    real(kind=dp), dimension(:, :), allocatable :: &
      dqdt, lap1psi, lap2psi
    integer :: i, j, imax, jmax

    imax = size(q, 1)
    jmax = size(q, 2)
    allocate(dqdt(imax, jmax), lap1psi(imax, jmax), lap2psi(imax, jmax))
    dqdt(:, :) = 0.0d0
    lap1psi(:, :) = 0.0d0
    lap2psi(:, :) = 0.0d0

    call mg_vcycle(qg_psi, q, qg_d, qg_f)
    do j=2, jmax-1
      do i=2, imax-1
        lap1psi(i, j) = q(i, j) + qg_f * qg_psi(i, j)
      end do
    end do
    do j=2, jmax-1
      do i=2, imax-1
        lap2psi(i, j) = fd_laplacian(lap1psi, i, j) * ddr
      end do
    end do
    do j=2, jmax-1
      do i=2, imax-1
        psix = (qg_psi(i+1, j) - qg_psi(i-1, j)) * d2r
        lap3psi = fd_laplacian(lap2psi, i, j) * ddr
        jac = fd_jacobian(qg_psi, q, i, j) * ddr
        dqdt(i, j) = -qg_beta * psix - qg_eps * jac - qg_a * lap3psi &
          + qg_tau0 * sin(math_tau * y(j))
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

  subroutine qg_load(un, vname, i, darray)
    integer, intent(in) :: un, i
    character(len=1), intent(in) :: vname
    real(kind=dp), dimension(:, :), allocatable, intent(inout) :: darray

    real(kind=sp), dimension(:, :), allocatable :: sarray

    integer :: imax, jmax
    character(len=11) :: fname

    imax = size(darray, 1)
    jmax = size(darray, 2)
    allocate(sarray(imax, jmax))
    write(unit=fname, fmt="(a1,i0.6,'.dat')") vname, i
    open(unit=un, file=fname, status="old", &
      action="read", access="direct", recl=size(sarray)*4)
    read(unit=un, rec=1) sarray
    close(unit=un)
    print *, minval(sarray), maxval(sarray)
    darray(:, :) = sarray(:, :)

  end subroutine qg_load

end module qg_module
