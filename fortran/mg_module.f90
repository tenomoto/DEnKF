module mg_module
  use :: fd_module, only: fd_laplacian
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  type vc_type
    real(kind=dp), dimension(:), allocatable :: array
    integer, dimension(:), allocatable :: head, tail, n
  end type

  logical, public :: mg_debug = .false.
  integer, dimension(3), public :: mg_itermax = [1, 1, 1]
  real(kind=dp), public :: mg_tol = 1.0d-5

  public mg_vcycle, mg_prolong_test, mg_vcycle_test

contains

  subroutine jacobi_step(p, q, res, d, f, niter, tol)
    real(kind=dp), dimension(:, :), intent(inout) :: p
    real(kind=dp), dimension(:, :), intent(in) :: q
    real(kind=dp), intent(out) :: res
    real(kind=dp), intent(in) :: d, f, tol
    integer, intent(inout) :: niter

    integer :: i, j, k, n, kmax
    real(kind=dp) :: dd, r
    real(kind=dp), dimension(:, :), allocatable :: pu

    kmax = niter
    n = size(p, 1)
    dd = d * d 
    allocate(pu(n, n))
    pu(:, :) = 0.0d0
    pu(2:n-1, 2:n-1) = p(2:n-1, 2:n-1)

    do k=1, kmax
      res = 0.0d0
      do j=2, n-1
        do i=2, n-1
          r = fd_laplacian(p, i, j) / dd - f * p(i, j) - q(i, j)
          pu(i, j) = p(i, j) + dd * r / (4.0d0 + dd * f)
          res =  res + r * r
        end do
      end do
      p(2:n-1, 2:n-1) = pu(2:n-1, 2:n-1)
      res = sqrt(res * dd)
      if (res < tol) then
        exit
        niter = k
      endif
    end do
    niter = k - 1
    deallocate(pu)

  end subroutine jacobi_step

  subroutine restrict(pin, pout)
    real(kind=dp), dimension(:, :), intent(in) :: pin
    real(kind=dp), dimension(:, :), intent(inout) :: pout

    pout(:, :) = pin(::2, ::2)

  end subroutine restrict

  subroutine prolong(pin, pout)
    real(kind=dp), dimension(:, :), intent(in) :: pin
    real(kind=dp), dimension(:, :), intent(inout) :: pout

    integer :: imax, jmax

    imax = size(pout, 1)
    jmax = size(pout, 2)
    pout(::2, ::2) = pin(:, :)
    pout(::2, 2:jmax-1:2) = 0.5 * (pout(::2, 1:jmax-2:2) + pout(::2, 3:jmax:2))
    pout(2:imax-1:2, ::2) = 0.5 * (pout(1:imax-2:2, ::2) + pout(3:imax:2, ::2))
    pout(2:imax-1:2, 2:jmax-1:2) = 0.5 * (pout(2:imax-1:2, 1:jmax-2:2) + pout(2:imax-1:2, 3:jmax:2))

  end subroutine prolong

  subroutine mg_prolong_test

    integer, parameter :: &
      n = 5, m = 2*(n - 1) + 1, u = 51
    real(kind=dp), dimension(n, n) :: p = 0.0d0
    real(kind=dp), dimension(m, m) :: q = 0.0d0
    integer :: i, j
    real(kind=dp) :: x, y, d

    d = 2.0d0 / (n - 1)
    open(unit=u, file="prolong.dat", access="direct", recl=size(q)*8, &
      action="write", status="replace")
    y = -1.0d0 + d
    do j=2, n-1
      x = -1.0d0 + d
      do i=2, n-1
        p(i, j) = exp(-(x**2 + y**2))
        x = x + d
      end do
      y = y + d 
    end do
    call prolong(p, q)
    write(unit=u, rec=1) q
    close(unit=u)

  end subroutine mg_prolong_test
        
  subroutine vc_init(p, n)
    type(vc_type), intent(inout) :: p
    integer, intent(in) :: n

    integer :: nn, nlev, pmax, i

    nn = n
    nlev = nint(log10(real(n, dp))/log10(2.0d0))
    allocate(p % head(nlev), p % tail(nlev), p % n(nlev))
    pmax = 0
    do i=1, nlev
      p % n(i) = nn
      p % head(i) = pmax + 1
      pmax = pmax + nn**2
      p % tail(i) = pmax
      nn = nn / 2 + 1
    end do
    allocate(p % array(pmax))

  end subroutine vc_init

  subroutine vc_get(p, l, pl)
    type(vc_type), target, intent(inout) :: p
    integer, intent(in) :: l
    real(kind=dp), dimension(:, :), pointer, intent(out) :: pl

    pl(1:p % n(l), 1:p % n(l)) => p % array(p % head(l):p % tail(l))

  end subroutine vc_get

  subroutine mg_vcycle(p0, q0, d, f)
    real(kind=dp), dimension(:, :), intent(inout) :: p0
    real(kind=dp), dimension(:, :), intent(in) :: q0
    real(kind=dp), intent(in) :: d, f

    integer :: n, nlev, niter, l
    real(kind=dp) :: res, h
    real(kind=dp), dimension(:, :), pointer :: pin, pout
    real(kind=dp), dimension(:, :), pointer :: qin, qout
    type(vc_type) :: p, q

    integer, parameter :: u = 51
    
    n = size(p0, 1)
    call vc_init(p, n)
    call vc_init(q, n)
    call vc_get(p, 1, pin)
    call vc_get(q, 1, qin)
    pin(:, :) = 0.0d0
    qin(:, :) = 0.0d0
    pin(2:n-1, 2:n-1) = p0(2:n-1, 2:n-1)
    qin(2:n-1, 2:n-1) = q0(2:n-1, 2:n-1)
    nlev = size(p % head)
    h = d
    niter = mg_itermax(1)
    call jacobi_step(pin, qin, res, d, f, niter, mg_tol)
    if (mg_debug) then
      print *, "pre: res=", res, " niter=", niter
    end if
    do l=1, nlev-1
      call vc_get(p, l, pin)
      call vc_get(q, l, qin)
      call vc_get(p, l + 1, pout)
      call vc_get(q, l + 1, qout)
      call restrict(pin, pout)
      call restrict(qin, qout)
      h = h * 2
      niter = mg_itermax(2)
      call jacobi_step(pout, qout, res, h, f, niter, mg_tol)
      if (mg_debug) then
        print *, "restrict: ", l, " res=", res, " niter=", niter
      end if
    end do
    do l=nlev, 2, -1
      call vc_get(p, l, pin)
      call vc_get(p, l - 1, pout)
      call vc_get(q, l - 1, qout)
      call prolong(pin, pout)
      h = h / 2
      niter = mg_itermax(3)
      call jacobi_step(pout, qout, res, h, f, niter, mg_tol)
      if (mg_debug) then
        print *, "prolong: ", " res=", res, " niter=", niter
      end if
    end do
    p0(2:n-1, 2:n-1) = pout(2:n-1, 2:n-1)

  end subroutine mg_vcycle

  subroutine mg_vcycle_test
    integer, parameter :: n = 129, u = 51

    real(kind=dp), parameter :: &
      d = 1.0d0 / (n - 1), f = 0.0d0

    integer :: i, j
    real(kind=dp) :: x, y
    real(kind=dp), dimension(n, n) :: ptrue, p, q

    mg_itermax = [1, 1, 10000]
    mg_tol = 1.0d-4
    y = 0.0d0
    do j=1, n
      x = 0.0d0
      do i=1, n
        ptrue(i, j) = (x**2 - x**4) * (y**4 - y**2)
        q(i, j) = -2.0d0 * ( &
          (1.0d0 - 6.0d0 * x**2) * y**2 * (1.0d0 - y**2) + &
          (1.0d0 - 6.0d0 * y**2) * x**2 * (1.0d0 - x**2))
        p(i, j) = 0.0d0
        x = x + d
      end do
      y = y + d
    end do
    print *, "ptrue min=", minval(ptrue), " max=", maxval(ptrue)
    print *, "q min=", minval(q), " max=", maxval(q)
    call mg_vcycle(p, q, d, f)

  end subroutine mg_vcycle_test

end module mg_module
