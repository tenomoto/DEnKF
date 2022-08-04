program foo
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer :: n = 129, m, l, i, ri, d

  type v_cycle_type
    real(kind=dp), dimension(:), allocatable :: data
    integer, dimension(:), allocatable :: head, tail
  end type

  type(v_cycle_type) :: p

  l = nint(log10(real(n, dp))/log10(2.0d0))
  allocate(p % head(l), p % tail(l))
  m = 0
  print *, 0, m
  do i=1, l
    p % head(i) = m + 1
    d = 2**(l - i + 1) + 1
    m = m + d ** 2
    p % tail(i) = m
  end do
  print *, p % head
  print *, p % tail
  print *, m
  allocate(p % data(m))
  do i=1, l
    d = 2**(l - i + 1) + 1
    print *, shape(reshape(p % data(p % head(i):p % tail(i)), [d, d]))
  end do
  deallocate(p % head, p % tail)
  deallocate(p % data)
end program foo
