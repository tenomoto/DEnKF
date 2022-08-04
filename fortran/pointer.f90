program pointer
  implicit none

  real, dimension(:), allocatable, target :: x
  real, dimension(:, :), pointer :: y

  allocate(x(6))
  x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  print *, "x=", x
  y(1:2, 1:2) => x(2:5)
  y(1, 2) = 0.0
  print *, "y="
  print *, y(1, :)
  print *, y(2, :)
  call mysub(y)
  print *, "y="
  print *, y(1, :)
  print *, y(2, :)
  y(1:2, 1:3) => x(1:6)
  print *, "y="
  print *, y(1, :)
  print *, y(2, :)
  call mysub(y)
  print *, "y="
  print *, y(1, :)
  print *, y(2, :)
  call mysub(y)
  deallocate(x)

contains

  subroutine mysub(a)
    real, dimension(:, :), intent(inout) :: a

    a = 2 * a
  end subroutine mysub

end program pointer

