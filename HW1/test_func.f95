program test_func

  use set_precision, only : prec
  use set_constants, only : zero, set_derived_constants
  use set_inputs
  use subroutines, only : newton_safe
  use exact_q1d_type

  implicit none

  integer :: i
  type(exact_q1d_t) :: soln

  call set_derived_constants
  call set_derived_inputs
  ! x, A, M, rho, u, p
  call allocate_exact_q1d(soln)

  call solve_exact_q1d(soln)
  write(*,*) '|    x   |    A    |         M         |        rho        |         u        |         p        |'
  do i = 1,iq
    write(*,'(F9.4,F9.4,F20.14,F20.14,F20.14,F20.14)') xq(i), Aq(i), soln%M(i), soln%rho(i), soln%u(i), soln%p(i)
  end do
  call deallocate_exact_q1d(soln)

end program test_func
