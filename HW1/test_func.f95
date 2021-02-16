program test_func

  use set_precision, only : prec
  use set_constants, only : zero, set_derived_constants
  use set_inputs
  use subroutines, only : newton_safe
  use exact_q1d_type

  implicit none

  integer :: i
  character(len=100) :: header_str
  type(exact_q1d_t) :: soln

  call set_derived_constants
  ! call read_in
  call set_derived_inputs
  ! x, A, M, rho, u, p
  call allocate_exact_q1d(soln)

  call solve_exact_q1d(soln)
  write(header_str,*) '|    x   |    A    |         M         |'// &
  &  '        rho        |         u        |         p        |'
  100 format(2(F9.4),4(F20.14))
  write(*,*) trim(adjustl(header_str))
  do i = 1,iq
    write(*,100) xq(i), Aq(i), soln%M(i), soln%rho(i), soln%u(i), soln%p(i)
  end do
  open(30,file='out.dat',status='unknown')
  write(30,*) 'TITLE = "Quasi-1D Nozzle Flow Output"'
  write(30,*) 'variables="x(m)""A(m^2)""M""rho(kg/m^3)""u(m/s)""p(kPa)"'
  write(30,*) trim(adjustl(header_str))
  do i = 1,iq
    write(30,100) xq(i), Aq(i), soln%M(i), soln%rho(i), soln%u(i), soln%p(i)
  end do
  close(30)
  call deallocate_exact_q1d(soln)

end program test_func
