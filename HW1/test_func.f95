program test_func

  use set_precision, only : prec
  use set_constants, only : zero
  use set_inputs
  use subroutines, only : newton_safe, f, df
  use exact_q1d_type

  implicit none

  integer :: i

  call set_derived_constants
  call set_derived_inputs



end program test_func
