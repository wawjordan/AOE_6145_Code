module exact_q1d_type

  use set_precision, only : prec
  use set_constants, only : half, one, two
  use set_inputs, g, gp1, gm1, A, Astar, iSS

  implicit none

  type exact_q1d_t

    real(prec), allocatable, dimension(:) :: M
    real(prec), allocatable, dimension(:) :: T
    real(prec), allocatable, dimension(:) :: rho
    real(prec), allocatable, dimension(:) :: u
    real(prec), allocatable, dimension(:) :: p

  end type exact_q1d_t

contains

  subroutine allocate_exact_q1d( soln )
    use set constants, only : zero, one
    use set_inputs, only : iq
    implicit none
    type(exact_q1d_t), intent(inout) :: soln
    integer :: i1, i2

    i1 = 1
    i2 = iq

    allocate( soln%M(i1:i2), soln%T(i1:i2), soln%rho(i1,i2) &
              soln%u(i1,i2), soln%p(i1:i2) )
    soln%M = zero
    soln%T = zero
    soln%rho = zero
    soln%u = zero
    soln%p = zero
  end subroutine allocate_exact_q1d

  subroutine deallocate_exact_q1d( soln )
    implicit none
    type(exact_q1d_t), intent(inout) :: soln
    deallocate( soln%M, soln%T, soln%rho, soln%u, soln%p )
  end subroutine deallocate_exact_q1d

  function f(M, A)
    real(prec) :: f
    real(prec), intent(in) :: M
    real(prec), intent(in) :: A
    f = ( (two/gp1)*(one + half*gm1*M**2) )**(gp1/gm1) - (A/Astar)**2
  end function f

  function df(M, A)
    real(prec) :: df
    real(prec), intent(in) :: M
    real(prec), intent(in) :: A
    df = two*M*( ( (two/gp1)*(one + half*gm1*M**2) )**(two/gm1) - (A/Astar)**2 )
  end subroutine df


end module exact_q1d_type
