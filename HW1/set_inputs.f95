!============================== set_inputs ================================80
module set_inputs

  use set_precision, only : prec
  use set_constants, only : zero, one, two, half

  implicit none

  private

  public :: iq, iSS, max_newton_iter, newton_tol
  public :: p0, T0, Astar

  integer :: max_newton_iter
  integer :: iq
  integer :: iSS
  real(prec) :: newton_tol
  real(prec) :: p0
  real(prec) :: T0
  real(prec) :: Astar
  real(prec) :: g = 1.4_prec
  real(prec) :: gp1 = zero
  real(prec) :: gm1 = zero
  real(prec) :: Ru = 8314.0_prec
  real(prec) :: Mair = 28.96_prec
  real(prec) :: R = zero
  real(prec) :: a0 = zero
  real(prec) :: rho0 = zero

  contains

  subroutine set_derived_inputs

    implicit none

    R = Ru/Mair
    gp1 = g + one
    gm1 = g - one
    a0 = sqrt(g*R*T0)
    rho0 = p0/(R*T0)

    write(*,*) 'R     = ', R, ' [J/(kmol*K)]'
    write(*,*) 'gamma = ', g
    write(*,*) 'a_0   = ', a0, ' [m/s]'
    write(*,*) 'rho_0 = ', rho0, ' [kg/m^3]'
    write(*,*) 'P_0   = ', p0, ' [kPa]'
    write(*,*) 'T_0   = ', T0, ' [K]'
    write(*,*) 'A*    = ', Astar, ' [m^2]'

  end subroutine set_derived_inputs



end module set_inputs
!comment
