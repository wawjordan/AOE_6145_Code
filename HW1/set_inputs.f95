!============================== set_inputs ================================80
module set_inputs

  use set_precision, only : prec
  use set_constants, only : zero, one, two, half, pi

  implicit none

  private

  public :: iq, iSS, max_newton_iter, newton_tol
  public :: p0, T0, Astar

  integer :: max_newton_iter = 20
  integer :: iq = 5
  integer :: iSS = 1
  real(prec) :: newton_tol = 1.0e-15_prec
  real(prec) :: p0 = 300.0_prec
  real(prec) :: T0 = 600.0_prec
  real(prec) :: Astar = 0.2_prec
  real(prec) :: g = 1.4_prec
  real(prec) :: gp1 = zero
  real(prec) :: gm1 = zero
  real(prec) :: Ru = 8314.0_prec
  real(prec) :: Mair = 28.96_prec
  real(prec) :: R = zero
  real(prec) :: a0 = zero
  real(prec) :: rho0 = zero
  real(prec) :: xmin = -one
  real(prec) :: xmax = one
  real(prec), dimension(:), allocatable :: xq
  real(prec), dimension(:), allocatable :: A

  contains
  function area(xq)
    real(prec) :: area
    real(prec), intent(in) :: xq
    area = 0.2_prec + 0.4_prec*( one + sin(pi*(xq-0.5_prec)))
  end function area

  subroutine set_derived_inputs

    implicit none
    integer :: i, i1, i2
    reak(prec), external :: area
    i = 1
    i1 = 1
    i2 = iq
    allocate(xq(i1:iq))
    allocate(A(i1:iq))
    do i = i1,iq
      xq(i) = xmin + i*(xmax-xmin)/iq
    end do
    do i = i1,iq
      A(i) = area(xq(i))
    end do
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
