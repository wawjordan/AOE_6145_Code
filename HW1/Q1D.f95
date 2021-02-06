!=================================== Q1D =====================================80
module Q1D

  use set_precision, only : prec
  use set_constants, only : zero, one, two, half
  use set_inputs,    only : g, gp1, gm1, R
  use subroutines,   only : newton_safe
  implicit none

  private

  public M, rho, u, p, T

  integer    :: k
  integer    :: N
  real(prec), dimension (:), allocatable :: M
  real(prec), dimension (:), allocatable :: psi
  real(prec), dimension (:), allocatable :: rho
  real(prec), dimension (:), allocatable :: u
  real(prec), dimension (:), allocatable :: p
  real(prec), dimension (:), allocatable :: T


end module Q1D
!=============================================================================80
