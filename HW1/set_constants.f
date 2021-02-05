!============================== set_constants ================================80
module set_constants

  use set_precision, only : Prec

  implicit none

  private

  public :: zero, one, two, three, four, six
  public :: half, third, fourth, fifth, sixth, tenth

  real(Prec), parameter :: zero   = 0.0_Prec
  real(Prec), parameter :: tenth  = 0.1_Prec
  real(Prec), parameter :: sixth  = 1.0_Prec/6.0_Prec
  real(Prec), parameter :: fifth  = 0.2_Prec
  real(Prec), parameter :: fourth = 0.25_Prec
  real(Prec), parameter :: third  = 1.0_Prec/3.0_Prec
  real(Prec), parameter :: half   = 0.5_Prec
  real(Prec), parameter :: one    = 1.0_Prec
  real(Prec), parameter :: two    = 2.0_Prec
  real(Prec), parameter :: three  = 3.0_Prec
  real(Prec), parameter :: four   = 4.0_Prec
  real(Prec), parameter :: six    = 6.0_Prec

end module set_constants
!comment
