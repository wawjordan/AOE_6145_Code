!================================ variables ==================================80
module variables

  use set_precision, only : prec
  use set_constants, only : zero, one, two, half

  implicit none

  private

  public :: xq, Aq, Mq, Tq, uq

  real(prec), dimension(:), allocatable :: xq
  real(prec), dimension(:), allocatable :: Aq
  real(prec), dimension(:), allocatable :: Mq
  real(prec), dimension(:), allocatable :: Tq
  real(prec), dimension(:), allocatable :: rhoq
  real(prec), dimension(:), allocatable :: uq
  real(prec), dimension(:), allocatable :: pq



end module variables
