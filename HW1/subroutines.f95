!=============================== subroutines =================================80
module subroutines

  implicit none

  private

  public newton_safe

  contains

  subroutine newton_safe( f, df, a, b, x, xk, e)

    use set_precision, only : prec
    use set_constants, only : half, one
    use set_inputs

    implicit none

    integer :: k = 1
    real(prec), external :: f, df
    real(prec), intent(in) :: a, b
    real(prec)             :: a2, b2
    real(prec), intent(out) :: x
    real(prec), intent(out), dimension(max_newton_iter+1) :: xk
    real(prec), intent(out), dimension(max_newton_iter) :: e

    a2 = a
    b2 = b
    xk(1) = a2
    xk(2) = xk(1) - f(xk(1))/df(xk(1))
    e(1) = abs(xk(k)-xk(k-1))/abs(xk(k))
    do k = 2, max_newton_iter
      if ( (e(k-1) < newton_tol).and.(k <= max_newton_iter) ) then
        goto 200  !exit iteration loop if converged
      else
        continue
      endif
      xk(k+1) = xk(k) - f(xk(k))/df(xk(k))
      if ( (xk(k+1) < a2).or.(xk(k+1) > b2) ) then
        xk(k+1) = a2 + half*(b2-a2)
        if ( sign(one,f(a2)) == sign(one,f(xk(k+1))) ) then
          a2 = xk(k+1)
        else
          b2 = xk(k+1)
        endif
      endif
      e(k) = abs(xk(k)-xk(k-1))/abs(xk(k))
    end do

    200 continue
    x = xk(k+1)

  end subroutine newton_safe

end module subroutines
!=============================================================================80
