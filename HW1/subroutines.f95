module subroutines

  implicit none

  private

  public newton_safe, newton_safe2, newton_safe3

  contains

  !============================== newton_safe ================================80
  !>
  !! Description: Safe-guarded 1D Newton's method.
  !!
  !! Inputs:      f:    Single variable function f(x).
  !!              df:   Derivative of f(x).
  !!              a:    Beginning of bracketing interval.
  !!              b:    End of bracketing interval.
  !!
  !! Outputs:     x:    Zero of f(x).
  !!              xk:   Array of iterates.
  !!              e:    Array of approximate relative iterative errors.
  !<
  !===========================================================================80
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
    ! xk(2) = b2
    e(1) = abs(xk(2)-xk(1))/abs(xk(2))
    do k = 2, max_newton_iter
      if (e(k-1) < newton_tol) then
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
    x = xk(k)

  end subroutine newton_safe


  !============================== newton_safe2 ===============================80
  !>
  !! Description: Safe-guarded 1D Newton's method, but with hard-coded input
  !!              parameter 'A1' for function f(x) because I'm dumb.
  !!
  !! Inputs:      A1:   Area.
  !!              f:    two variable function f(x,A1).
  !!              df:   Derivative of f(x) w.r.t. x.
  !!              a:    Beginning of bracketing interval.
  !!              b:    End of bracketing interval.
  !!
  !! Outputs:     x:    Zero of f(x).
  !!              xk:   Array of iterates.
  !!              e:    Array of approximate relative iterative errors.
  !<
  !===========================================================================80
  subroutine newton_safe2( A1, f1, df1, a, b, x, xk, e)

    use set_precision, only : prec
    use set_constants, only : half, one, zero
    use set_inputs

    implicit none

    integer :: k = 1
    real(prec), external :: f1, df1
    real(prec), intent(in) :: a, b, A1
    real(prec)             :: a2, b2
    real(prec), intent(out) :: x
    real(prec), intent(out), dimension(max_newton_iter+1) :: xk
    real(prec), intent(out), dimension(max_newton_iter) :: e

    xk = zero

    a2 = a
    b2 = b
    xk(1) = a2
    xk(2) = xk(1) - f1(xk(1),A1)/df1(xk(1),A1)

    if ( (xk(2) < a2).or.(xk(2) > b2) ) then
      xk(2) = a2 + half*(b2-a2)
      if ( sign(one,f1(a2,A1)) == sign(one,f1(xk(2),A1)) ) then
        a2 = xk(2)
      else
        b2 = xk(2)
      endif
    endif

    e(1) = abs(xk(2)-xk(1))/abs(xk(2))
    do k = 2, max_newton_iter
      if (e(k-1) < newton_tol) then
        goto 200  !exit iteration loop if converged
      else
        continue
      endif
      xk(k+1) = xk(k) - f1(xk(k),A1)/df1(xk(k),A1)
      if ( (xk(k+1) < a2).or.(xk(k+1) > b2) ) then
        xk(k+1) = a2 + half*(b2-a2)
        if ( sign(one,f1(a2,A1)) == sign(one,f1(xk(k+1),A1)) ) then
          a2 = xk(k+1)
        else
          b2 = xk(k+1)
        endif
      endif
      e(k) = abs(xk(k)-xk(k-1))/abs(xk(k))
    end do

    200 continue
    x = xk(k)

  end subroutine newton_safe2

  !============================== newton_safe3 ===============================80
  !>
  !! Description: Safe-guarded 1D Newton's method, but with hard-coded input
  !!              parameter 'A1' and reduced memory overhead because I'm dumb
  !!              and hard-coded array sizes in a previous version of the code.
  !!
  !! Inputs:      A1:   Area.
  !!              f:    two variable function f(x,A1).
  !!              df:   Derivative of f(x) w.r.t. x.
  !!              a:    Beginning of bracketing interval.
  !!              b:    End of bracketing interval.
  !!
  !! Outputs:     x:    Zero of f(x).
  !!              xk:   Array of last two iterates.
  !!              e:    approximate relative iterative error.
  !<
  !===========================================================================80
  subroutine newton_safe3( A1, f1, df1, a, b, x, xk, e)

    use set_precision, only : prec
    use set_constants, only : half, one, zero
    use set_inputs

    implicit none

    integer :: k = 1
    real(prec), external :: f1, df1
    real(prec), intent(in) :: a, b, A1
    real(prec)             :: a2, b2
    real(prec), intent(out) :: x
    real(prec) :: x_new
    real(prec), intent(out), dimension(2) :: xk
    real(prec), intent(out) :: e

    xk = zero

    a2 = a
    b2 = b
    xk(1) = a2
    xk(2) = xk(1) - f1(xk(1),A1)/df1(xk(1),A1)

    if ( (xk(2) < a2).or.(xk(2) > b2) ) then
      xk(2) = a2 + half*(b2-a2)
      if ( sign(one,f1(a2,A1)) == sign(one,f1(xk(2),A1)) ) then
        a2 = xk(2)
      else
        b2 = xk(2)
      endif
    endif

    e = abs(xk(2)-xk(1))/abs(xk(2))
    do k = 2, max_newton_iter
      if (e < newton_tol) then
        goto 200  !exit iteration loop if converged
      else
        continue
      endif
      x_new = xk(2) - f1(xk(2),A1)/df1(xk(2),A1)
      if ( (x_new < a2).or.(x_new > b2) ) then
        x_new = a2 + half*(b2-a2)
        if ( sign(one,f1(a2,A1)) == sign(one,f1(x_new,A1)) ) then
          a2 = x_new
        else
          b2 = x_new
        endif
      endif
      xk(1) = xk(2)
      xk(2) = x_new
      e = abs(xk(2)-xk(1))/abs(xk(2))
    end do

    200 continue
    x = xk(2)

  end subroutine newton_safe3

end module subroutines
!=============================================================================80
