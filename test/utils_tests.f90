module utils_tests
!! Utilities for tests.
   use pbepack_kinds
   use stdlib_optval, only: optval
   implicit none

contains

   pure real(rk) function aconst(xa, xb, y) result(res)
      !! Constant aggregation kernel for 1D system
      real(rk), intent(in) :: xa
         !! internal coordinate of particle a
      real(rk), intent(in) :: xb
         !! internal coordinate of particle b
      real(rk), intent(in) :: y(:)
         !! environment vector
      res = ONE
   end function

   pure real(rk) function asum(xa, xb, y) result(res)
      !! Sum aggregation kernel for 1D system
      real(rk), intent(in) :: xa
         !! internal coordinate of particle a
      real(rk), intent(in) :: xb
         !! internal coordinate of particle b
      real(rk), intent(in) :: y(:)
         !! environment vector
      res = (xa + xb)
   end function

   pure real(rk) function aprod(xa, xb, y) result(res)
      !! Product aggregation kernel for 1D system
      real(rk), intent(in) :: xa
         !! internal coordinate of particle a
      real(rk), intent(in) :: xb
         !! internal coordinate of particle b
      real(rk), intent(in) :: y(:)
         !! environment vector
      res = xa*xb
   end function

   pure real(rk) function bconst(x, y) result(res)
   !! Constant breakage kernel for 1D system
      real(rk), intent(in) :: x
         !! internal coordinate of particle
      real(rk), intent(in) :: y(:)
         !! environment vector
      res = ONE
   end function

   pure real(rk) function bsquare(x, y) result(res)
   !! Square breakage kernel for 1D system
      real(rk), intent(in) :: x
         !! internal coordinate of particle
      real(rk), intent(in) :: y(:)
         !! environment vector
      res = x**2
   end function

   pure real(rk) function duniform(xd, xo, y, m) result(res)
   !! Uniform daughter distribution kernel for 1D system
      real(rk), intent(in) :: xd
         !! internal coordinate of daughter particle
      real(rk), intent(in) :: xo
         !! internal coordinate of original particle
      real(rk), intent(in) :: y(:)
         !! environment vector
      integer, intent(in), optional :: m
         !! moment of \( x \) conserved during breakage
      integer :: m_
      m_ = optval(m, 1)
      res = (TWO*m_/xo**m_)*xd**(m_ - 1)
   end function

   pure real(rk) function gconst(x, y) result(res)
   !! Constant growth rate for 1D system
      real(rk), intent(in) :: x
         !! internal coordinate
      real(rk), intent(in) :: y(:)
         !! environment vector
      res = ONE
   end function

   pure real(rk) function glinear(x, y) result(res)
   !! Linear growth rate for 1D system
      real(rk), intent(in) :: x
         !! internal coordinate
      real(rk), intent(in) :: y(:)
         !! environment vector
      res = x
   end function

   elemental real(rk) function expo1d(x, x0, n0) result(res)
   !! 1D Exponential distribution.
      real(rk), intent(in) :: x
         !! random variable
      real(rk), intent(in) :: x0
         !! mean value
      real(rk), intent(in) :: n0
         !! initial number of particles
      res = (n0/x0)*exp(-x/x0)
   end function expo1d

   elemental real(rk) function solution_case1(x, x0, n0, a0, t) result(res)
   !! Anylytical solution for IC='expo1d' and constant aggregation frequency.
      real(rk), intent(in) :: x
         !! random variable
      real(rk), intent(in) :: x0
         !! mean value
      real(rk), intent(in) :: n0
         !! initial number of particles
      real(rk), intent(in) :: a0
         !! aggregation frequency
      real(rk), intent(in) :: t
         !! time

      real(rk) :: tau

      tau = n0*a0*t
      res = 4*n0/(x0*(tau + 2._rk)**2)*exp(-2*(x/x0)/(tau + 2._rk))

   end function solution_case1

end module utils_tests
