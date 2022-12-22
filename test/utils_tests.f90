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

   pure real(rk) function duniform(xd, xo, y, moment) result(res)
   !! Uniform daughter distribution kernel for 1D system
      real(rk), intent(in) :: xd
         !! internal coordinate of daughter particle
      real(rk), intent(in) :: xo
         !! internal coordinate of original particle
      real(rk), intent(in) :: y(:)
         !! environment vector
      integer, intent(in), optional :: moment
         !! moment of \( x \) conserved during breakage
      integer :: m
      m = optval(moment, 1)
      res = (TWO*m/xo**m)*xd**(m - 1)
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
      real(rk), intent(in), optional :: x0
         !! mean value
      real(rk), intent(in), optional :: n0
         !! initial number of particles
      real(rk) :: x0_, n0_
      x0_ = optval(x0, ONE)
      n0_ = optval(n0, ONE)
      res = (n0_/x0_)*exp(-x/x0_)
   end function expo1d

   pure function aconst_moments(t, x0, n0, a0) result(res)
   !! Moment solution for IC='expo1d' and constant aggregation frequency.
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in), optional :: x0
         !! mean value
      real(rk), intent(in), optional :: n0
         !! initial number of particles
      real(rk), intent(in), optional :: a0
         !! aggregation frequency
      real(rk) :: res(0:1)

      real(rk) :: x0_, n0_, a0_

      x0_ = optval(x0, ONE)
      n0_ = optval(n0, ONE)
      a0_ = optval(a0, ONE)
      res(0) = 2*n0_/(TWO + a0_*n0_*t)
      res(1) = n0_*x0_

   end function

   pure function asum_moments(t, x0, n0, a1) result(res)
   !! Moment solution for IC='expo1d' and sum aggregation frequency.
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in), optional :: x0
         !! mean value
      real(rk), intent(in), optional :: n0
         !! initial number of particles
      real(rk), intent(in), optional :: a1
         !! aggregation frequency factor
      real(rk) :: res(0:1)

      real(rk) :: x0_, n0_, a1_

      x0_ = optval(x0, ONE)
      n0_ = optval(n0, ONE)
      a1_ = optval(a1, ONE)

      res(0) = n0_*exp(-a1_*n0_*x0_*t)
      res(1) = n0_*x0_

   end function

end module utils_tests
