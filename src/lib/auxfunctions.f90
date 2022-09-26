module auxfunctions
!! Collection of auxiliary functions
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   public :: PI
   public :: delta_kronecker, heaviside, smooth_switch
   public :: gauss, lognorm, expo1d, expo2dmm, expo2dmw
   public :: average, stddev
   public :: islocalmax

   integer, parameter :: rk = real64
   real(rk), parameter :: PI = 4*atan(1._rk), ZERO = 0._rk, ONE = 1._rk

contains

   pure real(rk) function delta_kronecker(i, j)
    !! Delta kronecker /( \delta_{i,j} /).
      integer, intent(in) :: i
        !! Integer i
      integer, intent(in) :: j
        !! Integer j

      if (i == j) then
         delta_kronecker = ONE
      else
         delta_kronecker = ZERO
      end if

   end function delta_kronecker

   elemental real(rk) function heaviside(x)
      !! Heaviside function, \( H(x) \)
      real(rk), intent(in) :: x
        !! Independent variable

      if (x > ZERO) then
         heaviside = ONE
      else
         heaviside = ZERO
      end if

   end function heaviside

   elemental real(rk) function smooth_switch(x, xc, a, f1, f2)
      !! This function transforms the discontinuous function f(x) = If(x<xc,f1,f2) into a
      !! continuous function around xc.
      real(rk), intent(in) :: x
        !! independent variable
      real(rk), intent(in) :: xc
        !! critical value of x
      real(rk), intent(in) :: a
        !! parameter that controls the transition rate
      real(rk), intent(in) :: f1
        !! lower limit of f(x)
      real(rk), intent(in) :: f2
        !! upper limit of f(x)

      smooth_switch = f1 + 0.5_rk*(tanh((x - xc)/a) + ONE)*(f2 - f1)

   end function smooth_switch

   elemental real(rk) function gauss(x, mean, sigma)
      !! Gauss distribution (normal distribution).
      real(rk), intent(in) :: x
        !! ramdom variable
      real(rk), intent(in) :: mean
        !! mean value
      real(rk), intent(in) :: sigma
        !! standard deviation

      gauss = exp(-(x - mean)**2/(2*sigma**2))/(sigma*sqrt(2*Pi))

   end function gauss

   elemental real(rk) function lognorm(x, mean, sigma)
      !! Log-normal distribution.
      real(rk), intent(in) :: x
        !! ramdom variable
      real(rk), intent(in) :: mean
        !! mean value
      real(rk), intent(in) :: sigma
        !! standard deviation

      real(rk) :: meanLN, sigmaLN

      sigmaLN = sqrt(log(ONE + (sigma/mean)**2))
      meanLN = log(mean) - 0.5_rk*sigmaLN**2

      lognorm = ONE/(x*sqrt(2*Pi*sigmaLN**2)) &
                *exp(-(log(x) - meanLN)**2/(2*sigmaLN**2))

   end function lognorm

   elemental real(rk) function expo1d(x, x0, N0)
    !! 1D Exponential distribution.
      real(rk), intent(in) :: x
        !! random variable
      real(rk), intent(in) :: x0
        !! mean value
      real(rk), intent(in) :: N0
        !! integral of distribution

      expo1d = N0/(x0)*exp(-x/x0)

   end function expo1d

   elemental real(rk) function expo2dmm(m1, m2, m10, m20, N0)
      !! 2D Exponential distribution, with mass mass of each component as random variables.
      real(rk), intent(in) :: m1
        !! mass of component 1
      real(rk), intent(in) :: m2
        !! mass of component 2
      real(rk), intent(in) :: m10
        !! mean mass of component 1
      real(rk), intent(in) :: m20
        !! mean mass of component 2
      real(rk), intent(in) :: N0
        !! integral of distribution

      expo2dmm = N0/(m10*m20)*exp(-m1/m10 - m2/m20)

   end function expo2dmm

   elemental real(rk) function expo2dmw(m, w1, m10, m20, N0)
      !! 2D Exponential distribution, with total mass and composition of component 1 as
      !! independent variables.
      real(rk), intent(in) :: m
        !! total mass
      real(rk), intent(in) :: w1
        !! mass fraction of component 1
      real(rk), intent(in) :: m10
        !! mean mass of component 1
      real(rk), intent(in) :: m20
        !! mean mass of component 2
      real(rk), intent(in) :: N0
        !! integral of distribution

      real(rk) :: m1, m2

      m1 = m*w1
      m2 = m - m1

      expo2dmw = m*expo2dmm(m1, m2, m10, m20, N0)

   end function expo2dmw

   pure real(rk) function average(x, weight)
      !! This function computes the (weighted) average value of a vector x.
      real(rk), intent(in) :: x(:)
        !! vector or samples
      real(rk), intent(in), optional :: weight(:)
        !! optional vector of weights

      if (.not. present(weight)) then
         average = sum(x)/max(ONE, real(size(x), rk))
      else
         average = sum(x*weight)/sum(weight)
      end if

   end function average

   pure real(rk) function stddev(x)
      !! Standard deviation of a sample.
      real(rk), intent(in) :: x(:)
        !! vector of samples
      real(rk) :: xmean
      integer :: N

      N = size(x)

      if (N < 2) then
         stddev = ZERO
         return
      else
         xmean = average(x)
         stddev = sqrt(sum((x - xmean)**2)/(real(N, rk) - ONE))
      end if

   end function stddev

   pure logical function islocalmax(i, x)
      !! This function verifies if the ith element of vector x is a local maximum.
      integer, intent(in) :: i
        !! index of location to be checked
      real(rk), intent(in) :: x(:)
        !! vector
      real(rk), parameter :: tol = 1e-2_rk
      logical :: cond(2)

      if ((i < 3) .OR. (size(x) - i < 2)) then
         islocalmax = .false.
         return
      end if

      cond(1) = ((((x(i - 2) < x(i - 1)) .and. x(i - 1) < x(i)) &
                  .and. x(i) > x(i + 1)) .and. x(i + 1) > x(i + 2))
      cond(2) = x(i)/maxval(x) > tol

      islocalmax = all(cond)

   end function islocalmax

end module auxfunctions
