module pbepack_kinds
!! Aux module to define real precision.
   use, intrinsic :: iso_fortran_env, only: real32, real64
   implicit none
   private

   public :: rk, dp, ONE, ZERO, TWO, HALF

! #ifdef REAL32
!    integer, parameter :: rk = real32
! #elif REAL64
   integer, parameter :: rk = real64
! #else
!    integer, parameter :: rk = real64
! #endif
   integer, parameter :: dp = real64
   real(rk), parameter :: ZERO = 0._rk, ONE = 1._rk, TWO = 2._rk
   real(rk), parameter :: HALF = ONE/2

end module pbepack_kinds
