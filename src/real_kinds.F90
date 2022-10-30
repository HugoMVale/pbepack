module real_kinds
!! Aux module to define real precision.
   use, intrinsic :: iso_fortran_env, only: real32, real64, real128
   implicit none
   private

   public :: rk, ONE, ZERO, HALF

! #ifdef REAL32
!    integer, parameter :: rk = real32
! #elif REAL64
   integer, parameter :: rk = real64
! #elif REAL128
!    integer, parameter :: rk = real128
! #else
!    integer, parameter :: rk = real64
! #endif
   real(rk), parameter :: ZERO = 0._rk, ONE = 1._rk
   real(rk), parameter :: HALF = ONE/2

end module real_kinds
