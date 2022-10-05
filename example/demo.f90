program demo
   use, intrinsic :: iso_fortran_env, only: real64
   use pbetools
   use grid, only: grid1
   use auxfunctions, only: expo1d
   implicit none

   integer, parameter :: rk = real64
   integer, parameter :: nc = 100
   type(pbe) :: mypbe
   type(grid1) :: gx

   ! Define the spatial grid
   ! In this example, we use a linear grid, but any smooth grid can be used
   call gx%linear(xmin=0._rk, xmax=1e2_rk, ncells=nc, name="my grid")

   call mypbe%init(grid=gx, agg=afun) !, break = , daughter = , nuc = )
   print *, mypbe%ncells
   print *, mypbe%grid%scl

   ! call mypbe%integrate(ic=, t0, tout, sol)

   ! call mypbe%eval(u, t, y)
contains

   elemental real(rk) function ic1d(self, x)
      !! Initial condition
      class(pbe), intent(in) :: self
      real(rk), intent(in) :: x
      ic1d = expo1d(x, x0=1._rk, n0=1._rk)
   end function

   pure real(rk) function afun(xa, xb, t, y)
   !! Particle aggregation frequency
      real(rk), intent(in) :: xa
         !! internal coordinate of particle a
      real(rk), intent(in) :: xb
         !! internal coordinate of particle b
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), parameter :: a1 = 1._rk
      afun = a1*xa*xb
   end function

end program demo

