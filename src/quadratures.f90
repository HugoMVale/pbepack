module quadratures
   use real_kinds, only: rk
   use grids, only: grid1
   use stdlib_optval, only: optval
   implicit none
   private

   public :: quad1_r1

   abstract interface
      pure real(rk) function fx(x)
      !! Interface integrand
         import rk
         real(rk), intent(in) :: x
            !! argument
      end function
   end interface

   ! interface quad1
   !    module procedure :: quad1_r1
   ! end interface

contains

   pure function cellavg1(fnc, grid) result(res)
      procedure(fx) :: fnc
         !! function f(x) to average over grid cells
      type(grid1), intent(in) :: grid
         !! grid object
      real(rk), dimension(grid%ncells) :: res

      real(rk) :: fedges(0:grid%ncells), fcenter(grid%ncells)
      integer :: i, nc

      nc = grid%ncells

      do concurrent(i=0:nc)
         fedges(i) = fnc(grid%edges(i))
      end do

      do concurrent(i=1:nc)
         fcenter(i) = fnc(grid%center(i))
      end do

      res = (4*fcenter + fedges(0:nc - 1) + fedges(1:nc))/6

   end function cellavg1

   pure function quad1(fnc, xedges, order, average) result(res)
   !! Quadrature of \( f(x) \) over user-supplied grid.
      procedure(fx) :: fnc
         !! function f(x) to integrate
      real(rk), intent(in) :: xedges(0:)
         !! edges of the integration grid
      integer, intent(in) :: order
         !! order of the integration (order = 2 or 3)
      logical, intent(in), optional :: average
         !! flag to compute average of f(x) rather than integral
      real(rk) :: res(size(xedges) - 1)

      real(rk) :: fvalues(0:size(xedges) - 1)
      integer :: i, nc

      nc = size(xedges) - 1

      ! Evaluate f(x) at cell edges
      do concurrent(i=0:nc)
         fvalues(i) = fnc(xedges(i))
      end do

      ! Compute 2nd order approx (trapezium rule)
      res = (fvalues(0:nc - 1) + fvalues(1:nc))/2

      ! Compute 3rd order approx (Simpson's rule)
      select case (order)
      case (2)
      case (3)
         do concurrent(i=1:nc)
            fvalues(i) = fnc((xedges(i - 1) + xedges(i))/2)
         end do
         res = (res + 2*fvalues(1:nc))/3
      case default
         error stop "Invalid 'order'"
      end select

      ! Convert integral-average to integral
      if (.not. (optval(average, .false.))) then
         res = res*(xedges(1:nc) - xedges(0:nc - 1))
      end if

   end function quad1

   pure real(rk) function quad1_r1(fnc, a, b, order) result(res)
   !! Integral of \( f(x) \) over interval [a,b].
      procedure(fx) :: fnc
         !! function f(x) to integrate
      real(rk), intent(in) :: a
         !! lower limit of the integration interval
      real(rk), intent(in) :: b
         !! upper limit of the integration interval
      integer, intent(in) :: order
         !! order of the integration (order = 2 or 3)

      ! Compute 2nd order approx (trapezium rule)
      res = (fnc(a) + fnc(b))/2

      ! Compute 3rd order approx (Simpson's rule)
      if (order == 3) then
         res = (res + 2*fnc((a + b)/2))/3
      end if
      res = res*(b - a)

   end function quad1_r1

end module quadratures
