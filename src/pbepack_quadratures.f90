module pbepack_quadratures
   use pbepack_kinds
   use hrweno_grids, only: grid1
   use stdlib_optval, only: optval
   implicit none
   private

   public :: quadgrid1

   abstract interface
      pure real(rk) function fx(x)
      !! Interface integrand
         import rk
         real(rk), intent(in) :: x
            !! argument
      end function
   end interface

contains

   pure function quadgrid1(fnc, grid, average) result(res)
   !! Cell integral/average of \( f(x) \) over grid, using Simpson's 1/3 rule.
      procedure(fx) :: fnc
         !! function f(x) to integrate/average over grid cells
      type(grid1), intent(in) :: grid
         !! `grid1` object
      logical, intent(in), optional :: average
         !! flag to compute cell-average instead of cell-integral
      real(rk) :: res(grid%ncells)

      real(rk) :: fedges(0:grid%ncells), fcenter(grid%ncells)
      integer :: i

      associate (nc => grid%ncells)

         ! Evaluate f(x) at grid edges
         do concurrent(i=0:nc)
            fedges(i) = fnc(grid%edges(i))
         end do

         ! Evaluate f(x) at grid centers
         do concurrent(i=1:nc)
            fcenter(i) = fnc(grid%center(i))
         end do

         ! Cell-average using Simpson's 1/3 rule
         res = (4*fcenter + fedges(0:nc - 1) + fedges(1:nc))/6

         ! Convert cell-average to cell-integral
         if (.not. (optval(average, .false.))) then
            res = res*grid%width
         end if

      end associate

   end function quadgrid1

end module pbepack_quadratures
