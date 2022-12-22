module pbepack_quadratures
   use pbepack_kinds
   use hrweno_grids, only: grid1
   use stdlib_optval, only: optval
   implicit none
   private

   public :: quadgrid1, evalmoment

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
   !! Cell-integral or cell-average of \( f(x) \) over grid, using Simpson's 1/3 rule.
   !! ```
   !!   f(x) ^
   !!        |                * * * *
   !!        |          * * *++++++++++ *
   !!        |        *                  *
   !!        |     ++*++++++               *
   !!        |     *                    ++++*+++++
   !!        |                                 * *
   !!       -|----|----.----|-----.----|----.-----|---->
   !!               x_{i-i}      x_i      x_{i+1}
   !! ```
      procedure(fx) :: fnc
         !! function \( f(x) \) to integrate/average over grid cells
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

   pure real(rk) function evalmoment(u, grid, order, normalize) result(res)
   !! \(m\)-th order moment of \( u(x,t) \):
   !! $$
   !! M_m(t) = \int _0^\infty x^m u(x,t)dx = {\sum_i \bar{u_i}(t) \Delta x_i x_i^m}
   !! $$
      real(rk), intent(in) :: u(:)
         !! cell-average number density, \( \bar{u_i} \)
      type(grid1), intent(in) :: grid
         !! `grid1` object
      integer, intent(in) :: order
         !! order of the moment
      logical, intent(in), optional :: normalize
         !! if `true`, the result will be normalized by the 0-th moment

      ! Check `order`
      if (order < 0) then
         error stop "Invalid 'order'. Valid range: order >=0."
      end if

      ! Compute moment
      res = sum(u*grid%width*grid%center**order)

      ! Normalize if required
      if (optval(normalize, .false.)) then
         res = res/sum(u*grid%width)
      end if

   end function evalmoment

end module pbepack_quadratures
