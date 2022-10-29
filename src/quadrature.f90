module quadrature
   use real_kinds, only: rk
   use grids, only: grid1
   implicit none
   private

   abstract interface
      pure real(rk) function fx(x)
      !! Interface initial condition of PBE
         import rk
         real(rk), intent(in) :: x
            !! internal coordinate
      end function
   end interface

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

end module quadrature
