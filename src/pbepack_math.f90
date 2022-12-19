module pbepack_math
!! Auxiliary math functions.
   use pbepack_kinds
   implicit none
   public

contains

   elemental real(rk) function delta_kronecker(i, j) result(res)
   !! Delta kronecker \( \delta_{i,j} \).
      integer, intent(in) :: i, j
         !! integer

      res = ZERO
      if (i == j) res = ONE

   end function delta_kronecker

   elemental real(rk) function boxcar(x, a, b, height) result(res)
   !! Boxcar function.
      real(rk), intent(in) :: x
         !! argument
      real(rk), intent(in) :: a, b
         !! interval limit \( a < b \)
      real(rk), intent(in) :: height
         !! pulse height

      res = ZERO
      if (x > a .and. x < b) res = height

   end function boxcar

end module pbepack_math
