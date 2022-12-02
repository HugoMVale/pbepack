module pbepack_math
!! This module implements some auxiliary math tools.
   use pbepack_kinds
   implicit none
   public

contains

   pure real(rk) function delta_kronecker(i, j) result(res)
   !! Delta kronecker \( \delta_{i,j} \).
      integer, intent(in) :: i, j
         !! integers

      res = ZERO
      if (i == j) res = ONE

   end function delta_kronecker

end module pbepack_math
