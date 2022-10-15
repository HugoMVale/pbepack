module basetypes
   use grid, only: grid1
   implicit none
   private

   public :: base, pbeterm

   type, abstract :: base
   !! Base abstract class.
      character(:), allocatable :: name
         !! object name
      character(:), allocatable :: msg
         !! error message
      integer :: ierr = 0
         !! error status
      logical :: verbose = .false.
         !! verbose flag
   end type base

   type, extends(base), abstract :: pbeterm
   !! Abstract 1D PBE term class (e.g., aggregation, growth, etc.)
      type(grid1), pointer :: grid => null()
        !! pointer to grid object
      logical :: inited = .false.
   end type pbeterm

end module basetypes
