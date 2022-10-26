module basetypes
   use real_kinds, only: rk
   use grid, only: grid1
   implicit none
   private

   public :: base, pbeterm, particleterm

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
         !! flag initialization
   end type pbeterm

   type, extends(pbeterm), abstract :: particleterm
   !! Abstract 1D PBE particle term class (aggregation, breakage)
      integer :: moment
         !! moment of 'x' to be conserved upon aggregation
      real(rk), allocatable :: source(:)
         !! source(+) term
      real(rk), allocatable :: sink(:)
         !! sink(-) term
   end type particleterm

end module basetypes
