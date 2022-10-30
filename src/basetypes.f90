module basetypes
   use real_kinds, only: rk
   use grids, only: grid1
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
   !! Abstract 1D PBE term class (e.g., aggregation, growth, etc.).
      type(grid1), pointer :: grid => null()
         !! pointer to grid object
      logical :: inited = .false.
         !! flag initialization
   contains
      procedure, pass(self) :: set_grid
   end type pbeterm

   type, extends(pbeterm), abstract :: particleterm
   !! Abstract 1D PBE particle term class (aggregation, breakage).
      integer :: moment
         !! moment of 'x' to be conserved upon aggregation
      real(rk), allocatable :: source(:)
         !! source(+) term
      real(rk), allocatable :: sink(:)
         !! sink(-) term
   contains
      procedure, pass(self) :: set_moment
   end type particleterm

contains

   subroutine set_grid(self, grid)
   !! Setter method for grid.
      class(pbeterm), intent(inout) :: self
         !! object
      type(grid1), intent(in), target :: grid
         !! grid1 object

      if (grid%ncells > 0) then
         self%grid => grid
      else
         self%msg = "Invalid 'grid'."
         self%ierr = 1
         error stop self%msg
      end if

   end subroutine

   pure subroutine set_moment(self, moment)
   !! Setter method for moment.
      class(particleterm), intent(inout) :: self
         !! object
      integer, intent(in) :: moment
         !! moment of 'x' to be conserved upon aggregation

      if (moment > 0) then
         self%moment = moment
      else
         self%msg = "Invalid 'moment'. Valid range: moment >= 1."
         self%ierr = 1
         error stop self%msg
      end if

   end subroutine

end module basetypes
