module pbepack_basetypes
!!   This module implements a number of (abstract) derived data types which form the basis for
!! the other derived types in this package.
   use pbepack_kinds
   use hrweno_grids, only: grid1
   use stdlib_optval, only: optval
   implicit none
   private

   public :: base, pbeterm, particleterm

   type, abstract :: base
   !! Abstract `base` class.
      character(:), allocatable :: name
         !! object name
      character(:), allocatable :: msg
         !! error message
      integer :: ierr = 0
         !! error code
      logical :: inited = .false.
         !! initialization flag
      logical :: verbose = .false.
         !! verbose flag
   contains
      procedure, pass(self) :: set_name
      procedure, pass(self) :: error_msg
      procedure, pass(self) :: check_inited
   end type base

   type, extends(base), abstract :: pbeterm
   !! Abstract 1D PBE term class (e.g., aggregation, growth, etc.).
      type(grid1), pointer :: grid => null()
         !! pointer to grid object
      real(rk), allocatable :: udot(:)
         !! net rate of change, \( d\bar{u}/dt \)
   contains
      procedure, pass(self) :: set_grid
      procedure, pass(self) :: pbeterm_allocations
   end type pbeterm

   type, extends(pbeterm), abstract :: particleterm
   !! Abstract 1D PBE particle term class (aggregation and breakage).
      integer :: moment = 1
         !! moment of \(x\) to be preserved upon aggregation/breakage
      real(rk), allocatable :: udot_birth(:)
         !! rate of birth
      real(rk), allocatable :: udot_death(:)
         !! rate of death
   contains
      procedure, pass(self) :: set_moment
      procedure, pass(self) :: particleterm_allocations
   end type particleterm

contains

   pure subroutine error_msg(self, msg)
   !! Error method.
      class(base), intent(inout) :: self
         !! object
      character(*), intent(in) :: msg
         !! message

      self%msg = msg
      self%ierr = 1
      error stop "Object: "//self%name//". Message: "//self%msg

   end subroutine

   pure subroutine set_name(self, name, default)
   !! Setter method for name.
      class(base), intent(inout) :: self
         !! object
      character(*), intent(in), optional :: name
         !! name
      character(*), intent(in), optional :: default
         !! default name

      self%name = optval(name, optval(default, ""))

   end subroutine

   pure subroutine check_inited(self)
   !! Check initialization method.
      class(base), intent(inout) :: self
         !! object

      if (.not. self%inited) then
         call self%error_msg("Object not initialized.")
      end if

   end subroutine

   subroutine set_grid(self, grid)
   !! Setter method for grid.
      class(pbeterm), intent(inout) :: self
         !! object
      type(grid1), intent(in), target :: grid
         !! 'grid1' object

      if (grid%ncells > 1) then
         self%grid => grid
      else
         call self%error_msg("Invalid 'grid'.")
      end if

   end subroutine

   pure subroutine set_moment(self, moment)
   !! Setter method for moment.
      class(particleterm), intent(inout) :: self
         !! object
      integer, intent(in) :: moment
         !! moment of \( x \) to be preserved upon aggregation/breakage

      if (moment > 0) then
         self%moment = moment
      else
         call self%error_msg("Invalid 'moment'. Valid range: moment >= 1.")
      end if

   end subroutine

   pure subroutine pbeterm_allocations(self)
   !! Allocator for arrays of `pbeterm` class.
      class(pbeterm), intent(inout) :: self
         !! object

      if (associated(self%grid)) then
         allocate (self%udot(self%grid%ncells))
      else
         call self%error_msg("Allocation failed due to missing grid.")
      end if

   end subroutine

   pure subroutine particleterm_allocations(self)
   !! Allocator for arrays of `particleterm` class.
      class(particleterm), intent(inout) :: self
         !! object

      ! Call parent method
      call self%pbeterm_allocations()

      ! Do own allocations
      if (associated(self%grid)) then
         allocate (self%udot_birth(self%grid%ncells), self%udot_death(self%grid%ncells))
      else
         call self%error_msg("Allocation failed due to missing grid.")
      end if

   end subroutine

end module pbepack_basetypes
