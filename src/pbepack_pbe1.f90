module pbepack_pbe1
   use pbepack_kinds
   use pbepack_basetypes, only: base
   use pbepack_agg1, only: aggterm
   use pbepack_break1, only: breakterm
   use hrweno_grids, only: grid1
   implicit none
   private

   public :: pbe1

   type, extends(base) :: pbe1
   !! 1D PBE class.
      type(aggterm), pointer :: agg => null()
        !! aggregation object
      type(breakterm), pointer :: break => null()
         !! breakage object
      procedure(ic1_t), nopass, pointer :: ic => null()
        !! initial condition
   contains
      procedure, pass(self) :: eval => pbe_eval
      !procedure, pass(self) :: integrate => pbe_integrate
   end type pbe1

   abstract interface
      pure real(rk) function ic1_t(x)
      !! Interface initial condition of PBE
         import rk
         real(rk), intent(in) :: x
          !! internal coordinate
      end function
   end interface

   interface pbe1
      module procedure :: pbe1_init
   end interface pbe1

contains

   type(pbe1) function pbe1_init(grid, agg, break, ic, name) result(self)
   !! Initialize 'pbe1' object.
      type(grid1), intent(in), target :: grid
         !! grid object
      type(aggterm), intent(in), target, optional :: agg
         !! aggterm object
      type(breakterm), intent(in), target, optional :: break
         !! breakterm object
      procedure(ic1_t), optional :: ic
        !! initial condition, \( f_0(x) \)
      character(*), intent(in), optional :: name
        !! pbe name

      call self%set_grid(grid)
      if (present(agg)) self%agg => agg
      if (present(break)) self%break => break
      if (present(ic)) self%ic => ic
      call self%set_name(name)

   end function pbe1_init

   pure subroutine pbe_eval(self, np, y, result)
   !! Evaluate rate of aggregation at a given instant.
      class(pbe1), intent(inout) :: self
         !! object
      real(rk), intent(in) :: np(:)
         !! vector(ncells) with number of particles in cell 'i'
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), intent(out), optional :: result(:)
         !! vector(ncells) with net rate of change (birth-death)

      result = ZERO
      if (associated(self%agg)) then
         call self%agg%eval(np, y)
         result = result + self%agg%result
      end if

      if (associated(self%break)) then
         call self%break%eval(np, y)
         result = result + self%break%result
      end if

      ! ! Net rate
      ! if (present(result)) result = result_

   end subroutine pbe_eval

end module pbepack_pbe1
