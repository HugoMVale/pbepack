module pbepack_pbe1
   use pbepack_kinds
   use pbepack_basetypes, only: pbeterm
   use pbepack_agg1, only: aggterm, afnc_t
   use pbepack_break1, only: breakterm, bfnc_t, dfnc_t
   use pbepack_growth1, only: growthterm, gfnc_t
   use hrweno_grids, only: grid1
   use stdlib_optval, only: optval
   implicit none
   private

   public :: pbe1

   type, extends(pbeterm) :: pbe1
   !! 1D PBE class.
      type(aggterm) :: agg
        !! aggregation object
      type(breakterm) :: break
         !! breakage object
      type(growthterm) :: growth
         !! growth object
      procedure(ic_t), nopass, pointer :: ic => null()
        !! initial condition
   contains
      procedure, pass(self) :: eval => pbe_eval
      procedure, pass(self) :: integrate => pbe_integrate
   end type pbe1

   abstract interface
      pure real(rk) function ic_t(x)
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

   type(pbe1) function pbe1_init(grid, gfnc, afnc, bfnc, dfnc, moment, update_a, &
                                 update_b, update_d, name) result(self)
   !! Initialize 'pbe1' object.
      type(grid1), intent(in), target :: grid
         !! grid object
      procedure(gfnc_t), optional :: gfnc
         !! growth rate  function, \( g(x,y) \)
      procedure(afnc_t), optional :: afnc
         !! aggregation frequency function, \( a(x,x',y) \)
      procedure(bfnc_t), optional :: bfnc
         !! breakage frequency function, \( b(x,y) \)
      procedure(dfnc_t), optional :: dfnc
         !! daughter distribution function, \( d(x,x',y) \)
      integer, intent(in), optional :: moment
         !! moment of \( x \) to be preserved upon aggregation/breakage (> 0)
      logical, intent(in), optional :: update_a
         !! flag to select if \( a(x,x',y) \) is to be reevaluated at each step
      logical, intent(in), optional :: update_b
         !! flag to select if \( b(x,y) \) is to be reevaluated at each step
      logical, intent(in), optional :: update_d
         !! flag to select if \( d(x,x',y) \) is to be reevaluated at each step
      character(*), intent(in), optional :: name
         !! name
      integer :: moment_

      ! Grid
      call self%set_grid(grid)

      ! Growth term
      if (present(gfnc)) then
         self%growth = growthterm(grid, gfnc, k=3, name=name//":growth")
      end if

      ! Aggregation
      moment_ = optval(moment, 1)
      if (present(afnc)) then
         self%agg = aggterm(grid, afnc, moment_, optval(update_a, .true.), name=name//":agg")
      end if

      ! Breakage
      if (present(bfnc) .neqv. present(dfnc)) then
         call self%error_msg("Arguments 'bfnc' and 'dfnc' must _both_ be present or absent.")
      else if (present(bfnc) .and. present(dfnc)) then
         self%break = breakterm(grid, bfnc, dfnc, moment_, optval(update_b, .true.), &
                                optval(update_d, .true.), name=name//":break")
      end if

      ! Source

      call self%set_name(name)

   end function pbe1_init

   pure subroutine pbe_eval(self, np, y, udot)
   !! Evaluate rate of aggregation at a given instant.
      class(pbe1), intent(inout) :: self
         !! object
      real(rk), intent(in) :: np(:)
         !! vector(ncells) with number of particles in cell 'i'
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), intent(out), optional :: udot(:)
         !! vector(ncells) with net rate of change (birth-death)

      udot = ZERO
      if (self%agg%inited) then
         call self%agg%eval(np, y)
         udot = udot + self%agg%udot
      end if

      if (self%break%inited) then
         call self%break%eval(np, y)
         udot = udot + self%break%udot
      end if

      if (self%growth%inited) then
         call self%growth%eval(np, y)
         udot = udot + self%growth%udot
      end if

      ! ! Net rate
      ! if (present(result)) result = result_

   end subroutine pbe_eval

   pure subroutine pbe_integrate(self, np, y, result)
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
      if (self%agg%inited) then
         call self%agg%eval(np, y)
         result = result + self%agg%udot
      end if

      if (self%break%inited) then
         call self%break%eval(np, y)
         result = result + self%break%udot
      end if

      if (self%growth%inited) then
         call self%growth%eval(np, y)
         result = result + self%growth%udot
      end if

      ! ! Net rate
      ! if (present(result)) result = result_

   end subroutine pbe_integrate

end module pbepack_pbe1
