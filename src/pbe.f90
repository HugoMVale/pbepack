module pbe
   use real_kinds, only: rk
   use basetypes, only: base
   use grid, only: grid1
   use agg1, only: aggterm
   use stdlib_optval, only: optval
   implicit none
   private

   public :: pbe1

   type, extends(base) :: pbe1
   !! 1D PBE class.
      integer :: ncells = 0
        !! number of cells
      type(grid1), pointer :: grid => null()
        !! grid object
      !class(pbeterm), allocatable :: term(:)
        !! vector of pbe terms
      procedure(ic1_t), nopass, pointer :: ic => null()
        !! initial condition
   contains
      !procedure, pass(self) :: eval => pbe_eval
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

   type(pbe1) function pbe1_init(grid, ic, name) result(self)
   !! Initialize 'pbe1' object.
      type(grid1), intent(in), target :: grid
        !! grid object
      procedure(ic1_t), optional :: ic
        !! initial condition, \( f_0(x) \)
      character(*), intent(in), optional :: name
        !! pbe name
      self%ncells = grid%ncells
      self%grid => grid
      if (present(ic)) self%ic => ic
      self%name = optval(name, "")
   end function pbe1_init

end module pbe
