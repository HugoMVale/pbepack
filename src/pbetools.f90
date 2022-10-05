module pbetools
   use, intrinsic :: iso_fortran_env, only: real64
   use grid, only: grid1
   ! use auxfunctions, only: delta_kronecker
   ! use combtypes, only: comblist, combarray
   implicit none
   private

   public :: pbe, pbeterm

   integer, parameter :: rk = real64

   type :: pbe
   !! 1D PBE class.
      character(:), allocatable :: name
        !! pbe name
      character(:), allocatable :: msg
        !! error message
      integer :: ierr = 0
        !! error status
      integer :: ncells = 0
        !! number of cells
      type(grid1), pointer :: grid
        !! grid object
      class(pbeterm), allocatable :: term(:)
        !! vector of pbe terms
      procedure(ic1d), pointer :: ic => null()
        !! initial condition
   contains
      procedure, pass(self):: init => pbe_init
      !procedure, pass(self) :: eval => pbe_eval
      !procedure, pass(self) :: integrate => pbe_integrate
   end type pbe

   type, abstract :: pbeterm
   !! Abstract 1D PBE term class (e.g., aggregation, growth, etc.)
      character(:), allocatable :: name
        !! variable name
      character(:), allocatable :: msg
        !! error message
      integer :: ierr = 0
        !! error status
      integer :: ncells = 0
        !! number of cells
      type(grid1), pointer :: grid => null()
        !! pointer to grid
   end type pbeterm

   abstract interface
      elemental real(rk) function ic1d(self, x)
      !! aggregation kernel for 1D system
         import rk, pbe
         class(pbe), intent(in) :: self
         real(rk), intent(in) :: x
      end function

      pure real(rk) function aggfun1d(xa, xb, t, y)
      !! aggregation kernel for 1D system
         import :: rk
         real(rk), intent(in) :: xa
          !! internal coordinate of particle a
         real(rk), intent(in) :: xb
          !! internal coordinate of particle b
         real(rk), intent(in) :: t
          !! time
         real(rk), intent(in) :: y(:)
          !! environment vector
      end function
   end interface

contains

   pure subroutine pbe_init(self, grid, agg)
   !! Initialize 'pbe' object.
      class(pbe), intent(inout) :: self
        !! object
      type(grid1), intent(inout), target :: grid
        !! grid
      procedure(aggfun1d) :: agg
        !! aggregation frequency function, a(x,x',t,y)
      self%ncells = grid%ncells
      self%grid => grid

   end subroutine

end module pbetools
