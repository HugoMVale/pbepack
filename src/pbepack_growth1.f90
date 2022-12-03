module pbepack_growth1
!! This module implements derived types and procedures to compute the growth term
!! for 1D PBEs.
   use pbepack_kinds
   use pbepack_basetypes, only: pbeterm
   use hrweno_grids, only: grid1
   use hrweno_weno, only: weno
   implicit none
   private

   public :: growthterm

   type, extends(pbeterm) :: growthterm
   !! Growth term class.
      private
      procedure(gfnc1_t), nopass, pointer :: gfnc => null()
         !! growth rate function
      type(weno) :: wenorec
         !! WENO reconstruction
      real(rk), allocatable :: uleft(:)
         !! left reconstruction of number density function
      real(rk), allocatable :: uright(:)
         !! right reconstruction of number density function
      real(rk), allocatable :: flux_edges(:)
         !! flux at grid edges
   contains
      procedure, pass(self), public :: eval => growthterm_eval
      procedure, pass(self), public :: init2 => growthterm_init2
   end type growthterm

   abstract interface
      pure real(rk) function gfnc1_t(x, y)
      !! Growth rate for 1D system
         import :: rk
         real(rk), intent(in) :: x
            !! internal coordinate
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function
   end interface

   interface growthterm
      module procedure :: growthterm_init
   end interface growthterm

contains

   type(growthterm) function growthterm_init(gfnc, grid, name) result(self)
   !! Initialize 'growthterm' object.
      procedure(gfnc1_t) :: gfnc
         !! growth rate  function, \( g(x,y) \)
      type(grid1), intent(in), optional :: grid
         !! 'grid1' object
      character(*), intent(in), optional :: name
         !! name

      self%gfnc => gfnc
      call self%set_name(name)
      if (present(grid)) call self%init2(grid)

   end function growthterm_init

   subroutine growthterm_init2(self, grid)
   !! Initialize(2) 'growthterm' object.
      class(growthterm), intent(inout) :: self
         !! object
      type(grid1), intent(in) :: grid
         !! 'grid1' object

      call self%set_grid(grid)
      call self%pbeterm_allocations()
      associate (nc => self%grid%ncells, gx => self%grid)
         if (gx%scl == "linear") then
            self%wenorec = weno(ncells=nc, k=3)
         else
            self%wenorec = weno(ncells=nc, k=3, xedges=gx%edges)
         end if
         allocate (self%uleft(nc), self%uright(nc), self%flux_edges(0:nc))
      end associate
      self%inited = .true.

   end subroutine growthterm_init2

   pure subroutine growthterm_eval(self, u, y, udot)
   !! Evaluate rate of aggregation at a given instant.
      class(growthterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: u(:)
         !! vector(ncells) with cell-average number density
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), intent(out), optional :: udot(:)
         !! vector(ncells) with net rate of change

      integer :: i

      associate (nc => self%grid%ncells, uleft => self%uleft, uright => self%uright)
         ! Get reconstructed values at cell boundaries
         call self%wenorec%reconstruct(u, uleft, uright)

         ! ! Fluxes at interior cell boundaries
         ! ! One can use the Lax-Friedrichs or the Godunov method
         ! do concurrent(i=1:nc - 1)
         !    !fedges(i) = lax_friedrichs(flux, vr(i), vl(i+1), gx%r(i), t, alpha = 1._rk)
         !    fedges(i) = godunov(flux, vr(i), vl(i + 1), [gx%right(i)], t)
         ! end do

         ! ! Apply problem-specific flux constraints at domain boundaries
         ! fedges(0) = 0
         ! fedges(nc) = fedges(nc - 1)

         ! ! Evaluate du/dt
         ! udot = -(fedges(1:) - fedges(:nc - 1))/gx%width
      end associate

   contains

      pure real(rk) function gflux(u, x, t) result(res)
         real(rk), intent(in) :: u, x(:), t

         res = self%gfnc(x(1), y)*u

      end function

   end subroutine growthterm_eval

end module pbepack_growth1
