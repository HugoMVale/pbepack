module pbepack_growth1
!! This module implements derived types and procedures to compute the growth term
!! for 1D PBEs.
   use pbepack_kinds
   use pbepack_basetypes, only: pbeterm
   use hrweno_grids, only: grid1
   use hrweno_weno, only: weno
   use hrweno_fluxes, only: godunov
   implicit none
   private

   public :: growthterm, gfnc_t

   type, extends(pbeterm) :: growthterm
   !! Growth term class.
      private
      procedure(gfnc_t), nopass, pointer :: gfnc => null()
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
   end type growthterm

   abstract interface
      pure real(rk) function gfnc_t(x, y)
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

   type(growthterm) function growthterm_init(grid, gfnc, k, name) result(self)
   !! Initialize 'growthterm' object.
      type(grid1), intent(in) :: grid
         !! 'grid1' object
      procedure(gfnc_t) :: gfnc
         !! growth rate  function, \( g(x,\mathbold{y}) \)
      integer, intent(in) :: k
         !! 2*(k-1) order of the WENO reconstruction (1 <= k <= 3)
      character(*), intent(in), optional :: name
         !! name

      call self%set_grid(grid)
      self%gfnc => gfnc
      call self%set_name(name)
      call self%pbeterm_allocations()
      associate (nc => self%grid%ncells, gx => self%grid)
         if (gx%scale == "linear") then
            self%wenorec = weno(nc, k=k)
         else
            self%wenorec = weno(nc, k=k, xedges=gx%edges)
         end if
         allocate (self%uleft(nc), self%uright(nc), self%flux_edges(0:nc))
      end associate
      self%inited = .true.

   end function growthterm_init

   pure subroutine growthterm_eval(self, u, y, udot)
   !! Evaluate rate of particle growth at a given instant.
      class(growthterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: u(:)
         !! cell-average number density, \( \bar{u} \)
      real(rk), intent(in) :: y(:)
         !! environment vector, \( \mathbold{y} \)
      real(rk), intent(out), optional :: udot(:)
         !! net rate of change, \( d\bar{u}/dt \)

      integer :: i

      associate (nc => self%grid%ncells, gx => self%grid, uleft => self%uleft, &
                 uright => self%uright, flux_edges => self%flux_edges)
         ! Get reconstructed values at cell boundaries
         call self%wenorec%reconstruct(u, uleft, uright)

         ! Fluxes at interior cell boundaries
         do concurrent(i=1:nc - 1)
            flux_edges(i) = godunov(gflux, uright(i), uleft(i + 1), [gx%right(i)], ZERO)
         end do

         ! Apply problem-specific flux constraints at domain boundaries
         flux_edges(0) = 0
         flux_edges(nc) = flux_edges(nc - 1)

         ! Evaluate du/dt
         self%udot = -(flux_edges(1:) - flux_edges(:nc - 1))/gx%width
         if (present(udot)) udot = self%udot

      end associate

   contains

      pure real(rk) function gflux(u_, x_, t_) result(res)
         real(rk), intent(in) :: u_, x_(:), t_

         res = self%gfnc(x_(1), y)*u_

      end function

   end subroutine growthterm_eval

end module pbepack_growth1
