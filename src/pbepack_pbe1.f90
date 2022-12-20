module pbepack_pbe1
!! Derived types and procedures to solve 1D PBEs.
   use pbepack_kinds
   use pbepack_basetypes, only: base, pbeterm
   use pbepack_agg1, only: aggterm, afnc_t
   use pbepack_break1, only: breakterm, bfnc_t, dfnc_t
   use pbepack_growth1, only: growthterm, gfnc_t
   use pbepack_quadratures, only: quadgrid1
   use hrweno_grids, only: grid1
   use stdlib_optval, only: optval
   use odepack_mod, only: lsoda_class
   implicit none
   private

   public :: pbe, pbesol

   type, extends(pbeterm) :: pbe
   !! 1D PBE class.
      type(aggterm) :: agg
        !! aggregation object
      type(breakterm) :: break
         !! breakage object
      type(growthterm) :: growth
         !! growth object
   contains
      procedure, pass(self) :: eval => pbe_eval
      procedure, pass(self) :: integrate => pbe_integrate
   end type pbe

   type, extends(base) :: pbesol
   !! 1D PBE solution class.
      real(rk), allocatable :: t(:)
        !! vector(npoints) of time points
      real(rk), allocatable :: u(:, :)
        !! array(ncells,npoints) of /(\bar{u}_i(t_j)/)
      integer :: nfev = 0
         !! number of evaluations of /(d\bar{u}/dt/)
      logical :: success = .false.
         !! flag indicating if integration was successfull
   end type pbesol

   abstract interface
      pure real(rk) function ic_t(x)
      !! Interface initial condition of 1D PBE.
         import rk
         real(rk), intent(in) :: x
            !! internal coordinate
      end function
   end interface

   interface pbe
      module procedure :: pbe_init
   end interface pbe

   interface pbesol
      module procedure :: pbesol_init
   end interface pbesol

contains

   type(pbe) function pbe_init(grid, gfnc, afnc, bfnc, dfnc, moment, update_a, &
                               update_b, update_d, name) result(self)
   !! Initialize `pbe1` object.
      type(grid1), intent(in), target :: grid
         !! `grid1` object
      procedure(gfnc_t), optional :: gfnc
         !! growth rate  function, \( g(x,y) \)
      procedure(afnc_t), optional :: afnc
         !! aggregation frequency function, \( a(x,x',y) \)
      procedure(bfnc_t), optional :: bfnc
         !! breakage frequency function, \( b(x,y) \)
      procedure(dfnc_t), optional :: dfnc
         !! daughter distribution function, \( d(x,x',y) \)
      integer, intent(in), optional :: moment
         !! moment of \( x \) to be preserved upon aggregation/breakage (default=1)
      logical, intent(in), optional :: update_a
         !! flag to select if \( a(x,x',y) \) is to be reevaluated at each step (default=true)
      logical, intent(in), optional :: update_b
         !! flag to select if \( b(x,y) \) is to be reevaluated at each step (default=true)
      logical, intent(in), optional :: update_d
         !! flag to select if \( d(x,x',y) \) is to be reevaluated at each step (default=true)
      character(*), intent(in), optional :: name
         !! name (default="pbe")
      integer :: moment_

      ! Set properties
      call self%set_name(name, "pbe")
      call self%set_grid(grid)

      ! Init growth term
      if (present(gfnc)) then
         self%growth = growthterm(grid, gfnc, k=3, name=name//"-growth")
      end if

      ! Init aggregation term
      moment_ = optval(moment, 1)
      if (present(afnc)) then
         self%agg = aggterm(grid, afnc, moment_, optval(update_a, .true.), name=name//"-agg")
      end if

      ! Init breakage term
      if (present(bfnc) .neqv. present(dfnc)) then
         call self%error_msg("Arguments 'bfnc' and 'dfnc' must both be present or absent.")
      else if (present(bfnc) .and. present(dfnc)) then
         self%break = breakterm(grid, bfnc, dfnc, moment_, optval(update_b, .true.), &
                                optval(update_d, .true.), name=name//"-break")
      end if

      ! Init source term

      ! Allocate stuff
      call self%pbeterm_allocations()
      self%inited = .true.

   end function pbe_init

   pure subroutine pbe_eval(self, u, y, udot)
   !! Evaluate total rate of change at a given instant.
      class(pbe), intent(inout) :: self
         !! object
      real(rk), intent(in) :: u(:)
         !! cell-average number density, \( \bar{u} \)
      real(rk), intent(in) :: y(:)
         !! environment vector, \( y \)
      real(rk), intent(out), optional :: udot(:)
         !!  total rate of change, \( d\bar{u}/dt \)

      call self%check_inited()

      associate (udot_ => self%udot)

         udot_ = ZERO
         if (self%agg%inited) then
            call self%agg%eval(u, y)
            udot_ = udot_ + self%agg%udot
         end if

         if (self%break%inited) then
            call self%break%eval(u, y)
            udot_ = udot_ + self%break%udot
         end if

         if (self%growth%inited) then
            call self%growth%eval(u, y)
            udot_ = udot_ + self%growth%udot
         end if

         if (present(udot)) udot = udot_

      end associate

   end subroutine pbe_eval

   type(pbesol) function pbe_integrate(self, ic, times, atol, rtol, verbose) result(res)
   !! Integrate PBE using LSODA as ODE solver.
      class(pbe), intent(inout) :: self
         !! object
      procedure(ic_t) :: ic
         !! initial condition, \( u_0(x) \)
      real(rk), intent(in) :: times(:)
         !! time sequence for which output is wanted; the first value of `times` must be the
         !! initial time
      real(rk), intent(in), optional :: atol
         !! absolute tolerance (default=1e-6)
      real(rk), intent(in), optional :: rtol
         !! relative tolerance (default=1e-5)
      logical, intent(in), optional :: verbose
         !! verbose flag (default=false)

      type(lsoda_class) :: ode
      real(rk) :: u(self%grid%ncells), t, tout
      integer :: i, npoints, istate, itask

      ! Check inputs
      npoints = size(times)
      if (npoints < 2) then
         call self%error_msg("Invalid 'times' dimension. Valid range: size(times) >= 2.")
      end if

      ! Evaluate initial condition
      u = quadgrid1(ic, self%grid, average=.true.)

      ! Init ode solver object
      call ode%initialize(rhs, self%grid%ncells, istate=istate)
      if (istate < 0) then
         call self%error_msg("Initialization of ODE solver failed: "//ode%error_message)
      end if

      ! Init pbesolution object
      res = pbesol(self%grid%ncells, npoints)

      ! Integrate
      itask = 1
      istate = 1
      t = times(1)
      do i = 1, npoints
         tout = times(i)
         call ode%integrate(u, t, tout, optval(rtol, 1e-5_rk), [optval(atol, 1e-6_rk)], &
                            itask, istate)
         if (istate < 0) then
            call self%error_msg("Integration failed: "//ode%error_message)
         end if
         res%t(i) = t
         res%u(:, i) = u
      end do
      res%success = .true.

   contains

      subroutine rhs(this, neq, t_, v, vdot, ierr)
         class(lsoda_class), intent(inout) :: this
         integer, intent(in) :: neq
         real(dp), intent(in) :: t_, v(neq)
         real(dp), intent(out) :: vdot(neq)
         integer, intent(out) :: ierr
         call self%eval(v, [ZERO], vdot)
         ierr = 0
      end subroutine

   end function pbe_integrate

   pure type(pbesol) function pbesol_init(ncells, npoints) result(res)
   !! Initialize `pbesolution` object.
      integer, intent(in) :: ncells
         !! number of grid cells
      integer, intent(in) :: npoints
         !! number of time points

      allocate (res%t(npoints), res%u(ncells, npoints))
      res%inited = .true.

   end function

end module pbepack_pbe1
