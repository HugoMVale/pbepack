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
      procedure, pass(self), public :: eval => pbe_eval
      procedure, pass(self), public :: integrate => pbe_integrate
   end type pbe

   type, extends(base) :: pbesol
   !! 1D PBE solution class.
      real(rk), allocatable :: t(:)
         !! vector(npoints) of time points
      real(rk), allocatable :: u(:, :)
         !! array(ncells,npoints) of \( \bar{u}_i(t_k) \)
      real(rk), allocatable :: y(:, :)
         !! array(nenvs,npoints) of \( {y}_j(t_k) \)
      integer :: nfev = 0
         !! number of evaluations of \( d\bar{\textbf{u}}/dt \)
      logical :: success = .false.
         !! flag indicating if integration was successfull
   end type pbesol

   abstract interface
      pure real(rk) function u0fnc_t(x)
      !! Initial condition of 1D PBE.
         import rk
         real(rk), intent(in) :: x
            !! internal coordinate
      end function

      pure function ydotfnc_t(t, y)
      !! Derivative environment vector.
         import rk
         real(rk), intent(in) :: t
            !! time
         real(rk), intent(in) :: y(:)
            !! environment vector
         real(rk) :: ydotfnc_t(size(y))
      end function

   end interface

   interface pbe
      module procedure :: pbe_init
   end interface pbe

   interface pbesol
      module procedure :: pbesol_init
   end interface pbesol

contains

   type(pbe) function pbe_init(grid, g, a, b, d, moment, update_a, update_b, update_d, &
                               name) result(self)
   !! Initialize `pbe` object.
      type(grid1), intent(in), target :: grid
         !! `grid1` object
      procedure(gfnc_t), optional :: g
         !! growth rate  function, \( g(x,\textbf{y}) \)
      procedure(afnc_t), optional :: a
         !! aggregation frequency function, \( a(x,x',\textbf{y}) \)
      procedure(bfnc_t), optional :: b
         !! breakage frequency function, \( b(x,\textbf{y}) \)
      procedure(dfnc_t), optional :: d
         !! daughter distribution function, \( d(x,x',\textbf{y}) \)
      integer, intent(in), optional :: moment
         !! moment of \( x \) to be preserved upon aggregation/breakage (default=1)
      logical, intent(in), optional :: update_a
         !! if `true`, \( a(x,x',\textbf{y}) \) is reevaluated at each step (default=true)
      logical, intent(in), optional :: update_b
         !! if `true`, \( b(x,\textbf{y}) \) is reevaluated at each step (default=true)
      logical, intent(in), optional :: update_d
         !! if `true`, \( d(x,x',\textbf{y}) \) is reevaluated at each step (default=true)
      character(*), intent(in), optional :: name
         !! name (default="pbe")
      integer :: moment_

      ! Set properties
      call self%set_name(name, "pbe")
      call self%set_grid(grid)

      ! Init growth term
      if (present(g)) then
         self%growth = growthterm(grid, g, k=3, name=name//"-growth")
      end if

      ! Init aggregation term
      moment_ = optval(moment, 1)
      if (present(a)) then
         self%agg = aggterm(grid, a, moment_, optval(update_a, .true.), name=name//"-agg")
      end if

      ! Init breakage term
      if (present(b) .neqv. present(d)) then
         call self%error_msg("Arguments 'b' and 'd' must both be present or absent.")
      else if (present(b) .and. present(d)) then
         self%break = breakterm(grid, b, d, moment_, optval(update_b, .true.), &
                                optval(update_d, .true.), name=name//"-break")
      end if

      ! Init source term

      ! Allocate stuff
      call self%pbeterm_allocations()
      self%inited = .true.

   end function pbe_init

   pure type(pbesol) function pbesol_init(npoints, ncells, nenvs) result(res)
   !! Initialize `pbesol` object.
      integer, intent(in) :: npoints
         !! number of time points, i.e. `size(times)`
      integer, intent(in) :: ncells
         !! number of grid cells, i.e. `size(u)`
      integer, intent(in) :: nenvs
         !! number of environment variables, i.e. `size(y)`

      allocate (res%t(npoints), res%u(ncells, npoints), res%y(nenvs, npoints))
      res%inited = .true.

   end function

   pure subroutine pbe_eval(self, u, y, udot)
   !! Evaluate total rate of change at a given instant.
      class(pbe), intent(inout) :: self
         !! object
      real(rk), intent(in) :: u(:)
         !! cell-average number density, \( \bar{u_i} \)
      real(rk), intent(in) :: y(:)
         !! environment vector, \( y_j \)
      real(rk), intent(out), optional :: udot(:)
         !! total rate of change, \( d\bar{u_i}/dt \)

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

   type(pbesol) function pbe_integrate(self, times, u0fnc, y0, ydotfnc, atol, rtol, &
                                       verbose) result(res)
   !! Integrate PBE using LSODA as ODE solver.
      class(pbe), intent(inout) :: self
         !! object
      real(rk), intent(in) :: times(:)
         !! time sequence for which output is wanted; the first value of `times` must be the
         !! initial time
      procedure(u0fnc_t) :: u0fnc
         !! initial condition of number density, \( u(x,t=0) \)
      real(rk), intent(in), optional :: y0(:)
         !! initial condition of environment vector, \( y_j(t=0) \)
      procedure(ydotfnc_t), optional :: ydotfnc
         !! derivative of environment vector, \( dy_j/dt \)
      real(rk), intent(in), optional :: atol
         !! absolute tolerance (default=1e-15)
      real(rk), intent(in), optional :: rtol
         !! relative tolerance (default=1e-5)
      logical, intent(in), optional :: verbose
         !! verbose flag (default=false)

      type(lsoda_class) :: ode
      real(rk), allocatable :: z(:)
      real(rk) :: t, tout
      integer :: i, ncells, nenvs, npoints, istate, itask

      ! Check `times`
      npoints = size(times)
      if (npoints < 2) then
         call self%error_msg("Invalid 'times' dimension. Valid range: size(times) >= 2.")
      end if

      ! Check `y0` and `ydotfnc`
      if (present(y0)) then
         nenvs = size(y0)
      else
         nenvs = 0
      end if
      if (present(ydotfnc) .and. nenvs == 0) then
         call self%error_msg("Dimension mismatch between 'y0' and 'ydotfnc'.")
      end if

      ! Allocate joint state vector
      ncells = self%grid%ncells
      allocate (z(ncells + nenvs))

      ! Init ode solver object
      call ode%initialize(zdotfnc, size(z), istate=istate)
      if (istate < 0) then
         call self%error_msg("Initialization of ODE solver failed: "//ode%error_message)
      end if

      ! Init pbesolution object
      res = pbesol(npoints, ncells, nenvs)

      ! Evaluate and assign initial condition
      z(:ncells) = quadgrid1(u0fnc, self%grid, average=.true.)
      if (nenvs > 0) z(ncells + 1:) = y0

      ! Integrate
      itask = 1
      istate = 1
      t = times(1)
      do i = 1, npoints
         tout = times(i)
         call ode%integrate(z, t, tout, &
                            optval(rtol, 1e-5_rk), [optval(atol, 1e-15_rk)], &
                            itask, istate)
         if (istate < 0) then
            call self%error_msg("Integration failed: "//ode%error_message)
         end if
         ! Store integration results in pbesol object
         res%t(i) = t
         res%u(:, i) = z(:ncells)
         if (nenvs > 0) res%y(:, i) = z(ncells + 1:)
      end do
      call ode%info(nfe=res%nfev)
      res%success = .true.

   contains

      subroutine zdotfnc(this, neq, t_, z_, zdot_, ierr)
         class(lsoda_class), intent(inout) :: this
         integer, intent(in) :: neq
         real(dp), intent(in) :: t_, z_(neq)
         real(dp), intent(out) :: zdot_(neq)
         integer, intent(out) :: ierr

         if (nenvs == 0) then
            call self%eval(z_, VOIDREAL, zdot_)
         else
            call self%eval(z_(:ncells), z_(ncells + 1:), zdot_(:ncells))
            if (present(ydotfnc)) then
               zdot_(ncells + 1:) = ydotfnc(t_, z_(ncells + 1:))
            end if
         end if

         ierr = 0

      end subroutine

   end function pbe_integrate

end module pbepack_pbe1
