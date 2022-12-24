module pbepack_flow
!! Derived types and procedures to compute the *flow* term for 1D PBEs.
   use pbepack_kinds, only: rk, EPS
   use pbepack_basetypes, only: pbeterm
   use hrweno_grids, only: grid1
   implicit none
   private

   public :: flowterm, uinfnc_t

   type, extends(pbeterm) :: flowterm
   !! Flow term class.
      private
      procedure(uinfnc_t), nopass, pointer :: uinfnc => null()
         !! inlet number density function
   contains
      procedure, pass(self), public :: eval => flowterm_eval
   end type flowterm

   abstract interface
      pure real(rk) function uinfnc_t(x, t)
      !! Growth rate for 1D system
         import :: rk
         real(rk), intent(in) :: x
            !! internal coordinate
         real(rk), intent(in) :: t
            !! time
      end function
   end interface

   interface flowterm
      module procedure :: flowterm_init
   end interface flowterm

contains

   type(flowterm) function flowterm_init(grid, name) result(self)
   !! Initialize `flowterm` object.
   !!$$
   !! \left ( \dot{u} \right )_{\mathrm{flow}} =
   !! \frac{q_{\mathrm{in}} u_{\mathrm{in}}(x,t) - q_{\mathrm{out}}u(x,t)}{V}
   !! - u(x,t) \frac{\dot{V}}{V}
   !!$$
      type(grid1), intent(in) :: grid
         !! `grid1` object
      character(*), intent(in), optional :: name
         !! name (default="growth-term")

      ! Set properties
      call self%set_name(name, "flow-term")
      call self%set_grid(grid)

      ! Allocate stuff
      call self%pbeterm_allocations()
      self%inited = .true.

   end function flowterm_init

   pure subroutine flowterm_eval(self, u, uin, qin, qout, v, vdot, udot)
   !! Evaluate rate of change due to flow at a given instant \(t\).
      class(flowterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: u(:)
         !! cell-average number density, \( \bar{u_i}(t) \)
      real(rk), intent(in) :: uin(:)
         !! inlet cell-average number density, \( \bar{u_{\mathrm{in},i}}(t) \)
      real(rk), intent(in) :: qin
         !! inlet flowrate, \( q_{\mathrm{in}}(t) \)
      real(rk), intent(in) :: qout
         !! outlet flowrate, \( q_{\mathrm{out}}(t) \)
      real(rk), intent(in) :: v
         !! system volume, \( V(t) \)
      real(rk), intent(in) :: vdot
         !! volume derivative, \( \dot{V}(t) \)
      real(rk), intent(out), optional :: udot(:)
         !! net rate of change, \( \dot{u_i}(t) \)

      call self%check_inited()

      self%udot = (qin*uin - qout*u)/(v + EPS) - u/(v + EPS)*vdot

      if (present(udot)) udot = self%udot

   end subroutine flowterm_eval

end module pbepack_flow
