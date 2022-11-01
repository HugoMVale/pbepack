module break1
!! Thi modules implements derived types and procedures to compute the breakage term
!! for 1D PBEs.
   use real_kinds
   use basetypes, only: particleterm
   use grids, only: grid1
   use stdlib_optval, only: optval

   implicit none
   private

   public :: breakterm

   type, extends(particleterm) :: breakterm
   !! Breakage term class.
      procedure(bfnc1_t), nopass, pointer :: bfnc => null()
         !! breakage frequency
      procedure(dfnc1_t), nopass, pointer :: dfnc => null()
         !! daughter distribution
      real(rk), allocatable :: b(:)
         !! vector of breakage frequencies
      real(rk), allocatable :: d(:, :)
         !! matrix of daughter probabilities
      real(rk), allocatable :: sink(:)
         !! vector of sink rates
   contains
      procedure, pass(self) :: eval => breakterm_eval
   end type breakterm

   abstract interface
      pure real(rk) function bfnc1_t(x, y)
      !! Breakage frequency for 1D system
         import :: rk
         real(rk), intent(in) :: x
            !! internal coordinate of particle
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function

      pure real(rk) function dfnc1_t(xd, xo, y)
      !! Daughter distribution for 1D system
         import :: rk
         real(rk), intent(in) :: xd
            !! internal coordinate of daughter particle
         real(rk), intent(in) :: xo
            !! internal coordinate of original particle
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function
   end interface

   interface breakterm
      module procedure :: breakterm_init
   end interface breakterm

contains

   type(breakterm) function breakterm_init(bfnc, dfnc, moment, grid, name) result(self)
   !! Initialize 'breakterm' object.
      procedure(bfnc1_t) :: bfnc
         !! breakage frequency, \( b(x,y) \)
      procedure(dfnc1_t) :: dfnc
         !! daughter distribution, \( d(x,x',y) \)
      integer, intent(in) :: moment
         !! moment of 'x' to be conserved upon aggregation
      type(grid1), intent(in), target, optional :: grid
         !! grid1 object
      character(*), intent(in), optional :: name
         !! object name

      self%bfnc => bfnc
      self%dfnc => dfnc

      call self%set_moment(moment)

      if (present(grid)) then
         call self%set_grid(grid)
         associate (nc => self%grid%ncells)
            allocate (self%b(nc), self%d(nc, nc), self%sink(nc))
         end associate
         self%inited = .true.
         self%msg = "Initialization completed successfully"
      else
         self%msg = "Missing 'grid'."
      end if

      self%name = optval(name, "")

   end function breakterm_init

   pure subroutine breakterm_eval(self, np, y, result, birth, death)
   !! Evaluate the rate of breakage at a given instant, using the technique described in
   !! Section 3.3 of Kumar and Ramkrishna (1996). The birth/source term is computed according
   !! to a slightly different procedure: the summation is done from the perspective of the
   !! original particles (the ones that break up), which allows for a simpler and more
   !! efficient calculation of the particle split fractions.
      class(breakterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: np(:)
         !! vector(ncells) with number of particles in cell 'i'
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), intent(out), optional :: result(:)
         !! vector(ncells) with net rate of change (birth-death)
      real(rk), intent(out), optional :: birth(:)
         !! vector(ncells) with birth term
      real(rk), intent(out), optional :: death(:)
         !! vector(ncells) with death term

      real(rk) :: weight(2)
      integer :: i, k

      associate (nc => self%grid%ncells, x => self%grid%center, b => self%b, &
                 result_ => self%result, sink => self%sink)

         ! Evaluate breakage frequency for all particle sizes
         do concurrent(i=1:nc)
            b(i) = self%bfnc(x(i), y)
         end do

         ! Death/sink term
         sink = b*np
         if (present(death)) death = sink

         ! Birth/source term
         result_ = ZERO
         do k = 2, nc
            do i = 1, k - 1
               weight = daughter_split(i, k)
               result_(i) = result_(i) + weight(1)*sink(k)
               result_(i + 1) = result_(i + 1) + weight(2)*sink(k)
            end do
         end do

         if (present(birth)) birth = result_
         result_ = result_ - sink

         if (present(result)) result = result_

      end associate

   contains

      pure function daughter_split(i_, k_) result(res)
      !! Evaluate the fraction of cell k assigned to cells i and (i+1).
      !! Based on Equation 27, but improved to consider the contribution to cell 1 arising
      !! from the integral between the left boundary of the grid domain and the center of cell
      !! 1. For this small region, we can only chose to preserve number or mass (not both). The
      !! current algorithm implements mass conservation.
         integer, intent(in) :: i_, k_
         real(rk) :: res(2)
         real(rk) :: d0, dm

         associate (x => self%grid%center, xo => self%grid%center(k_), m => self%moment)
            d0 = dfnc_integral(x(i_), x(i_ + 1), xo, 0)
            dm = dfnc_integral(x(i_), x(i_ + 1), xo, m)
            res(1) = (d0*x(i_ + 1)**m - dm)/(x(i_ + 1)**m - x(i_)**m)
            res(2) = d0 - res(1)
            if (i_ == 1) then
               !res(1) = res(1) + dfnc_integral(self%grid%left(1), x(1), xo, 0)
               res(1) = res(1) + dfnc_integral(self%grid%left(1), x(1), xo, m)/x(1)**m
            end if
         end associate

      end function daughter_split

      pure real(rk) function dfnc_integral(xa, xb, xo, m_) result(res)
      !! Numerical approximation to \( \int_{x_a}^{x_b} x^m dfnc(x,x_o,y)dx \).
      !! Equation 28.
      !! @todo: The quadrature is done with Simpson's rule, but perhaps a 2nd order method
      !! such as the trapezoidal rule would work just fine.
         real(rk), intent(in) :: xa, xb, xo
         integer, intent(in) :: m_

         real(rk) :: xc
         xc = (xa + xb)/2
         res = ( &
               xa**m_*self%dfnc(xa, xo, y) + &
               xb**m_*self%dfnc(xb, xo, y) + &
               xc**m_*self%dfnc(xc, xo, y)*4 &
               )*(xb - xa)/6
      end function dfnc_integral

   end subroutine breakterm_eval

end module break1
