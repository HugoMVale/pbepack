module break1
!! Thi modules implements derived types and procedures to compute the breakage term
!! for 1D PBEs.
   use real_kinds
   use basetypes, only: particleterm
   use grids, only: grid1
   use quadratures, only: quad1_r1
   use stdlib_optval, only: optval

   implicit none
   private

   public :: breakterm

   type, extends(particleterm) :: breakterm
   !! Breakage term class.
      procedure(bf1_t), nopass, pointer :: bf => null()
         !! breakage frequency
      procedure(df1_t), nopass, pointer :: df => null()
         !! daughter distribution
      real(rk), allocatable :: b(:)
         !! vector of breakage frequencies
      real(rk), allocatable :: d(:, :)
         !! matrix of daughter probabilities
   contains
      procedure, pass(self) :: eval => breakterm_eval
   end type breakterm

   abstract interface
      pure real(rk) function bf1_t(x, y)
      !! Breakage frequency for 1D system
         import :: rk
         real(rk), intent(in) :: x
            !! internal coordinate of particle
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function

      pure real(rk) function df1_t(xd, xo, y)
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

   type(breakterm) function breakterm_init(bf, df, moment, grid, name) result(self)
   !! Initialize 'breakterm' object.
      procedure(bf1_t) :: bf
         !! breakage frequency, \( b(x,y) \)
      procedure(df1_t) :: df
         !! daughter distribution, \( d(x,x',y) \)
      integer, intent(in) :: moment
         !! moment of 'x' to be conserved upon aggregation
      type(grid1), intent(in), target, optional :: grid
         !! grid1 object
      character(*), intent(in), optional :: name
         !! object name

      self%bf => bf
      self%df => df

      call self%set_moment(moment)

      if (present(grid)) then
         call self%set_grid(grid)
         associate (nc => self%grid%ncells)
            allocate (self%b(nc), self%d(nc, nc), self%source(nc), self%sink(nc))
         end associate
         self%inited = .true.
         self%msg = "Initialization completed successfully"
      else
         self%msg = "Missing 'grid'."
      end if

      self%name = optval(name, "")

   end function breakterm_init

   pure subroutine breakterm_eval(self, np, y, total, source, sink)
   !! Evaluate rate of breakage at a given instant.
      class(breakterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: np(:)
         !! vector(ncells) with number of particles in cell 'i'
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), intent(out), optional :: total(:)
         !! vector(ncells) with sum of birth and death terms
      real(rk), intent(out), optional :: source(:)
         !! vector(ncells) with source (birth, +) term
      real(rk), intent(out), optional :: sink(:)
         !! vector(ncells) with sink (death, -) term

      integer :: i, k

      associate (nc => self%grid%ncells, b => self%b, &
                 source_ => self%source, sink_ => self%sink)

         ! Evaluate breakage frequency for all particle sizes
         do concurrent(i=1:nc)
            b(i) = self%bf(self%grid%center(i), y)
         end do

         ! Sink/death term
         sink_ = b*np

         ! Source/birth term
         source_ = ZERO
         !do concurrent(i=1:nc)
         do i = 1, nc
            do k = i, nc
               source_(i) = source_(i) + weight()*sink_(k)
            end do
         end do

         if (present(source)) source = source_
         if (present(sink)) sink = sink_
         if (present(total)) total = source_ + sink_

      end associate

   contains

      pure real(rk) function weight()
      !! Weight fraction of particle k assigned to particle i.

         real(rk) :: term(2)

         associate (x => self%grid%center, m => self%moment)
            if (i == k) then
               term(1) = ZERO
            else
               term(1) = (dfint(i, 0)*x(i + 1)**m - dfint(i, m))/ &
                         (x(i + 1)**m - x(i)**m)
            end if

            if (i == 1) then
               term(2) = ZERO
            else
               term(2) = (dfint(i - 1, 0)*x(i - 1)**m - dfint(i - 1, m))/ &
                         (x(i - 1)**m - x(i)**m)
            end if
         end associate

         weight = sum(term)

      end function weight

      pure real(rk) function dfint(i_, m_) result(res)
         integer, intent(in) :: i_, m_
         real(rk) :: xa, xb, xm
         associate (x => self%grid%center)
            xa = x(i_)
            xb = x(i_ + 1)
            xm = (xa + xb)/2
            res = ( &
                  xa**m_*self%df(xa, x(k), y) + &
                  xb**m_*self%df(xb, x(k), y) + &
                  xm**m_*self%df(xm, x(k), y)*4 &
                  )*(xb - xa)/6
         end associate

      end function dfint

   end subroutine breakterm_eval

end module break1
