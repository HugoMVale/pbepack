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
      procedure(bf1_t), nopass, pointer :: bf => null()
         !! breakage frequency
      procedure(df1_t), nopass, pointer :: df => null()
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
            allocate (self%b(nc), self%d(nc, nc), self%sink(nc))
         end associate
         self%inited = .true.
         self%msg = "Initialization completed successfully"
      else
         self%msg = "Missing 'grid'."
      end if

      self%name = optval(name, "")

   end function breakterm_init

   impure subroutine breakterm_eval(self, np, y, result, birth, death)
   !! Evaluate rate of breakage at a given instant.
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

      integer :: i, k

      associate (nc => self%grid%ncells, x => self%grid%center, b => self%b, &
                 result_ => self%result, sink => self%sink)

         ! Evaluate breakage frequency for all particle sizes
         do concurrent(i=1:nc)
            b(i) = self%bf(x(i), y)
         end do

         ! Death/sink term
         sink = b*np
         if (present(death)) death = sink

         ! Birth/source term
         result_ = ZERO
         !do concurrent(i=1:nc)
         do i = 1, nc
            do k = i, nc
               result_(i) = result_(i) + weight()*sink(k)
            end do
         end do
         if (present(birth)) birth = result_
         result_ = result_ - sink

         if (present(result)) result = result_

      end associate

   contains

      impure real(rk) function weight()
      !! Weight fraction of particle k assigned to particle i.

         real(rk) :: term(2)

         associate (x => self%grid%center, xc => self%grid%center(i), m => self%moment)
            if (i == 1) then
               term(1) = dfint(self%grid%left(i), x(i), 0)
            else
               term(1) = (dfint(x(i - 1), x(i), 0)*x(i - 1)**m - dfint(x(i - 1), x(i), m))/ &
                         (x(i - 1)**m - x(i)**m)
            end if

            if (i == k) then
               term(2) = ZERO
            else
               term(2) = (dfint(x(i), x(i + 1), 0)*x(i + 1)**m - dfint(x(i), x(i + 1), m))/ &
                         (x(i + 1)**m - x(i)**m)
            end if

            !print *, i, k, term

         end associate

         weight = sum(term)

      end function weight

      pure real(rk) function dfint(xa, xb, m_) result(res)
         real(rk), intent(in) :: xa, xb
         integer, intent(in) :: m_
         real(rk) :: xm
         xm = (xa + xb)/2
         associate (xk => self%grid%center(k))
            res = ( &
                  xa**m_*self%df(xa, xk, y) + &
                  xb**m_*self%df(xb, xk, y) + &
                  xm**m_*self%df(xm, xk, y)*4 &
                  )*(xb - xa)/6
         end associate

      end function dfint

   end subroutine breakterm_eval

end module break1
