module break1
!! Thi modules implements the derived types and procedures to compute the breakage term
!! for 1D PBEs.
   use real_kinds
   use basetypes, only: particleterm
   use grids, only: grid1
   use auxfunctions, only: delta_kronecker
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
      procedure(bf1_t), optional :: bf
         !! breakage frequency, \( b(x,y) \)
      procedure(df1_t), optional :: df
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

      if (present(name)) self%name = name

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

      real(rk) :: sumbi
      integer :: i, k

      associate (gx => self%grid, nc => self%grid%ncells, b => self%b, d => self%d, &
                 source_ => self%source, sink_ => self%sink)

         ! Evaluate breakage frequency for all particle sizes
         do concurrent(i=1:nc)
            b(i) = self%bf(gx%center(i), y)
         end do

         !

         ! Source/birth term
         do concurrent(i=1:nc)

            sumbi = 0._rk
            do k = i, nc

               sumbi = sumbi + b(k)*np(k)

            end do
            source_(i) = sumbi

         end do

         ! Sink/death term
         sink_ = np*b

         if (present(source)) source = source_
         if (present(sink)) sink = sink_
         if (present(total)) total = source_ + sink_

      end associate

   end subroutine breakterm_eval

end module break1
