module break1
   use real_kinds, only: rk
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
      !procedure, pass(self) :: eval => breakterm_eval
   end type breakterm

   abstract interface
      pure real(rk) function bf1_t(x, t, y)
      !! Breakage frequency for 1D system
         import :: rk
         real(rk), intent(in) :: x
            !! internal coordinate of particle
         real(rk), intent(in) :: t
            !! time
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function

      pure real(rk) function df1_t(xd, xo, t, y)
      !! Daughter distribution for 1D system
         import :: rk
         real(rk), intent(in) :: xd
            !! internal coordinate of daughter particle
         real(rk), intent(in) :: xo
            !! internal coordinate of original particle
         real(rk), intent(in) :: t
            !! time
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function
   end interface

   interface breakterm
      module procedure :: breakterm_init
   end interface breakterm

contains

   type(breakterm) function breakterm_init(bf, df, moment, grid) result(self)
   !! Initialize 'breakterm' object.
      procedure(bf1_t), optional :: bf
         !! breakage frequency, \( b(x,t,y) \)
      procedure(df1_t), optional :: df
         !! daughter distribution, \( d(x,x',t,y) \)
      integer, intent(in) :: moment
         !! moment of 'x' to be conserved upon aggregation
      type(grid1), intent(in), target, optional :: grid
         !! grid1 object

      self%bf => bf
      self%df => df

      call self%set_moment(moment)

      if (present(grid)) then
         call self%set_grid(grid)
         allocate (self%b(grid%ncells), self%d(grid%ncells, grid%ncells))
         self%inited = .true.
         self%msg = "Initialization completed successfully"
      else
         self%msg = "Missing 'grid'."
      end if

   end function breakterm_init

   ! pure subroutine aggterm_eval(self, np, t, y, total, birth, death)
   ! !! Evaluate rate of aggregation at a given instant.
   !    class(aggterm), intent(inout) :: self
   !       !! object
   !    real(rk), intent(in) :: np(:)
   !       !! vector(ncells) with number of particles in cell 'i'
   !    real(rk), intent(in) :: t
   !       !! time
   !    real(rk), intent(in) :: y(:)
   !       !! environment vector
   !    real(rk), intent(out), optional :: total(:)
   !       !! vector(ncells) with sum of birth and death terms
   !    real(rk), intent(out), optional :: birth(:)
   !       !! vector(ncells) with birth term
   !    real(rk), intent(out), optional :: death(:)
   !       !! vector(ncells) with death term

   !    real(rk) :: sumbi, weight
   !    integer:: i, j, k, n

   !    associate (gx => self%grid, nc => self%grid%ncells, array_comb => self%array_comb, &
   !               birth_ => self%birth, death_ => self%death, a => self%a)

   !       ! Evaluate afun for all particle combinations
   !       do concurrent(j=1:nc, k=1:nc)
   !          a(j, k) = self%af(gx%center(j), gx%center(k), t, y)
   !       end do

   !       ! Birh term
   !       do concurrent(i=1:nc)

   !          sumbi = 0._rk
   !          do n = 1, array_comb(i)%size()

   !             j = array_comb(i)%ia(n)
   !             k = array_comb(i)%ib(n)
   !             weight = array_comb(i)%weight(n)

   !             sumbi = sumbi + (1._rk - 0.5_rk*delta_kronecker(j, k))*weight*a(j, k)*np(j)*np(k)

   !          end do
   !          birth_(i) = sumbi

   !       end do

   !       ! Death term
   !       death_ = np*matmul(a, np)

   !       if (present(birth)) birth = birth_
   !       if (present(death)) death = death_
   !       if (present(total)) total = birth_ + death_

   !    end associate

   ! end subroutine aggterm_eval

end module break1
