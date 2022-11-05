module agg1
!! This module implements derived types and procedures to compute the aggregation term
!! for 1D PBEs.
   use real_kinds
   use basetypes, only: particleterm
   use grids, only: grid1
   use aggtypes, only: comblist, combarray
   use stdlib_optval, only: optval

   implicit none
   private

   public :: aggterm

   type, extends(particleterm) :: aggterm
   !! Aggregation term class.
      procedure(afnc1_t), nopass, pointer :: afnc => null()
         !! aggregation frequency
      real(rk), allocatable :: a(:, :)
         !! matrix of aggregation frequencies
      type(combarray), allocatable, private :: array_comb(:)
         !! array of particle combinations and weights for birth term
      logical :: update_a = .true.
         !! flag to select if matrix 'a' should be updated at each step.
      logical :: empty_a = .true.
         !! flag indicating state of matrix 'a'.
   contains
      procedure, pass(self) :: eval => aggterm_eval
      procedure, pass(self), private :: aggterm_allocations
      procedure, pass(self), private :: compute_combinations
      procedure, pass(self), private :: compute_a
   end type aggterm

   abstract interface
      pure real(rk) function afnc1_t(xa, xb, y)
      !! Aggregation frequency for 1D system
         import :: rk
         real(rk), intent(in) :: xa
            !! internal coordinate of particle a
         real(rk), intent(in) :: xb
            !! internal coordinate of particle b
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function
   end interface

   interface aggterm
      module procedure :: aggterm_init
   end interface aggterm

contains

   type(aggterm) function aggterm_init(afnc, moment, grid, update_a, name) result(self)
   !! Initialize 'aggterm' object.
      procedure(afnc1_t) :: afnc
         !! aggregation frequency, \( a(x,x',y) \)
      integer, intent(in) :: moment
         !! moment of \( x \) to be conserved upon aggregation
      type(grid1), intent(in), optional :: grid
         !! grid1 object
      logical, intent(in), optional :: update_a
         !! flag to select if matrix 'a' should be updated at each step
      character(*), intent(in), optional :: name
         !! object name

      self%afnc => afnc

      call self%set_moment(moment)

      if (present(update_a)) self%update_a = update_a

      self%name = optval(name, "")

      if (present(grid)) then
         call self%set_grid(grid)
         call self%aggterm_allocations()
         call self%compute_combinations()
         self%inited = .true.
         self%msg = "Initialization completed successfully"
      else
         self%msg = "Missing 'grid'."
      end if

   end function aggterm_init

   pure subroutine aggterm_eval(self, np, y, result, birth, death)
   !! Evaluate rate of aggregation at a given instant.
      class(aggterm), intent(inout) :: self
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

      real(rk) :: weight
      integer:: i, j, k, n

      associate (nc => self%grid%ncells, &
                 array_comb => self%array_comb, a => self%a, &
                 birth_ => self%birth, death_ => self%death, result_ => self%result)

         ! Evaluate aggregation frequency for all particle combinations
         if (self%update_a .or. self%empty_a) call self%compute_a(y)

         ! Birth term
         birth_ = ZERO
         do concurrent(i=1:nc)
            do n = 1, array_comb(i)%size()
               j = array_comb(i)%ia(n)
               k = array_comb(i)%ib(n)
               weight = array_comb(i)%weight(n)
               birth_(i) = birth_(i) + &
                           (ONE - HALF*delta_kronecker(j, k))*weight*a(j, k)*np(j)*np(k)
            end do
         end do
         if (present(birth)) birth = birth_

         ! Death term
         death_ = np*matmul(a, np)
         if (present(death)) death = death_

         ! Net rate
         result_ = birth_ - death_
         if (present(result)) result = result_

      end associate

   end subroutine aggterm_eval

   pure subroutine compute_a(self, y)
   !! Compute/update array of aggregation frequencies.
      class(aggterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: y(:)
         !! environment vector

      integer :: i, j

      ! The array is symmetric
      ! Can be improved with symmetric matrix from Lapack
      associate (nc => self%grid%ncells, x => self%grid%center)
         do j = 1, nc
            do i = j, nc
               self%a(i, j) = self%afnc(x(i), x(j), y)
            end do
         end do
         do j = 2, nc
            do i = 1, j - 1
               self%a(i, j) = self%a(j, i)
            end do
         end do
      end associate

      if (self%empty_a) self%empty_a = .false.

   end subroutine compute_a

   impure subroutine compute_combinations(self)
   !! Precompute particle combinations leading to birth and respective weights.
      class(aggterm), intent(inout) :: self
         !! object

      type(comblist) :: list_comb
      real(rk) :: vcenter(self%grid%ncells), vnew, vc, vl, vr, weight
      integer :: i, j, k
      logical :: found

      associate (gx => self%grid, nc => self%grid%ncells, m => self%moment)

         ! Aux vector with 'm'-th power of cell centers
         vcenter = gx%center**m

         ! Search for all combinations of particles i,j leading to the
         ! birth of a new particle in the region of "influence" of cell k
         call list_comb%new
         do i = 1, nc

            ! Center of current cell
            vc = vcenter(i)

            ! Left and right cells, with protection for domain bounds
            if (i == 1) then
               vl = gx%left(i)**m
            else
               vl = vcenter(i - 1)
            end if

            if (i == nc) then
               vr = gx%right(i)**m
            else
               vr = vcenter(i + 1)
            end if

            do j = 1, i
               do k = j, i

                  ! Coordinate of the new (born) particle
                  vnew = vcenter(j) + vcenter(k)

                  ! Check if new particle falls in region of "influence" of cell k
                  if ((vl < vnew) .and. (vnew <= vc)) then
                     found = .true.
                     weight = (vnew - vl)/(vc - vl)
                  else if ((vc < vnew) .and. (vnew < vr)) then
                     found = .true.
                     weight = (vr - vnew)/(vr - vc)
                  else
                     found = .false.
                  end if

                  if (found) then
                     call list_comb%append(j, k, weight)
                  end if

               end do
            end do

            ! Allocate substructure of correct size
            call self%array_comb(i)%alloc(list_comb%size())

            ! Copy list content to array, then clear list
            call list_comb%toarray(self%array_comb(i))
            call list_comb%clear

         end do
      end associate

   end subroutine compute_combinations

   pure subroutine aggterm_allocations(self)
   !! Allocate arrays.
      class(aggterm), intent(inout) :: self
         !! object

      ! Call parent method
      call self%particleterm_allocations()

      ! Do own allocations
      if (associated(self%grid)) then
         associate (nc => self%grid%ncells)
            allocate (self%array_comb(nc), self%a(nc, nc))
         end associate
      else
         self%msg = "Allocation failed due to missing grid."
         self%ierr = 1
         error stop self%msg
      end if

   end subroutine aggterm_allocations

   pure real(rk) function delta_kronecker(i, j)
   !! Delta kronecker \( \delta_{i,j} \).
      integer, intent(in) :: i
       !! integer i
      integer, intent(in) :: j
       !! integer j

      if (i == j) then
         delta_kronecker = ONE
      else
         delta_kronecker = ZERO
      end if

   end function delta_kronecker

end module agg1
