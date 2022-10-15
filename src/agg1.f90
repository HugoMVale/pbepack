module agg1
   use real_kinds, only: rk
   use basetypes, only: pbeterm
   use grid, only: grid1
   use combtypes, only: comblist, combarray
   use auxfunctions, only: delta_kronecker
   use stdlib_optval, only: optval
   use omp_lib

   implicit none
   private

   public :: aggterm

   type, extends(pbeterm) :: aggterm
   !! Aggregation term class.
      procedure(aggf1), nopass, pointer :: af => null()
         !! aggregation frequency
      integer :: moment
         !! moment of 'x' to be conserved upon aggregation
      real(rk), allocatable :: birth(:)
         !! birth term
      real(rk), allocatable :: death(:)
         !! death term
      real(rk), allocatable :: a(:, :)
         !! array of aggregation frequencies
      type(combarray), allocatable, private :: array_comb(:)
         !! Array of particle combinations and weights for birth term
   contains
      procedure, pass(self) :: eval => aggterm_eval
      procedure, pass(self), private :: aggterm_combinations
   end type aggterm

   abstract interface
      pure real(rk) function aggf1(xa, xb, t, y)
      !! Aggregation kernel for 1D system
         import :: rk
         real(rk), intent(in) :: xa
            !! internal coordinate of particle a
         real(rk), intent(in) :: xb
            !! internal coordinate of particle b
         real(rk), intent(in) :: t
            !! time
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function
   end interface

   interface aggterm
      module procedure :: aggterm_init
   end interface aggterm

contains

   type(aggterm) function aggterm_init(af, moment, grid) result(self)
   !! Initialize 'aggterm' object.
      procedure(aggf1), optional :: af
         !! aggregation frequency
      integer, intent(in) :: moment
         !! moment of 'x' to be conserved upon aggregation
      type(grid1), intent(in), target, optional :: grid
         !! grid1 object

      ! Assign inputs
      self%af => af
      if (moment > 0) then
         self%moment = moment
      else
         self%msg = "Invalid input 'moment'. Valid range: moment >= 1."
         self%ierr = 1
         error stop self%msg
      end if

      if (present(grid)) then
         self%grid => grid
         call self%aggterm_combinations()
         self%inited = .true.
         self%msg = "Initialization completed successfully"
      else
         self%msg = "Missing 'grid'."
      end if

      ! Precompute combinations

   end function aggterm_init

   impure subroutine aggterm_combinations(self)
   !! Precompute particle combinations leading to birth and respective weights.
      class(aggterm), intent(inout) :: self
         !! object

      type(comblist) :: list_comb
      real(rk) :: vcenter(self%grid%ncells), vnew, vc, vl, vr, weight
      integer :: i, j, k
      logical :: found

      associate (gx => self%grid, nc => self%grid%ncells, m => self%moment)

         ! Allocate internal storage arrays
         allocate (self%array_comb(nc), self%birth(nc), self%death(nc), self%a(nc, nc))

         ! Aux vector with 'm'-th power of cell centers
         vcenter = gx%center**m

         ! Seach for all combinations of particles i,j leading to the
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

            ! Copy list content to aggcomb, then clear list
            call list_comb%toarray(self%array_comb(i))
            call list_comb%clear

         end do
      end associate

   end subroutine aggterm_combinations

   pure subroutine aggterm_eval(self, np, t, y, total, birth, death)
   !! Evaluate rate of aggregation at a given instant.
      class(aggterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: np(:)
         !! vector(ncells) with number of particles in cell 'i'
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), intent(out), optional :: total(:)
         !! vector(ncells) with sum of birth and death terms
      real(rk), intent(out), optional :: birth(:)
         !! vector(ncells) with birth term
      real(rk), intent(out), optional :: death(:)
         !! vector(ncells) with death term

      real(rk) :: sumbi, weight
      integer:: i, j, k, n

      associate (gx => self%grid, nc => self%grid%ncells, array_comb => self%array_comb, &
                 birth_ => self%birth, death_ => self%death, a => self%a)

         ! Evaluate afun for all particle combinations
         do concurrent(j=1:nc, k=1:nc)
            a(j, k) = self%af(gx%center(j), gx%center(k), t, y)
         end do

         ! Birh term
         do concurrent(i=1:nc)

            sumbi = 0._rk
            do n = 1, array_comb(i)%size()

               j = array_comb(i)%ia(n)
               k = array_comb(i)%ib(n)
               weight = array_comb(i)%weight(n)

               sumbi = sumbi + (1._rk - 0.5_rk*delta_kronecker(j, k))*weight*a(j, k)*np(j)*np(k)

            end do
            birth_(i) = sumbi

         end do

         ! Death term
         death_ = np*matmul(a, np)

         if (present(birth)) birth = birth_
         if (present(death)) death = death_
         if (present(total)) total = birth_ + death_

      end associate

   end subroutine aggterm_eval

end module agg1
