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
      procedure(af1_t), nopass, pointer :: af => null()
         !! aggregation frequency
      real(rk), allocatable :: a(:, :)
         !! matrix of aggregation frequencies
      type(combarray), allocatable, private :: array_comb(:)
         !! Array of particle combinations and weights for birth term
   contains
      procedure, pass(self) :: eval => aggterm_eval
      procedure, pass(self), private :: aggterm_combinations
   end type aggterm

   abstract interface
      pure real(rk) function af1_t(xa, xb, y)
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

   type(aggterm) function aggterm_init(af, moment, grid, name) result(self)
   !! Initialize 'aggterm' object.
      procedure(af1_t) :: af
         !! aggregation frequency, \( a(x,x',y) \)
      integer, intent(in) :: moment
         !! moment of 'x' to be conserved upon aggregation
      type(grid1), intent(in), optional :: grid
         !! grid1 object
      character(*), intent(in), optional :: name
         !! object name

      self%af => af

      call self%set_moment(moment)

      if (present(grid)) then
         call self%set_grid(grid)
         call self%aggterm_combinations()
         self%inited = .true.
         self%msg = "Initialization completed successfully"
      else
         self%msg = "Missing 'grid'."
      end if

      self%name = optval(name, "")

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
         allocate (self%array_comb(nc), self%source(nc), self%sink(nc), self%a(nc, nc))

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

            ! Copy list content to array, then clear list
            call list_comb%toarray(self%array_comb(i))
            call list_comb%clear

         end do
      end associate

   end subroutine aggterm_combinations

   pure subroutine aggterm_eval(self, np, y, total, source, sink)
   !! Evaluate rate of aggregation at a given instant.
      class(aggterm), intent(inout) :: self
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

      real(rk) :: weight
      integer:: i, j, k, n

      associate (gx => self%grid, nc => self%grid%ncells, array_comb => self%array_comb, &
                 source_ => self%source, sink_ => self%sink, a => self%a)

         ! Evaluate aggregation frequency for all particle combinations
         do concurrent(j=1:nc, k=1:nc)
            a(j, k) = self%af(gx%center(j), gx%center(k), y)
         end do

         ! Source/birth term
         source_ = ZERO
         do concurrent(i=1:nc)
            do n = 1, array_comb(i)%size()
               j = array_comb(i)%ia(n)
               k = array_comb(i)%ib(n)
               weight = array_comb(i)%weight(n)
               source_(i) = source_(i) + &
                            (ONE - HALF*delta_kronecker(j, k))*weight*a(j, k)*np(j)*np(k)
            end do
         end do

         ! Sink/death term
         sink_ = np*matmul(a, np)

         if (present(source)) source = source_
         if (present(sink)) sink = sink_
         if (present(total)) total = source_ + sink_

      end associate

   end subroutine aggterm_eval

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
