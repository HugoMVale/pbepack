module pbepack_agg
!! Derived types and procedures to compute the *aggregation* term for 1D PBEs.
   use pbepack_kinds, only: rk, ZERO, ONE, HALF
   use pbepack_lib, only: delta_kronecker
   use pbepack_algebra, only: spmatrix
   use pbepack_basetypes, only: particleterm
   use pbepack_aggtypes, only: comblist, combarray
   use hrweno_grids, only: grid1
   implicit none
   private

   public :: aggterm, afnc_t

   type, extends(particleterm) :: aggterm
   !! Aggregation term class.
      private
      procedure(afnc_t), nopass, pointer :: afnc => null()
         !! aggregation frequency function
      type(spmatrix) :: a
         !! matrix of aggregation frequencies
      type(combarray), allocatable :: array_comb(:)
         !! array of particle combinations and weights for birth term
      logical :: update_a = .true.
         !! flag to select if matrix **a** should be updated at each step.
      logical :: empty_a = .true.
         !! flag indicating state of matrix **a**.
   contains
      procedure, pass(self), public :: eval => aggterm_eval
      procedure, pass(self), private :: compute_combinations
      procedure, pass(self), private :: compute_a
   end type aggterm

   abstract interface
      pure real(rk) function afnc_t(xa, xb, y)
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

   type(aggterm) function aggterm_init(grid, a, moment, update_a, name) result(self)
   !! Initialize `aggterm` object.
   !!$$
   !! \left ( \dot{u} \right )_{\mathrm{aggregation}} =
   !! \int _0^{x/2} a(x-x',x',\mathbf{y})u(x-x',t)u(x',t)dx'
   !! -u(x,t)\int _0^\infty a(x,x',\mathbf{y})u(x',t)dx'
   !!$$
   !! where \(moment\) is the moment of \(x\) conserved upon aggregation. For example,
   !! if \(x\) denotes particle mass or volume, then \(moment=1\), whereas if \(x\) denotes
   !! particle radius or diameter, then \(moment=3\).
      type(grid1), intent(in) :: grid
         !! `grid1` object
      procedure(afnc_t) :: a
         !! aggregation frequency function, \( a(x,x',\textbf{y}) \)
      integer, intent(in), optional :: moment
         !! moment of \(x\) conserved during aggregation (default=1)
      logical, intent(in), optional :: update_a
         !! if `true`, \( a(x,x',\textbf{y}) \) is reevaluated at each step (default=true)
      character(*), intent(in), optional :: name
         !! name (default="agg-term")

      ! Set parameters
      call self%set_name(name, "agg-term")
      call self%set_grid(grid)
      self%afnc => a
      if (present(moment)) call self%set_moment(moment)
      if (present(update_a)) self%update_a = update_a

      ! Allocate and pre-compute stuff
      call self%particleterm_allocations()
      associate (nc => self%grid%ncells)
         allocate (self%array_comb(nc))
         self%a = spmatrix(nc)
      end associate
      call self%compute_combinations()
      self%inited = .true.

   end function aggterm_init

   pure subroutine aggterm_eval(self, u, y, udot, udot_birth, udot_death)
   !! Evaluate rate of aggregation at a given instant \(t\).
      class(aggterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: u(:)
         !! cell-average number density, \( \bar{u_i}(t) \)
      real(rk), intent(in) :: y(:)
         !! environment vector, \(y_j(t)\)
      real(rk), intent(out), optional :: udot(:)
         !! net rate of change (birth-death), \( \dot{u_i}(t) \)
      real(rk), intent(out), optional :: udot_birth(:)
         !! rate of birth
      real(rk), intent(out), optional :: udot_death(:)
         !! rate of death

      real(rk) :: weight, up(size(u))
      integer:: i, j, k, n

      call self%check_inited()

      associate (nc => self%grid%ncells, &
                 array_comb => self%array_comb, a => self%a, &
                 birth => self%udot_birth, death => self%udot_death)

         ! Evaluate aggregation frequency for all particle combinations
         if (self%update_a .or. self%empty_a) call self%compute_a(y)

         ! Number of particles in each cell
         up = u*self%grid%width

         ! Birth term
         birth = ZERO
         do concurrent(i=1:nc)
            do n = 1, array_comb(i)%size()
               j = array_comb(i)%ia(n)
               k = array_comb(i)%ib(n)
               weight = array_comb(i)%weight(n)
               birth(i) = birth(i) + &
                          (ONE - HALF*delta_kronecker(j, k))*weight*a%get(j, k)*up(j)*up(k)
            end do
         end do
         birth = birth/self%grid%width
         if (present(udot_birth)) udot_birth = birth

         ! Death term
         death = u*a%multvec(up)
         if (present(udot_death)) udot_death = death

         ! Net rate
         self%udot = birth - death
         if (present(udot)) udot = self%udot

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
      associate (nc => self%grid%ncells, x => self%grid%center)
         do concurrent(j=1:nc)
            do concurrent(i=1:j)
               call self%a%set(i, j, self%afnc(x(i), x(j), y))
            end do
         end do
      end associate

      if (self%empty_a) self%empty_a = .false.

   end subroutine compute_a

   impure subroutine compute_combinations(self)
   !! Precompute particle combinations leading to birth and respective weights.
   !! @note This procedure is impure because some 'ftlList' methods are impure. Otherwise,
   !! this procedure could be "pure".
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

end module pbepack_agg
