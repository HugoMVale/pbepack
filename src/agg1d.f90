module agg1d
!!   This module contains ...
   use, intrinsic :: iso_fortran_env, only: real64
   use grid, only: grid1
   use auxfunctions, only: delta_kronecker
   use combtypes, only: comblist, combarray
   use pbetools, only: pbeterm
   implicit none
   private

   public :: aggterm

   integer, parameter :: rk = real64

   type, extends(pbeterm) :: aggterm
      integer :: moment
      type(combarray), allocatable :: array_comb(:)
   contains
      procedure, pass(self) :: init => aggterm_init
      procedure, pass(self) :: eval => aggterm_eval
      procedure, pass(self), private :: aggterm_combinations
   end type

   abstract interface
      pure real(rk) function aggfun1d(xa, xb, t, y)
      !! aggregation kernel for 1D system
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

contains

   impure subroutine aggterm_init(self, grid, moment)
   !! Initialize 'aggterm' object.
      class(aggterm), intent(inout) :: self
         !! object
      type(grid1), intent(inout), target :: grid
         !! grid
      integer, intent(in) :: moment
         !! moment of 'x' to be conserved upon aggregation

      ! Clear object if required
      if (allocated(self%name)) deallocate (self%name)
      if (allocated(self%msg)) deallocate (self%msg)
      if (allocated(self%array_comb)) deallocate (self%array_comb)
      if (associated(self%grid)) nullify (self%grid)
      self%ierr = 1

      ! Assign inputs
      self%grid => grid
      self%moment = moment
      self%ncells = self%grid%ncells

      ! Precompute combinations
      allocate (self%array_comb(self%ncells))
      call self%aggterm_combinations()
      self%ierr = 0

   end subroutine aggterm_init

   pure subroutine aggterm_eval(self, np, afun, t, y, birth, death)
   !! Evaluate rate of aggregation at given instant
      class(aggterm), intent(inout) :: self
         !! object
      real(rk), intent(in) :: np(:)
         !! vector(N) with number of particles in cell 'i'
      procedure(aggfun1d) :: afun
         !! aggregation frequency function, a(x,x',t,y)
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), intent(out) :: birth(:)
         !! vector(N) with birth term
      real(rk), intent(out) :: death(:)
         !! vector(N) with death term

      real(rk) :: a(self%ncells, self%ncells), sumbi, weight
      integer:: i, j, k, n

      associate (gx => self%grid, nc => self%ncells, array_comb => self%array_comb)

         ! Evaluate afun for all particle combinations
         do concurrent(j=1:nc, k=1:nc)
            a(j, k) = afun(gx%center(j), gx%center(k), t, y)
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
            birth(i) = sumbi

         end do

         ! Death term
         death = np*matmul(a, np)

      end associate

   end subroutine aggterm_eval

   impure subroutine aggterm_combinations(self)
   !! Precompute particle combinations leading to birth and respective weights.
      class(aggterm), intent(inout) :: self
         !! object

      type(comblist) :: list_comb
      real(rk) :: vcenter(self%ncells), vnew, vc, vl, vr, weight
      integer :: i, j, k
      logical :: found

      associate (gx => self%grid, nc => self%ncells, m => self%moment)
         ! Aux vector with 'm'-th power of cell centers
         vcenter = gx%center**m

         ! Instantiate dynamic list
         call list_comb%new

         ! Seach for all combinations of particles i,j leading to the
         ! birth of a new particle in the region of "influence" of cell k
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

end module agg1d
