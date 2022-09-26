module aggbreak1d
!!   This module contains two total variation diminishing (TVD) high-order schemes for solving
!! initial value problems. It is very important to use TVD schemes for time integration. Even
!! with a TVD spacial discretization, if the time discretization is done by a non-TVD method,
!! the result may be oscillatory.
!!   Source: ICASE 97-65 by Shu, 1997.
   use, intrinsic :: iso_fortran_env, only: real64
   use grid, only: grid1
   use auxfunctions, only: delta_kronecker
   use combtypes, only: comblist, combarray
   implicit none
   private

   public :: agg1d, agg1d_init

   integer, parameter :: rk = real64

   abstract interface
      pure real(rk) function aggkernel1d(xa, xb, t, y)
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

   pure subroutine agg1d(np, gx, aggcomb, afun, t, y, birth, death)
      real(rk), intent(in) :: np(:)
        !! vector(N) with number of particles in cell 'i'
      type(grid1), intent(in) :: gx
        !! grid object
      type(combarray), intent(in) :: aggcomb(:)
         !! vector(N) of particle combinations leading to birth
      procedure(aggkernel1d) :: afun
        !! aggregation frequency function, a(x,x',t,y)
      real(rk), intent(in) :: t
        !! time
      real(rk), intent(in) :: y(:)
        !! environment vector
      real(rk), intent(out) :: birth(:)
        !! vector(N) with birth term
      real(rk), intent(out) :: death(:)
        !! vector(N) with death term

      real(rk) :: a(size(np), size(np)), sumbi, weight
      integer:: i, j, k, n, nc

      ! Evaluate afun for all particle combinations
      nc = size(np)
      do concurrent(j=1:nc, k=1:nc)
         a(j, k) = afun(gx%center(j), gx%center(k), t, y)
      end do

      ! Birh term
      do concurrent(i=1:nc)

         sumbi = 0._rk
         do n = 1, aggcomb(i)%alength()

            j = aggcomb(i)%ia(n)
            k = aggcomb(i)%ia(n)
            weight = aggcomb(i)%weight(n)

            sumbi = sumbi + (1._rk - 0.5_rk*delta_kronecker(j, k))*weight*a(j, k)*np(j)*np(k)

         end do
         birth(i) = sumbi

      end do

      ! Death term
      death = np*matmul(a, np)

   end subroutine agg1d

   impure subroutine agg1d_init(gx, m, aggcomb)
      type(grid1), intent(in) :: gx
         !! grid object
      integer, intent(in) :: m
         !! moment of 'x' to be conserved upon aggregation
      type(combarray), intent(out) :: aggcomb(:)
         !! vector(N) of particle combinations leading to birth

      type(comblist) :: list_comb
      real(rk) :: vcenter(gx%ncells), vnew, vc, vl, vr, weight
      integer :: i, j, k, nc
      logical :: found

      ! Aux vector with 'm'-th power of cell centers
      vcenter = gx%center**m

      ! Instantiate dynamic list
      call list_comb%new

      ! Seach for all combinations of particles i,j leading to the
      ! birth of a new particle in the region of "influence" of cell k
      nc = gx%ncells
      do i = 1, nc

         ! Center of current cell
         vc = vcenter(i)

         ! Left and right cells, with protection for domain bounds
         if (i == 1) then
            vl = gx%left(i)**m
         else if (i == nc) then
            vr = gx%right(i)**m
         else
            vl = vcenter(i - 1)
            vr = vcenter(i + 1)
         end if

         do j = 1, i
            do k = j, i

               ! Coordinate of the new (born) particle
               vnew = vcenter(j) + vcenter(k)

               ! Check if new particle falls in region of "influence" of cell k
               found = .true.
               if ((vl <= vnew) .and. (vnew <= vc)) then
                  weight = (vnew - vl)/(vc - vl)
               else if ((vc <= vnew) .and. (vnew <= vr)) then
                  weight = (vr - vnew)/(vr - vc)
               else
                  found = .false.
               end if

               if (found) then
                  call list_comb%append(j, k, weight)
               end if

            end do
         end do

         ! Allocate substructure of right size
         call aggcomb(i)%alloc(list_comb%llength())

         ! Copy list content to aggcomb, then clear list
         call list_comb%toarray(aggcomb(i))
         call list_comb%clear

      end do

   end subroutine agg1d_init

end module aggbreak1d
