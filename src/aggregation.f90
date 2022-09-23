module aggregation
!!   This module contains two total variation diminishing (TVD) high-order schemes for solving
!! initial value problems. It is very important to use TVD schemes for time integration. Even
!! with a TVD spacial discretization, if the time discretization is done by a non-TVD method,
!! the result may be oscillatory.
!!   Source: ICASE 97-65 by Shu, 1997.
   use, intrinsic :: iso_fortran_env, only: real64
   use ftlListIntModule
   use grid, only: grid1
   implicit none
   private

   public :: agg1, agg1_init, combseries

   integer, parameter :: rk = real64

   abstract interface
      pure real(rk) function aggkernel1(xa, xb, t, y)
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

   abstract interface
      pure real(rk) function aggkernel2(xa, xb, t, y)
        !! aggregation kernel for 2D system
         import :: rk
         real(rk), intent(in) :: xa(2)
            !! internal coordinates of particle a
         real(rk), intent(in) :: xb(2)
            !! internal coordinates of particle b
         real(rk), intent(in) :: t
            !! time
         real(rk), intent(in) :: y(:)
            !! environment vector
      end function
   end interface

   type :: combseries
      integer, allocatable :: ia(:)
      integer, allocatable :: ib(:)
      real(rk), allocatable :: weight(:)
   end type

contains

   pure subroutine agg1(np, gx, rate, t, y, birth, death)
      real(rk), intent(in) :: np(:)
        !! vector(N) with number of particles in cell 'i'
      type(grid1), intent(in) :: gx
        !! grid object
      procedure(aggkernel1) :: rate
        !! aggregation rate coefficient
      real(rk), intent(in) :: t
        !! time
      real(rk), intent(in) :: y(:)
        !! environment vector
      real(rk), intent(out) :: birth(:)
        !! vector(N) with birth term
      real(rk), intent(out) :: death(:)
        !! vector(N) with death term

      real(rk) :: beta(size(np), size(np))
      integer:: i, j, nc

      ! Evaluate rate for all particle combinations
      nc = size(np)
      do concurrent(i=1:nc, j=1:nc)
         beta(i, j) = rate(gx%center(i), gx%center(j), t, y)
      end do

      ! Birh term
      birth = 0

      ! Death term
      death = np*matmul(beta, np)

   end subroutine agg1

   impure subroutine agg1_init(gx, m, aggcomb)
      type(grid1), intent(in) :: gx
            !! grid object
      integer, intent(in) :: m
            !! moment of 'x' to be conserved upon aggregation
      type(combseries), intent(out) :: aggcomb(:)
            !! vector(N) of particle combinations leading to birth
      type(ftlListInt) :: intlist(2)
      type(ftlListIntIterator) :: iter(2)

      real(rk) :: xnew, xl, xr
      integer :: i, j, k, n, nc

      ! Instantiate lists
      call intlist(1)%New
      call intlist(2)%New

      ! Seach for particles born in the region of "influence" of cell k
      nc = gx%ncells
      do k = 1, nc
         n = 0
         do i = 1, k
            do j = 1, k
               xnew = (gx%center(i)**m + gx%center(j)**m)**(1.0_rk/m)

               if (k == 1) then
                  xl = gx%left(1)
               else
                  xl = gx%center(k - 1)
               end if

               if (k == nc) then
                  xr = gx%right(k)
               else
                  xr = gx%center(k + 1)
               end if

               if ((xl <= xnew) .and. (xnew <= xr)) then
                  n = n + 1
                  call intlist(1)%PushBack(i)
                  call intlist(2)%PushBack(j)
               end if

            end do
         end do

         ! Allocate substructure of right size
         allocate (aggcomb(k)%ia(intlist(1)%Size()))
         allocate (aggcomb(k)%ib(intlist(2)%Size()))

         ! Copy list content to aggcomb
         call intlist(1)%CopyToArray(aggcomb(k)%ia)
         call intlist(2)%CopyToArray(aggcomb(k)%ib)

         ! Clear list before next iteration
         call intlist(1)%Clear
         call intlist(2)%Clear

      end do

      ! Compue the weights for each particle

   end subroutine agg1_init

end module aggregation
