module test_agg1d
!! Test for module 'agg1d' using test-drive.
   use iso_fortran_env, only: real64, stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use agg1d, only: agg, agg_init
   use combtypes, only: combarray
   use grid, only: grid1
   implicit none
   private

   public :: collect_tests_agg1d

   integer, parameter :: rk = real64
   logical, parameter :: verbose = .false.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_agg1d(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mass conservation", test_mass_conservation) &
                  !new_unittest("calc_c", test_calc_c), &
                  !new_unittest("wenok with non-uniform grid", test_wenok_nonuniform) &
                  ]

   end subroutine

   subroutine test_mass_conservation(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 123
      type(grid1) :: gx
      type(combarray) :: aggcomb(nc)
      real(rk), dimension(nc) :: np, birth, death
      real(rk) :: t, y(0:0), t0, tend, sum_birth, sum_death
      integer :: m, scl

      call cpu_time(t0)

      ! Test linear and log grids
      do scl = 1, 2
         call gx%new(1._rk, 1e3_rk, nc, scl=scl)

         ! Test different moments
         do m = 1, 3
            call agg_init(gx, m, aggcomb)
            np = 0
            np(1:nc/2 - 1) = 1
            t = 0._rk
            y = 0
            call agg(np, gx, aggcomb, aconst, t, y, birth, death)
            sum_birth = sum(birth*gx%center**m)
            sum_death = sum(death*gx%center**m)
            call check(error, sum_birth, sum_death, rel=.true., thr=1e-14_rk)

            if (allocated(error) .or. verbose) then
               write (stderr, '(a12,es26.16e3)'), "sum_birth= ", sum_birth
               write (stderr, '(a12,es26.16e3)'), "sum_death= ", sum_death
            end if
            if (allocated(error)) return
         end do
      end do

      call cpu_time(tend)
      print '("Time = ",f8.5," seconds.")', (tend - t0)

   end subroutine test_mass_conservation

   ! subroutine test_mass_conservation(error)
   !    type(error_type), allocatable, intent(out) :: error

   !    integer, parameter :: nc = 100
   !    type(grid1) :: gx
   !    type(combarray) :: aggcomb(nc)
   !    real(rk), dimension(nc) :: np, birth, death
   !    real(rk) :: t, y(0:0), t0, tend
   !    integer :: i, m, scl

   !    do scl = 1, 2
   !       call gx%new(0._rk, 1e3_rk, nc, scl=scl)
   !       do m = 1, 3
   !          call agg_init(gx, 1, aggcomb)
   !          np = 0
   !          np(1:nc/2 - 1) = 1
   !          t = 0._rk
   !          y = 0
   !          call agg(np, gx, aggcomb, aconst, t, y, birth, death)
   !          print *, "dn/n0=", sum(birth - death)/sum(np)
   !          print *, "birth=", sum(birth*gx%center)
   !          print *, "death=", sum(death*gx%center)
   !       end do
   !    end do
   !    ! do i = 1, nc
   !    !    print *, "i= ", i
   !    !    print *, "ia=", aggcomb(i)%ia
   !    !    print *, "ib=", aggcomb(i)%ib
   !    !    print *, "we=", aggcomb(i)%weight
   !    ! end do

   !    call cpu_time(t0)
   !    call cpu_time(tend)

   !    print '("Time = ",f8.5," seconds.")', (tend - t0)

   !    !call agg1(np, gx, aggcomb, afun, t, y, birth, death)

   !    ! call check(error, 1, 1)
   !    ! real(rk), dimension(:), allocatable :: vext, vl, vr
   !    ! real(rk) :: eps, atol
   !    ! integer :: i, k, nc = 30

   !    ! ! Run check for each order
   !    ! eps = 1e-6_rk
   !    ! atol = 1e-9_rk
   !    ! do k = 1, 3

   !    !    ! Allocate arrays
   !    !    if (allocated(vext)) deallocate (vext)
   !    !    if (allocated(vl)) deallocate (vl)
   !    !    if (allocated(vr)) deallocate (vr)
   !    !    allocate (vext(1 - (k - 1):nc + (k - 1)), vl(nc), vr(nc))

   !    !    ! Set cell average value including ghost cells
   !    !    ! Just a reactangular pulse __|¯¯|__
   !    !    vext = 0
   !    !    vext(nc/3:2*nc/3) = 1

   !    !    ! Call procedure
   !    !    call wenok(k, eps, vext, vl, vr)

   !    !    ! Check error
   !    !    call check(error, vext(1:nc), vl, thr=atol)
   !    !    call check(error, vext(1:nc), vr, thr=atol)

   !    !    ! Detailed comparison for debugging
   !    !    if (allocated(error) .or. verbose) then
   !    !       write (stderr, '(2(a4),3(a26))'), "k", "i", "v(i)", "vl(i)", "vr(i)"
   !    !       do i = 1, nc
   !    !          write (stderr, '(2(i4),3(es26.16e3))'), k, i, vext(i), vl(i), vr(i)
   !    !       end do
   !    !    end if

   !    ! end do

   ! end subroutine test_1

   pure real(rk) function aconst(xa, xb, t, y) result(res)
   !! Constant aggregation kernel for 1D system
      real(rk), intent(in) :: xa
         !! internal coordinate of particle a
      real(rk), intent(in) :: xb
         !! internal coordinate of particle b
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in) :: y(:)
         !! environment vector
      res = 1._rk
   end function

end module test_agg1d
