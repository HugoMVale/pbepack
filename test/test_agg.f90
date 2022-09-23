module test_agg
!! Test for module 'weno' using test-drive.
   use iso_fortran_env, only: real64, stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use aggregation, only: agg1, agg1_init, combseries
   use grid, only: grid1
   implicit none
   private

   public :: collect_tests_agg

   integer, parameter :: rk = real64
   logical, parameter :: verbose = .false.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_agg(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("aggregation 1D", test_agg1_init) &
                  !new_unittest("calc_c", test_calc_c), &
                  !new_unittest("wenok with non-uniform grid", test_wenok_nonuniform) &
                  ]

   end subroutine

   subroutine test_agg1_init(error)
      type(error_type), allocatable, intent(out) :: error
      type(grid1) :: gx
      integer, parameter :: nc = 4
      type(combseries) :: aggcomb(nc)

      call gx%new(1._rk, 100._rk, nc)

      call agg1_init(gx, 3, aggcomb)

      ! call check(error, 1, 1)
      ! real(rk), dimension(:), allocatable :: vext, vl, vr
      ! real(rk) :: eps, atol
      ! integer :: i, k, nc = 30

      ! ! Run check for each order
      ! eps = 1e-6_rk
      ! atol = 1e-9_rk
      ! do k = 1, 3

      !    ! Allocate arrays
      !    if (allocated(vext)) deallocate (vext)
      !    if (allocated(vl)) deallocate (vl)
      !    if (allocated(vr)) deallocate (vr)
      !    allocate (vext(1 - (k - 1):nc + (k - 1)), vl(nc), vr(nc))

      !    ! Set cell average value including ghost cells
      !    ! Just a reactangular pulse __|¯¯|__
      !    vext = 0
      !    vext(nc/3:2*nc/3) = 1

      !    ! Call procedure
      !    call wenok(k, eps, vext, vl, vr)

      !    ! Check error
      !    call check(error, vext(1:nc), vl, thr=atol)
      !    call check(error, vext(1:nc), vr, thr=atol)

      !    ! Detailed comparison for debugging
      !    if (allocated(error) .or. verbose) then
      !       write (stderr, '(2(a4),3(a26))'), "k", "i", "v(i)", "vl(i)", "vr(i)"
      !       do i = 1, nc
      !          write (stderr, '(2(i4),3(es26.16e3))'), k, i, vext(i), vl(i), vr(i)
      !       end do
      !    end if

      ! end do

   end subroutine test_agg1_init

end module test_agg
