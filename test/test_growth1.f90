module test_growth1
    !! Test for module 'growth1' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds
   use pbepack_pbe1, only: pbe1
   use hrweno_grids, only: grid1
   use utils_tests, only: gconst
   use stdlib_strings, only: to_string
   implicit none
   private

   public :: collect_tests_growth1

   logical, parameter :: verbose = .true.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_growth1(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("moment conservation", test_moment_conservation) &
                  !new_unittest("analytical solution, case 1", test_case1) &
                  ]

   end subroutine

   subroutine test_moment_conservation(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 10*3
      type(grid1) :: gx
      type(pbe1) :: eq
      real(rk), dimension(nc) :: u, udot
      real(rk) :: y(0:0), delta_moment(0:1)
      integer :: i

      call gx%linear(1._rk, 6._rk, nc)

      eq = pbe1(grid=gx, gfnc=gconst, name="test_moment_conservation")
      u = ZERO; u(nc/3 + 1:2*nc/3) = ONE
      y = ZERO
      call eq%growth%eval(u, y, udot=udot)

      do i = lbound(delta_moment, 1), ubound(delta_moment, 1)
         delta_moment(i) = sum(udot*gx%width*gx%center**i)/sum(u*gx%width)
      end do

      call check(error, delta_moment(0), ZERO)
      call check(error, delta_moment(1), ONE, rel=.true., thr=1e-14_rk)

      if (allocated(error) .or. verbose) then
         print *
         write (stderr, '(a18,a24)') "scale type       =", gx%scale
         do i = lbound(delta_moment, 1), ubound(delta_moment, 1)
            write (stderr, '(a18,(es24.14e3))') "delta_moment("//to_string(i)//")  =", delta_moment(i)
         end do
      end if

   end subroutine test_moment_conservation

end module test_growth1
