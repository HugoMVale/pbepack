module test_pbe1
!! Test for module 'pbe' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds, only: rk
   use pbepack_agg1, only: aggterm
   use pbepack_break1, only: breakterm
   use pbepack_pbe1, only: pbe1
   use hrweno_grids, only: grid1
   use utils_tests
   use stdlib_strings, only: to_string
   implicit none
   private

   public :: collect_tests_pbe1

   logical, parameter :: verbose = .true.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_pbe1(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("PBE terms", test_pbeterms) &
                  !new_unittest("analytical solution, case 1", test_case1) &
                  !new_unittest("wenok with non-uniform grid", test_wenok_nonuniform) &
                  ]

   end subroutine

   subroutine test_pbeterms(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 42
      type(grid1) :: gx
      type(aggterm) :: agg
      type(breakterm) :: break
      type(pbe1) :: eq
      real(rk) :: np(nc), result(nc, 0:3)
      real(rk) :: y(0:0)
      integer :: moment, i

      ! Init linear grid
      call gx%linear(1._rk, 1e3_rk, nc)

      ! Init aggregation term
      moment = 1
      agg = aggterm(afnc=aconst, moment=moment, grid=gx, name="agg")

      ! Init breakage term
      break = breakterm(bfnc=bconst, dfnc=dfuni, moment=moment, grid=gx, name="break")

      ! Init pbe object
      eq = pbe1(grid=gx, agg=agg, break=break, name="a pbe with some terms")

      ! Evaluate pbe at a given point
      np = ZERO
      np(1:nc/2 - 1) = ONE
      y = ZERO
      result = ZERO
      call eq%eval(np, y, result(:, 0))
      call agg%eval(np, y, result(:, 1))
      call break%eval(np, y, result(:, 2))
      call check(error, result(:, 0), sum(result(:, 1:), dim=2))

   contains

      pure real(rk) function dfuni(xd, xo, y) result(res)
         real(rk), intent(in) :: xd, xo, y(:)
         res = duniform(xd, xo, y, moment)
      end function

   end subroutine test_pbeterms

end module test_pbe1
