module test_pbe1
!! Test for module 'pbe' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds, only: rk
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
      type(pbe1) :: eq
      real(rk) :: np(nc), result(nc, 0:3)
      real(rk) :: y(0:0)
      integer :: moment

      ! Init linear grid
      call gx%linear(1._rk, 1e3_rk, nc)

      ! Init pbe object
      moment = 1

      eq = pbe1(gx, gfnc=glinear, afnc=aconst, bfnc=bconst, dfnc=dfuni, moment=moment, &
                name="a pbe with several terms")

      ! Evaluate pbe at a given point
      np = ZERO
      np(1:nc/2 - 1) = ONE
      y = ZERO
      result = ZERO
      call eq%agg%eval(np, y, result(:, 1))
      call eq%break%eval(np, y, result(:, 2))
      call eq%growth%eval(np, y, result(:, 3))
      call eq%eval(np, y, result(:, 0))
      call check(error, result(:, 0), sum(result(:, 1:3), dim=2))

   contains

      pure real(rk) function dfuni(xd, xo, y) result(res)
         real(rk), intent(in) :: xd, xo, y(:)
         res = duniform(xd, xo, y, moment)
      end function

   end subroutine test_pbeterms

end module test_pbe1
