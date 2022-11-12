module test_pbe1
!! Test for module 'pbe' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds, only: rk
   use pbepack_agg1, only: aggterm
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
                  new_unittest("PBE with aggterm", test_agg) &
                  !new_unittest("analytical solution, case 1", test_case1) &
                  !new_unittest("wenok with non-uniform grid", test_wenok_nonuniform) &
                  ]

   end subroutine

   subroutine test_agg(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 100
      type(grid1) :: gx
      type(aggterm) :: agg
      type(pbe1) :: eq
      real(rk), dimension(nc) :: np, result
      real(rk) :: y(0:0), moment_birth_0, moment_birth_m, moment_death_0, moment_death_m
      integer :: moment, scl

      ! Init linear grid
      call gx%linear(1._rk, 1e3_rk, nc)

      ! Init aggregation term
      agg = aggterm(afnc=aconst, moment=1, grid=gx, name="agg")

      ! Init pbe object
      eq = pbe1(grid=gx, agg=agg, name="a pbe with aggregation only")

      ! Evaluate pbe at a given point
      np = 0
      np(1:nc/2 - 1) = 1
      y = 0._rk
      call eq%eval(np, y, result)

   end subroutine test_agg

end module test_pbe1
