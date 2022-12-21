module test_growth1
    !! Test for module 'growth1' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds
   use pbepack_pbe1, only: pbe
   use pbepack_quadratures, only: evalmoment
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

      integer, parameter :: ncells = 10*3
      type(grid1) :: gx
      type(pbe) :: mypbe
      real(rk), dimension(ncells) :: u, udot
      real(rk) :: delta_moment(0:1), unorm
      integer :: i

      call gx%linear(1._rk, 6._rk, ncells)

      mypbe = pbe(grid=gx, g=gconst, name="test_moment_conservation")
      u = ZERO; u(ncells/3 + 1:2*ncells/3) = ONE
      call mypbe%growth%eval(u, y=VOIDREAL, udot=udot)

      unorm = evalmoment(u, gx, 0)
      delta_moment(0) = evalmoment(udot, gx, 0)
      delta_moment(1) = evalmoment(udot, gx, 1)
      delta_moment = delta_moment/unorm

      call check(error, delta_moment, [ZERO, ONE])

      if (allocated(error) .or. verbose) then
         print *
         write (stderr, '(a18,a24)') "scale type       =", gx%scale
         do i = lbound(delta_moment, 1), ubound(delta_moment, 1)
            write (stderr, '(a18,(es24.14e3))') "delta_moment("//to_string(i)//")  =", delta_moment(i)
         end do
      end if

   end subroutine test_moment_conservation

end module test_growth1
