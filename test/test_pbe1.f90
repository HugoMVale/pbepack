module test_pbe1
!! Test for module 'pbe' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds, only: rk
   use pbepack_pbe1, only: pbe, pbesol
   use pbepack_math, only: boxcar
   use hrweno_grids, only: grid1
   use utils_tests
   use stdlib_strings, only: to_string
   use stdlib_math, only: linspace
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
                  new_unittest("PBE terms", test_pbeterms), &
                  new_unittest("PBE integration", test_pbeintegration) &
                  ]

   end subroutine

   subroutine test_pbeterms(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 42
      type(grid1) :: gx
      type(pbe) :: equation
      real(rk) :: u(nc), udot(nc, 0:4)
      real(rk) :: y(0:0)
      integer :: moment

      ! Init linear grid
      call gx%linear(1._rk, 1e3_rk, nc)

      ! Init pbe object
      moment = 1
      equation = pbe(gx, gfnc=glinear, afnc=aconst, bfnc=bconst, dfnc=dfuni, moment=moment, &
                     name="test_pbeterms")

      ! Evaluate rate of change at a given point
      u = ZERO; u(1:nc/2 - 1) = ONE
      y = ZERO
      call equation%eval(u, y, udot(:, 0))
      call equation%agg%eval(u, y, udot(:, 1))
      call equation%break%eval(u, y, udot(:, 2))
      call equation%growth%eval(u, y, udot(:, 3))
      udot(:, 4) = ZERO ! place holder for source
      call check(error, udot(:, 0), sum(udot(:, 1:), dim=2))

   contains

      pure real(rk) function dfuni(xd, xo, y) result(res)
         real(rk), intent(in) :: xd, xo, y(:)
         res = duniform(xd, xo, y, moment)
      end function

   end subroutine test_pbeterms

   subroutine test_pbeintegration(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 42/2
      type(grid1) :: gx
      type(pbe) :: equation
      type(pbesol) :: solution
      real(rk) :: u(nc), y0(1:-2) = ONE
      real(rk), allocatable :: times(:)

      ! Init linear grid
      call gx%linear(1._rk, 2e3_rk, nc)

      ! Init pbe object
      equation = pbe(gx, gfnc=glinear, afnc=aconst, bfnc=bconst, dfnc=dfuni, moment=1, &
                     name="test_pbeintegration")

      ! Integrate
      times = linspace(ZERO, 1e0_rk, 10)
      solution = equation%integrate(times, u0, y0, rtol=1e-5_rk, atol=1e-10_rk, verbose=.true.)
      !call check(error, result(:, 0), sum(result(:, 1:3), dim=2))

   contains

      pure real(rk) function dfuni(xd, xo, y) result(res)
         real(rk), intent(in) :: xd, xo, y(:)
         res = duniform(xd, xo, y, moment=1)
      end function

      pure real(rk) function u0(x) result(res)
         real(rk), intent(in) :: x
         res = boxcar(x, 1e1_rk, 5e2_rk, ONE)
      end function

   end subroutine test_pbeintegration

end module test_pbe1
