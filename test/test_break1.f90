module test_break1
!! Test for module 'break1d' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds
   use pbepack_lib
   use pbepack
   use utils_tests
   implicit none
   private

   public :: collect_tests_break1

   logical, parameter :: verbose = .false.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_break1(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("moment conservation", test_moment_conservation), &
                  new_unittest("case b(x)=x", test_blinear) &
                  ]

   end subroutine

   subroutine test_moment_conservation(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncells = 200
      type(grid1) :: gx
      type(pbe) :: mypbe
      real(rk), dimension(ncells) :: u, udot_birth, udot_death
      real(rk) :: moment_birth_0, moment_birth_m, moment_death_0, moment_death_m
      integer :: moment, scale

      ! Test linear and log grids
      do scale = 1, 2
         select case (scale)
         case (1)
            call gx%linear(0._rk, 3e2_rk, ncells)
         case (2)
            call gx%geometric(0._rk, 2e3_rk, 1.01_rk, ncells)
         end select

         ! Test different moments
         do moment = 1, 3

            mypbe = pbe(grid=gx, b=bconst, d=dfuni, moment=moment, update_b=.false., &
                        name="test_moment_conservation")
            u = ZERO; u(ncells) = ONE
            call mypbe%break%eval(u, y=VOIDREAL, udot_birth=udot_birth, udot_death=udot_death)

            moment_birth_0 = evalmoment(udot_birth, gx, 0)
            moment_death_0 = evalmoment(udot_death, gx, 0)
            moment_birth_m = evalmoment(udot_birth, gx, moment)
            moment_death_m = evalmoment(udot_death, gx, moment)

            call check(error, moment_birth_0, moment_death_0*2, rel=.true., thr=1e-2_rk)
            call check(error, moment_birth_m, moment_death_m, rel=.true., thr=1e-8_rk)

            if (allocated(error) .or. verbose) then
               print *
               write (stderr, '(a18,a24)') "scale type       =", gx%scale
               write (stderr, '(a18,i24)') "preserved moment =", moment
               write (stderr, '(a18,(es24.14e3))') "source/sink(0)   =", moment_birth_0/moment_death_0
               write (stderr, '(a18,(es24.14e3))') "source/sink("//to_string(moment)//")   =", moment_birth_m/moment_death_m
            end if

            if (allocated(error)) return

         end do
      end do

   contains

      pure real(rk) function dfuni(xd, xo, y) result(res)
         real(rk), intent(in) :: xd, xo, y(:)
         res = duniform(xd, xo, y, moment)
      end function

   end subroutine test_moment_conservation

   subroutine test_blinear(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 200
      type(grid1) :: gx
      type(pbe) :: equation
      type(pbesol) :: solution
      real(rk) :: u(nc), n0, x0
      real(rk), dimension(0:2) :: moment, moment_ref
      real(rk), allocatable :: times(:)
      integer :: i

      ! Init linear grid
      call gx%geometric(ZERO, 2e1_rk, 1.03_rk, nc)

      ! Init pbe object
      equation = pbe(gx, b=blinear, d=dfuni, name="test_blinear")

      ! Integrate
      times = linspace(ZERO, TWO, 2)
      solution = equation%integrate(times, u0)

      ! Compute solution moments
      do i = 0, 2
         moment(i) = evalmoment(solution%u(:, size(times)), gx, i)
      end do

      ! Compute reference moments
      n0 = evalmoment(solution%u(:, 1), gx, 0)
      x0 = evalmoment(solution%u(:, 1), gx, 1)/n0
      moment_ref = blinear_moments(times(size(times)), x0=x0, n0=n0)

      ! Check moments
      call check(error, moment(0), moment_ref(0), rel=.true., thr=1e-2_rk)
      call check(error, moment(1), moment_ref(1), rel=.true., thr=1e-5_rk)

      if (allocated(error) .or. verbose) then
         do i = 0, 2
            write (stderr, '(a18,(es24.14e3))') &
               "num./analyt.("//to_string(i)//") =", moment(i)/moment_ref(i)
         end do
         call solution%write("test_blinear")
      end if

   contains

      pure real(rk) function dfuni(xd, xo, y) result(res)
         real(rk), intent(in) :: xd, xo, y(:)
         res = duniform(xd, xo, y, 1)
      end function

   end subroutine test_blinear

   pure real(rk) function u0(x) result(res)
      real(rk), intent(in) :: x
      res = expo1d(x)
   end function

end module test_break1
