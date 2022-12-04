module test_break1
!! Test for module 'break1d' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds
   use pbepack_break1, only: breakterm
   use hrweno_grids, only: grid1
   use utils_tests, only: bconst, duniform
   use stdlib_strings, only: to_string
   implicit none
   private

   public :: collect_tests_break1

   logical, parameter :: verbose = .true.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_break1(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("moment conservation", test_moment_conservation) &
                  ! new_unittest("call from pbe", test_pbe_break)
                  !new_unittest("analytical solution, case 1", test_case1) &
                  ]

   end subroutine

   subroutine test_moment_conservation(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 200
      type(grid1) :: gx
      type(breakterm) :: break
      real(rk), dimension(nc) :: np, birth, death
      real(rk) :: y(0:0), moment_birth_0, moment_birth_m, moment_death_0, moment_death_m
      integer :: moment, scale

      ! Test linear and log grids
      do scale = 1, 2
         select case (scale)
         case (1)
            call gx%linear(0._rk, 3e2_rk, nc)
         case (2)
            call gx%geometric(0._rk, 2e3_rk, 1.01_rk, nc)
         end select

         ! Test different moments
         do moment = 1, 3
            break = breakterm(bfnc=bconst, dfnc=dfuni, moment=moment, grid=gx, &
                              update_b=.false.)
            np = ZERO
            np(nc) = ONE
            y = ZERO
            call break%eval(np, y, udot_birth=birth, udot_death=death)

            moment_birth_0 = sum(birth)
            moment_death_0 = sum(death)
            moment_birth_m = sum(birth*gx%center**moment)
            moment_death_m = sum(death*gx%center**moment)

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

end module test_break1
