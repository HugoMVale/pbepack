module test_break1
!! Test for module 'break1d' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_strings, only: to_string
   use real_kinds, only: rk
   use break1, only: breakterm
   use grids, only: grid1
   use utils_tests, only: bconst, duniform
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
                  new_unittest("mass conservation", test_mass_conservation) &
                  !new_unittest("analytical solution, case 1", test_case1) &
                  !new_unittest("wenok with non-uniform grid", test_wenok_nonuniform) &
                  ]

   end subroutine

   subroutine test_mass_conservation(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 500
      type(grid1) :: gx
      type(breakterm) :: break
      real(rk), dimension(nc) :: np, source, sink
      real(rk) :: y(0:0), moment_source_0, moment_source_m, moment_sink_0, moment_sink_m
      real(rk) :: t0, tend
      integer :: moment, scl

      call cpu_time(t0)

      ! Test linear and log grids
      do scl = 1, 2
         select case (scl)
         case (1)
            call gx%linear(0._rk, 1e2_rk, nc)
         case (2)
            call gx%geometric(0._rk, 2e3_rk, 1.01_rk, nc)
         end select

         ! Test different moments
         do moment = 1, 3
            break = breakterm(bf=bconst, df=dfuni, moment=moment, grid=gx)
            np = 0
            np(nc - 2) = 2
            y = 0._rk
            call break%eval(np, y, source=source, sink=sink)

            moment_source_0 = sum(source)
            moment_sink_0 = sum(sink)
            moment_source_m = sum(source*gx%center**moment)
            moment_sink_m = sum(sink*gx%center**moment)

            call check(error, moment_source_m, moment_sink_m, &
                       rel=.true., thr=1e-5_rk)
            call check(error, moment_source_0, moment_sink_0*2, &
                       rel=.true., thr=1e-2_rk)

            if (allocated(error) .or. verbose) then
               print *
               write (stderr, '(a18,3x,i1)') "scale type       =", gx%scl
               write (stderr, '(a18,3x,i1)') "preserved moment =", moment
               write (stderr, '(a18,(es24.14e3))') "source/sink(0)   =", &
                  moment_source_0/moment_sink_0
               write (stderr, '(a18,(es24.14e3))') "source/sink("//to_string(moment)//")   =", &
                  moment_source_m/moment_sink_m
            end if
            if (allocated(error)) return
         end do

      end do

      call cpu_time(tend)
      print *
      print '("Time = ",f8.5," seconds.")', (tend - t0)
      print *

   contains

      pure real(rk) function dfuni(xd, xo, y) result(res)
         real(rk), intent(in) :: xd, xo, y(:)
         res = duniform(xd, xo, y, moment)
      end function

   end subroutine test_mass_conservation

end module test_break1
