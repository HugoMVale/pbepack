module test_break1
!! Test for module 'break1d' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use real_kinds, only: rk
   use break1, only: breakterm
   use grids, only: grid1
   use utils_tests, only: bconst, bsquare, duniform
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

      integer, parameter :: nc = 200
      type(grid1) :: gx
      type(breakterm) :: break
      real(rk), dimension(nc) :: np, source, sink
      real(rk) :: y(0:0), t0, tend, sum_source, sum_sink
      integer :: m, scl

      call cpu_time(t0)

      ! Test linear and log grids
      do scl = 1, 2
         select case (scl)
         case (1)
            call gx%linear(1._rk, 1e3_rk, nc)
         case (2)
            call gx%log(1._rk, 1e3_rk, nc)
         end select

         ! Test different moments
         do m = 1, 2
            break = breakterm(bf=bsquare, df=duniform, moment=m, grid=gx)
            np = 0
            np(nc - 2) = 2
            y = 0._rk
            call break%eval(np, y, source=source, sink=sink)
            sum_source = sum(source*gx%center**m)
            sum_sink = sum(sink*gx%center**m)
            call check(error, sum_source, sum_sink, rel=.true., thr=1e-2_rk)

            if (allocated(error) .or. verbose) then
               write (stderr, '(a13,i1,a1,a2,es26.16e3)') "sum_source  (", m, ")", "=", sum_source
               write (stderr, '(a13,i1,a1,a2,es26.16e3)') "sum_sink    (", m, ")", "=", sum_sink
               write (stderr, '(a17,es26.16e3)') "source/sink (0) =", sum(source)/sum(sink)
            end if
            if (allocated(error)) return
         end do

      end do

      call cpu_time(tend)
      print '("Time = ",f8.5," seconds.")', (tend - t0)

   end subroutine test_mass_conservation

end module test_break1
