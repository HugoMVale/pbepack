module test_agg1
!! Test for module 'agg1d' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds
   use pbepack_pbe1, only: pbe
   use hrweno_grids, only: grid1
   use utils_tests
   use stdlib_strings, only: to_string
   implicit none
   private

   public :: collect_tests_agg1

   logical, parameter :: verbose = .true.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_agg1(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("moment conservation", test_moment_conservation) &
                  ]

   end subroutine

   subroutine test_moment_conservation(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 200
      type(grid1) :: gx
      type(pbe) :: eq
      real(rk), dimension(nc) :: u, birth, death
      real(rk) :: y(0:0), moment_birth_0, moment_birth_m, moment_death_0, moment_death_m
      integer :: moment, scale

      ! Test linear and log grids
      do scale = 1, 2
         select case (scale)
         case (1)
            call gx%linear(1._rk, 1e3_rk, nc)
         case (2)
            call gx%log(1._rk, 1e3_rk, nc)
         end select

         ! Test different moments
         do moment = 1, 3
            eq = pbe(grid=gx, afnc=aprod, moment=moment, update_a=.false., &
                     name="test_moment_conservation")
            u = ZERO; u(1:nc/2 - 1) = ONE
            y = ZERO
            call eq%agg%eval(u, y, udot_birth=birth, udot_death=death)

            moment_birth_0 = sum(birth*gx%width)
            moment_death_0 = sum(death*gx%width)
            moment_birth_m = sum(birth*gx%width*gx%center**moment)
            moment_death_m = sum(death*gx%width*gx%center**moment)

            call check(error, moment_birth_0, moment_death_0/2, rel=.true., thr=1e-14_rk)
            call check(error, moment_birth_m, moment_death_m, rel=.true., thr=1e-14_rk)

            if (allocated(error) .or. verbose) then
               print *
               write (stderr, '(a18,a24)') "scale type       =", gx%scale
               write (stderr, '(a18,i24)') "preserved moment =", moment
               write (stderr, '(a18,(es24.14e3))') "birth/death(0)   =", moment_birth_0/moment_death_0
               write (stderr, '(a18,(es24.14e3))') "birth/death("//to_string(moment)//")   =", moment_birth_m/moment_death_m
            end if
            if (allocated(error)) return

         end do
      end do

   end subroutine test_moment_conservation

   subroutine test_case1(error)
      type(error_type), allocatable, intent(out) :: error

      !    integer, parameter :: nc = 100
      !    type(grid1) :: gx
      !    type(combarray) :: aggcomb(nc)
      !    real(rk), dimension(nc) :: n, birth, death
      !    real(rk) :: t, y(0:0), t0, tend
      !    integer :: i, m, scl

      ! call cpu_time(t0)
      !    ! ! System parameters
      !    ! n0 = 1._rk
      !    ! x0 = 1._rk
      !    ! a0 = 1._rk

      !    ! ! Grid
      !    ! call gx%new(0._rk, 10*x0, nc, scl=1)

      !    ! ! Initial condition
      !    ! n = cell_average(expo1d, gx)
      !    ! time = 0

      !    ! Call solver

      ! call cpu_time(tend)
      ! print *
      ! print '("Time = ",f8.5," seconds.")', (tend - t0)
      ! print *

   end subroutine test_case1

end module test_agg1
