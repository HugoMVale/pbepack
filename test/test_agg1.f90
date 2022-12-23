module test_agg1
!! Test for module 'agg1d' using test-drive.
   use iso_fortran_env, only: stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pbepack_kinds
   use pbepack_pbe1, only: pbe, pbesol
   use pbepack_quadratures, only: evalmoment
   use hrweno_grids, only: grid1
   use utils_tests
   use stdlib_strings, only: to_string
   use stdlib_math, only: linspace
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
                  new_unittest("moment conservation", test_moment_conservation), &
                  new_unittest("case a(x,x')=1", test_aconst), &
                  new_unittest("case a(x,x')=x+x'", test_asum) &
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
            call gx%linear(1._rk, 1e3_rk, ncells)
         case (2)
            call gx%log(1._rk, 1e3_rk, ncells)
         end select

         ! Test different moments
         do moment = 1, 3

            mypbe = pbe(grid=gx, a=aprod, moment=moment, update_a=.false., &
                        name="test_moment_conservation")
            u = ZERO; u(1:ncells/2 - 1) = ONE
            call mypbe%agg%eval(u, y=VOIDREAL, udot_birth=udot_birth, udot_death=udot_death)

            moment_birth_0 = evalmoment(udot_birth, gx, 0)
            moment_death_0 = evalmoment(udot_death, gx, 0)
            moment_birth_m = evalmoment(udot_birth, gx, moment)
            moment_death_m = evalmoment(udot_death, gx, moment)

            call check(error, moment_birth_0, moment_death_0/2, rel=.true., thr=10*EPS)
            call check(error, moment_birth_m, moment_death_m, rel=.true., thr=10*EPS)

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

   subroutine test_aconst(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 100
      type(grid1) :: gx
      type(pbe) :: equation
      type(pbesol) :: solution
      real(rk) :: u(nc), n0, x0
      real(rk), dimension(0:1) :: moment, momoment_ref
      real(rk), allocatable :: times(:)
      integer :: i

      ! Init linear grid
      call gx%log(1e-2_rk, 1e2_rk, nc)

      ! Init pbe object
      equation = pbe(gx, a=aconst, name="test_aconst")

      ! Integrate
      times = linspace(ZERO, TWO, 2)
      solution = equation%integrate(times, u0)

      ! Compute solution moments
      do i = 0, 1
         moment(i) = evalmoment(solution%u(:, size(times)), gx, i)
      end do

      ! Compute reference moments
      n0 = evalmoment(solution%u(:, 1), gx, 0)
      x0 = evalmoment(solution%u(:, 1), gx, 1)/n0
      momoment_ref = aconst_moments(times(size(times)), x0=x0, n0=n0)

      ! Check moments
      call check(error, moment(0), momoment_ref(0), rel=.true., thr=1e-6_rk)
      call check(error, moment(1), momoment_ref(1), rel=.true., thr=100*EPS)

      if (allocated(error) .or. verbose) then
         do i = 0, 1
            write (stderr, '(a18,(es24.14e3))') &
               "num./analyt.("//to_string(i)//") =", moment(i)/momoment_ref(i)
         end do
         call solution%write("test_aconst")
      end if

   end subroutine test_aconst

   subroutine test_asum(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 200
      type(grid1) :: gx
      type(pbe) :: equation
      type(pbesol) :: solution
      real(rk) :: u(nc), n0, x0
      real(rk), dimension(0:1) :: moment, moment_ref
      real(rk), allocatable :: times(:)
      integer :: i

      ! Init linear grid
      call gx%log(1e-2_rk, 1e4_rk, nc)

      ! Init pbe object
      equation = pbe(gx, a=asum, name="test_asum")

      ! Integrate
      times = linspace(ZERO, TWO, 2)
      solution = equation%integrate(times, u0)

      ! Compute solution moments
      do i = 0, 1
         moment(i) = evalmoment(solution%u(:, size(times)), gx, i)
      end do

      ! Compute reference moments
      n0 = evalmoment(solution%u(:, 1), gx, 0)
      x0 = evalmoment(solution%u(:, 1), gx, 1)/n0
      moment_ref = asum_moments(times(size(times)), x0=x0, n0=n0)

      ! Check moments
      call check(error, moment(0), moment_ref(0), rel=.true., thr=1e-6_rk)
      call check(error, moment(1), moment_ref(1), rel=.true., thr=1e2*EPS)

      if (allocated(error) .or. verbose) then
         do i = 0, 1
            write (stderr, '(a18,(es24.14e3))') &
               "num./analyt.("//to_string(i)//") =", moment(i)/moment_ref(i)
         end do
         call solution%write("test_asum")
      end if

   end subroutine test_asum

   pure real(rk) function u0(x) result(res)
      real(rk), intent(in) :: x
      res = expo1d(x)
   end function

end module test_agg1
