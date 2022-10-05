module test_agg1d
!! Test for module 'agg1d' using test-drive.
   use iso_fortran_env, only: real64, stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use agg1d, only: aggterm
   use grid, only: grid1
   implicit none
   private

   public :: collect_tests_agg1d

   integer, parameter :: rk = real64
   logical, parameter :: verbose = .true.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_agg1d(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mass conservation", test_mass_conservation), &
                  new_unittest("analytical solution, case 1", test_case1) &
                  !new_unittest("calc_c", test_calc_c), &
                  !new_unittest("wenok with non-uniform grid", test_wenok_nonuniform) &
                  ]

   end subroutine

   subroutine test_mass_conservation(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 1000
      type(grid1) :: gx
      type(aggterm) :: agg
      real(rk), dimension(nc) :: np, birth, death
      real(rk) :: t, y(0:0), t0, tend, sum_birth, sum_death
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
         do m = 1, 3
            call agg%init(gx, m)
            np = 0
            np(1:nc/2 - 1) = 1
            t = 0._rk
            y = 0._rk
            call agg%eval(np, aconst, t, y, birth, death)
            sum_birth = sum(birth*gx%center**m)
            sum_death = sum(death*gx%center**m)
            call check(error, sum_birth, sum_death, rel=.true., thr=1e-14_rk)

            if (allocated(error) .or. verbose) then
               write (stderr, '(a12,es26.16e3)'), "sum_birth= ", sum_birth
               write (stderr, '(a12,es26.16e3)'), "sum_death= ", sum_death
            end if
            if (allocated(error)) return
         end do
      end do

      call cpu_time(tend)
      print '("Time = ",f8.5," seconds.")', (tend - t0)

   end subroutine test_mass_conservation

   subroutine test_case1(error)
      type(error_type), allocatable, intent(out) :: error

      !    integer, parameter :: nc = 100
      !    type(grid1) :: gx
      !    type(combarray) :: aggcomb(nc)
      !    real(rk), dimension(nc) :: n, birth, death
      !    real(rk) :: t, y(0:0), t0, tend
      !    integer :: i, m, scl

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

   end subroutine test_case1

   pure real(rk) function aconst(xa, xb, t, y) result(res)
   !! Constant aggregation kernel for 1D system
      real(rk), intent(in) :: xa
         !! internal coordinate of particle a
      real(rk), intent(in) :: xb
         !! internal coordinate of particle b
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), parameter :: a0 = 1._rk
      res = a0
   end function

   pure real(rk) function asum(xa, xb, t, y) result(res)
   !! Sum aggregation kernel for 1D system
      real(rk), intent(in) :: xa
         !! internal coordinate of particle a
      real(rk), intent(in) :: xb
         !! internal coordinate of particle b
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), parameter :: a1 = 1._rk
      res = a1*(xa + xb)
   end function

   pure real(rk) function aprod(xa, xb, t, y) result(res)
   !! Product aggregation kernel for 1D system
      real(rk), intent(in) :: xa
         !! internal coordinate of particle a
      real(rk), intent(in) :: xb
         !! internal coordinate of particle b
      real(rk), intent(in) :: t
         !! time
      real(rk), intent(in) :: y(:)
         !! environment vector
      real(rk), parameter :: a1 = 1._rk
      res = a1*xa*xb
   end function

   elemental real(rk) function expo1d(x, x0, n0)
   !! 1D Exponential distribution.
      real(rk), intent(in) :: x
         !! random variable
      real(rk), intent(in) :: x0
         !! mean value
      real(rk), intent(in) :: n0
         !! initial number of particles

      expo1d = (n0/x0)*exp(-x/x0)

   end function expo1d

   elemental real(rk) function solution_case1(x, x0, n0, a0, t) result(res)
   !! Anylytical solution for IC='expo1d' and constant aggregation frequency.
      real(rk), intent(in) :: x
         !! random variable
      real(rk), intent(in) :: x0
         !! mean value
      real(rk), intent(in) :: n0
         !! initial number of particles
      real(rk), intent(in) :: a0
         !! aggregation frequency
      real(rk), intent(in) :: t
         !! time

      real(rk) :: tau

      tau = n0*a0*t
      res = 4*n0/(x0*(tau + 2._rk)**2)*exp(-2*(x/x0)/(tau + 2._rk))

   end function solution_case1

end module test_agg1d
