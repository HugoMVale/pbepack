module pbepack_lib
!! Auxiliary functions.
   use pbepack_kinds
   use stdlib_optval, only: optval
   use stdlib_strings, only: to_string
   implicit none
   public

contains

   elemental real(rk) function delta_kronecker(i, j) result(res)
   !! Delta kronecker \( \delta_{i,j} \).
      integer, intent(in) :: i, j
         !! integer

      res = ZERO
      if (i == j) res = ONE

   end function delta_kronecker

   elemental real(rk) function boxcar(x, a, b, height) result(res)
   !! Boxcar function.
      real(rk), intent(in) :: x
         !! argument
      real(rk), intent(in) :: a, b
         !! interval limit \( a < b \)
      real(rk), intent(in) :: height
         !! pulse height

      res = ZERO
      if (x > a .and. x < b) res = height

   end function boxcar

   elemental real(rk) function expo1d(x, x0, n0) result(res)
   !! 1D Exponential distribution.
      real(rk), intent(in) :: x
         !! random variable
      real(rk), intent(in), optional :: x0
         !! mean value
      real(rk), intent(in), optional :: n0
         !! initial number of particles
      real(rk) :: x0_, n0_
      x0_ = optval(x0, ONE)
      n0_ = optval(n0, ONE)
      res = (n0_/x0_)*exp(-x/x0_)
   end function expo1d

   subroutine writearray(filename, x, fmt)
   !! Write content of real array (rank=1,2) to file.
      character(*), intent(in) :: filename
         !! file name
      real(rk), intent(in) :: x(..)
         !! array to be written
      character(*), intent(in), optional :: fmt
         !! format specifier (default="es16.5e3")

      character(:), allocatable :: fmt_
      integer :: i, funit

      fmt_ = optval(fmt, "es16.5e3")

      open (newunit=funit, file=filename, status="replace", action="write", position="rewind")

      select rank (x)
      rank (1)
         if (size(x, 1) > 0) then
            do i = 1, size(x, 1)
               write (funit, '('//fmt_//')') x(i)
            end do
         end if
      rank (2)
         if (size(x, 1) > 0 .and. size(x, 2) > 0) then
            do i = 1, size(x, 1)
               write (funit, '('//to_string(size(x, 2))//'('//fmt_//'))') x(i, :)
            end do
         end if
      end select

      close (funit)

   end subroutine writearray

   pure function buildfilename(folder, basename, suffix, extension) result(res)
   !! Build file name from its parts.
      character(*), intent(in), optional :: folder
         !! folder (default=".\")
      character(*), intent(in), optional :: basename
         !! basename (default="temp")
      character(*), intent(in), optional :: suffix
         !! suffix (default="")
      character(*), intent(in), optional :: extension
         !! extension (default="txt")
      character(:), allocatable :: res
      res = optval(folder, ".\")//optval(basename, "temp")
      if (present(suffix)) then
         res = res//"_"//suffix
      end if
      res = res//"."//optval(extension, "txt")
   end function

end module pbepack_lib
