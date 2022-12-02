module pbepack_algebra
!! This module implements special array types.
   use pbepack_kinds
   implicit none
   private

   public:: spmatrix

   type :: spmatrix
   !! Symmetric array class
      real(rk), allocatable :: ap(:)
         !! vector with array values in packed storage format
      integer :: n
         !! number of rows or columns
      character(1) :: uplo = "u"
         !! flag to specify whether the upper(u) or lower(l) triangle is supplied
   contains
      procedure, pass(self) :: multvec => spmatrix_multvec
      procedure, pass(self) :: get => spmatrix_get
      procedure, pass(self) :: set => spmatrix_set
   end type spmatrix

   interface spmatrix
      module procedure :: spmatrix_init
   end interface spmatrix

contains

   pure type(spmatrix) function spmatrix_init(n) result(res)
      !! Initialize 'spmatrix' object.
      integer, intent(in) :: n
         !! number of rows or columns

      if (n > 0) then
         res%n = n
      else
         error stop "Invalid 'n'. Valid range: n > 0"
      end if

      allocate (res%ap(n*(n + 1)/2))

   end function spmatrix_init

   elemental real(rk) function spmatrix_get(self, i, j) result(res)
   !! Get value of \( a_{i,j} \).
      class(spmatrix), intent(in) :: self
         !! object
      integer, intent(in) :: i, j
         !! index

      res = self%ap(i + j*(j - 1)/2)

   end function spmatrix_get

   pure subroutine spmatrix_set(self, i, j, x)
   !! Set value of \( a_{i,j} = x \).
      class(spmatrix), intent(inout) :: self
         !! object
      integer, intent(in) :: i, j
         !! index
      real(rk), intent(in) :: x
         !! value

      self%ap(i + j*(j - 1)/2) = x

   end subroutine spmatrix_set

   pure function spmatrix_multvec(self, x) result(y)
      !! Performs the matrix-vector operation:
      !! ```
      !! y:= A*x
      !! ```
      !! where x is an n element vector and A is a n*n symmetric matrix.
      !! Source: https://netlib.org/lapack/lug/node123.html
      class(spmatrix), intent(in) :: self
            !! symmetric array(n,n)
      real(rk), intent(in) :: x(:)
            !! vector(n)
      real(rk) :: y(size(x))

      integer :: i, j, k, kk
      real(rk) :: temp1, temp2

      if (self%n /= size(x)) error stop "Dimension mismatch."

      ! Implemented for upper triangle only
      associate (ap => self%ap)
         y = ZERO
         kk = 1
         do j = 1, self%n
            temp1 = x(j)
            temp2 = ZERO
            k = kk
            do i = 1, j - 1
               y(i) = y(i) + temp1*ap(k)
               temp2 = temp2 + ap(k)*x(i)
               k = k + 1
            end do
            y(j) = y(j) + temp1*ap(kk + j - 1) + temp2
            kk = kk + j
         end do
      end associate

   end function spmatrix_multvec

end module pbepack_algebra
