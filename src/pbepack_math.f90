module pbepack_math
!! This module implements some auxiliary math tools.
   use pbepack_kinds
   implicit none
   private

   public:: delta_kronecker, spmatrix

   type :: spmatrix
   !! Symmetric array class
      real(rk), allocatable :: ap(:)
         !! vector with array values in packed storage format
      integer :: n
         !! number of rows/columns
      character(1) :: uplo = "u"
         !! flag to specify whether the upper(u) or lower(l) triangle is supplied
   contains
      procedure, pass(self) :: multvec => spmatrix_multvec
   end type spmatrix

   interface spmatrix
      module procedure :: spmatrix_init
   end interface spmatrix

contains

   pure real(rk) function delta_kronecker(i, j) result(res)
   !! Delta kronecker \( \delta_{i,j} \).
      integer, intent(in) :: i, j
         !! integers

      if (i == j) then
         res = ONE
      else
         res = ZERO
      end if

   end function delta_kronecker

   pure type(spmatrix) function spmatrix_init(n) result(self)
   !! Initialize 'spmatrix' object.
      integer, intent(in) :: n
         !! number of rows/columns

      if (n > 0) then
         self%n = n
      else
         error stop "Invalid 'n'. Valid range: n > 0"
      end if

      allocate (self%ap(n*(n + 1)/2))

   end function spmatrix_init

   pure function spmatrix_multvec(self, x) result(y)
   !! Performs the matrix-vector operation:
   !! ```
   !! y:= A*x
   !! ```
   !! where x is an n element vector and A is a n*n symmetric matrix supplied in packed form.
      class(spmatrix), intent(in) :: self
        !! symmetric array supplied in packed storage format (upper triangle)
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

end module pbepack_math
