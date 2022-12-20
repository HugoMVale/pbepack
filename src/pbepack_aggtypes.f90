module pbepack_aggtypes
!!   Derived data types to store information about particle combinations, as required for
!! evaluating aggregation birth terms.
   use pbepack_kinds
   use ftlListIntModule
   use ftlListDoubleModule
   implicit none
   private

   public :: comblist, combarray

   type :: comblist
   !! List of particle combinations.
      type(ftlListInt) :: ia
         !! index of particle 'a'
      type(ftlListInt) :: ib
         !! index of particle 'b'
      type(ftlListDouble) :: weight
         !! weight of 'a+b' assigned to pivot
   contains
      procedure, pass(self) :: new => list_new
      procedure, pass(self) :: size => list_size
      procedure, pass(self) :: append => list_append
      procedure, pass(self) :: toarray => list_toarray
      procedure, pass(self) :: clear => list_clear
   end type comblist

   type :: combarray
   !! Array of particle combinations
      integer, allocatable :: ia(:)
         !! index of particle 'a'
      integer, allocatable :: ib(:)
         !! index of particle 'b'
      real(rk), allocatable :: weight(:)
         !! weight of 'a+b' assigned to pivot
   contains
      procedure, pass(self) :: alloc => array_alloc
      procedure, pass(self) :: size => array_size
   end type combarray

contains

   impure subroutine list_new(self)
   !! Constructor `comblist`.
      class(comblist), intent(inout) :: self
         !! object

      call self%ia%New
      call self%ib%New
      call self%weight%New

   end subroutine list_new

   pure integer function list_size(self) result(res)
   !! Size of `comblist`.
      class(comblist), intent(in) :: self
         !! object

      res = self%ia%Size()

   end function list_size

   impure subroutine list_append(self, ia, ib, weight)
   !! Append values to `comblist` object.
      class(comblist), intent(inout) :: self
         !! object
      integer :: ia, ib
         !! particle index
      real(rk) :: weight
         !! weight

      call self%ia%PushBack(ia)
      call self%ib%PushBack(ib)
      call self%weight%PushBack(weight)

   end subroutine list_append

   impure subroutine list_clear(self)
   !! Clear `comblist` object.
      class(comblist), intent(inout) :: self
         !! object

      call self%ia%Clear
      call self%ib%Clear
      call self%weight%Clear

   end subroutine list_clear

   impure subroutine list_toarray(self, array)
   !! Copy contents of `comblist` object to array of same type.
      class(comblist), intent(inout) :: self
         !! object
      type(combarray), intent(inout) :: array
         !! target array

      call self%ia%CopyToArray(array%ia)
      call self%ib%CopyToArray(array%ib)
      call self%weight%CopyToArray(array%weight)

   end subroutine list_toarray

   pure subroutine array_alloc(self, n)
   !! Allocate `combarray`.
      class(combarray), intent(inout) :: self
         !! object
      integer, intent(in) :: n
         !! array size

      allocate (self%ia(n))
      allocate (self%ib(n))
      allocate (self%weight(n))

   end subroutine array_alloc

   pure integer function array_size(self) result(res)
   !! Size of `combarray`.
      class(combarray), intent(in) :: self
         !! object

      res = size(self%ia)

   end function array_size

end module pbepack_aggtypes
