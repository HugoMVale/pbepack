module combtypes
!!   This module implements two derived data types to store information about particle
!! combinations, as required for evaluating aggregation birth terms.
   use real_kinds, only: rk
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
      procedure, pass(self) :: new
      procedure, pass(self) :: size => list_size
      procedure, pass(self) :: append
      procedure, pass(self) :: toarray
      procedure, pass(self) :: clear
   end type comblist

   type :: combarray
      ! Array of particle combinations
      integer, allocatable :: ia(:)
         !! index of particle 'a'
      integer, allocatable :: ib(:)
         !! index of particle 'b'
      real(rk), allocatable :: weight(:)
         !! weight of 'a+b' assigned to pivot
   contains
      procedure, pass(self) :: alloc
      procedure, pass(self) :: size => array_size
   end type combarray

contains

   impure subroutine new(self)
   !! Constructor comblist.
      class(comblist), intent(inout) :: self
         !! Object
      call self%ia%New
      call self%ib%New
      call self%weight%New
   end subroutine new

   pure integer function list_size(self) result(res)
   !! Size(comblist).
      class(comblist), intent(in) :: self
         !! Object
      res = self%ia%Size()
   end function list_size

   impure subroutine append(self, ia, ib, weight)
   !! Append values to comblist.
      class(comblist), intent(inout) :: self
         !! object
      integer :: ia
         !! index particle a
      integer :: ib
         !! index particle b
      real(rk) :: weight
         !! weight
      call self%ia%PushBack(ia)
      call self%ib%PushBack(ib)
      call self%weight%PushBack(weight)
   end subroutine append

   impure subroutine clear(self)
   !! Clear comblist.
      class(comblist), intent(inout) :: self
         !! object
      call self%ia%Clear
      call self%ib%Clear
      call self%weight%Clear
   end subroutine clear

   impure subroutine toarray(self, array)
   !! Copy contents of comblist to array of same type.
      class(comblist), intent(inout) :: self
         !! object
      type(combarray), intent(inout) :: array
         !! target array
      call self%ia%CopyToArray(array%ia)
      call self%ib%CopyToArray(array%ib)
      call self%weight%CopyToArray(array%weight)
   end subroutine toarray

   pure subroutine alloc(self, n)
   !! Allocate combarray.
      class(combarray), intent(inout) :: self
         !! Object
      integer, intent(in) :: n
         !! Array size
      allocate (self%ia(n))
      allocate (self%ib(n))
      allocate (self%weight(n))
   end subroutine alloc

   pure integer function array_size(self) result(res)
   !! Size(combarray).
      class(combarray), intent(in) :: self
         !! Object
      res = size(self%ia)
   end function array_size

end module
