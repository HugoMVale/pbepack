program demo
   use ftlListIntModule
   implicit none
   type(ftlListInt) :: intlist

   call intlist%New

   block
      integer :: it
      do it = 1, 10
         call intlist%PushBack(it*it)
      end do
      print *, intlist%Size()
   end block

   block
      integer, allocatable :: arr(:)
      integer :: i
      allocate (arr(intlist%Size()))
      call intlist%CopyToArray(arr)
      print '(*(g0:, ",", 1x))', arr
   end block

   call intlist%Delete
end program demo
