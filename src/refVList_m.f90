! Copyright (C) 2012 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later


! vlist template (vector list)
! template list based on array of pointers 
!
! insert/append insert *copies* of the data into the vlist
! therefore a copy constructor (and overloaded assignment) is required for the data type
! if deallocation of the data type is required set  call dp%destroy() appropriately to calling the destructor for pointer variable dp
!
! m4 -D_MODNAME_=IntVListModule -D_PRE_=int_ -D_TYPE_=integer -D_DA_= -D_USEMOD_= -D_IMP_= vlist_tmpl_m.m4 > int_vlist_m.f90
! m4 -D_MODNAME_=RefListModule -D_PRE_=refl_ -D_TYPE_=type(reference) -D_DA_="call dp%destroy()" -D_USEMOD_="use referenceADT" -D_IMP_="import reference"  vlist_tmpl_m.m4 > ref_vlist_m.f90
!
! careful: sort (with quick sort) requires in external "isSmaller" logical function 
! for increasing order (or "isGreater" for decreasing order)
! This function *must* have <= or >= semantics!!
!
module refVList_m

   use refADT_m
   implicit none
   
   type refl_data_pointer
      private
      type(reference), pointer :: p_data => null()
   end type 

   type refl_vlist
      private
      type(refl_data_pointer), pointer :: vector(:) => null()
      integer          :: cursize=0
      integer          :: maxSize=0
   contains
      procedure :: create => refl_create
      procedure :: destroy => refl_destroy
      procedure :: resize => refl_resize
      procedure :: isEmpty => refl_isEmpty
      procedure :: size => refl_size
      procedure :: elem => refl_elem
      procedure :: append => refl_append
      procedure :: insert => refl_insert
      procedure :: del => refl_del
      procedure :: delFirst => refl_delFirst
      procedure :: delLast => refl_delLast
      procedure :: sort => refl_sort
   end type

   interface assignment(=)
      module procedure refl_copy
   end interface

contains

   subroutine refl_create(vl,maxsize)
      class(refl_vlist)    :: vl
      integer, optional, intent(in) :: maxsize
      integer ms
      if (present(maxsize)) then 
         ms = maxsize
      else
         ms = 10
      end if
      allocate(vl%vector(ms))
      vl%maxsize = ms
      vl%cursize = 0
   end subroutine

   subroutine refl_copy(s,rhs)
      class(refl_vlist), intent(inout)  :: s
      type(refl_vlist), intent(in)  :: rhs
      type(reference), pointer :: rp
      integer i
      !DBG!call assert(associated(rhs%vector),"refl_copy: rhs not initialized")
      if (.not.associated(s%vector)) then 
         call refl_create(s,rhs%maxsize)
      else
         call refl_del(s)
      end if
      do i=1,rhs%cursize
         rp => rhs%vector(i)%p_data
         call refl_append(s,rp)
      end do
   end subroutine refl_copy

   subroutine refl_destroy(vl)
      class(refl_vlist)    :: vl
      call refl_del(vl)
      deallocate(vl%vector)
      vl%maxsize = 0
      vl%cursize = 0
   end subroutine

   subroutine refl_resize(vl,maxsize)
      class(refl_vlist)    :: vl
      integer, optional, intent(in) :: maxsize
      type(refl_data_pointer), pointer :: vp(:)
      integer ms,i
      real, parameter :: grow_factor=1.3
      if (present(maxsize)) then 
         ms = max(maxsize,vl%cursize)
      else
         ms = max(int(vl%maxSize*grow_factor),vl%maxSize+1)
      end if
      allocate(vp(ms))
      do i=1,min(ms,vl%cursize)
         vp(i)%p_data => vl%vector(i)%p_data
      end do
      deallocate(vl%vector)
      vl%vector => vp
      vl%maxsize = ms
   end subroutine refl_resize

   logical function refl_isEmpty(vl)
      class(refl_vlist) :: vl
      refl_isEmpty = vl%cursize == 0
   end function refl_isEmpty

   integer function refl_size(vl)
      class(refl_vlist) :: vl
      refl_size = vl%cursize
   end function refl_size

   function refl_elem(vl,i)
      class(refl_vlist) :: vl
      integer, intent(in) :: i
      type(reference), pointer :: refl_elem
      refl_elem => vl%vector(i)%p_data
   end function

   subroutine refl_append(vl,d)
      class(refl_vlist) :: vl
      type(reference), intent(in) :: d
      vl%cursize = vl%cursize + 1
      if (vl%cursize > vl%maxsize) call refl_resize(vl)
      allocate(vl%vector(vl%cursize)%p_data)
      vl%vector(vl%cursize)%p_data = d
   end subroutine refl_append

   subroutine refl_insert(vl,p,d)
      class(refl_vlist) :: vl
      integer, intent(in) :: p   ! position where to insert
      type(reference), intent(in) :: d   ! data
      integer i
      vl%cursize = vl%cursize + 1
      if (vl%cursize > vl%maxsize) call refl_resize(vl)
      do i=vl%cursize,p+1,-1
         vl%vector(i)%p_data => vl%vector(i-1)%p_data
      end do
      nullify(vl%vector(p)%p_data)
      allocate(vl%vector(p)%p_data)
      vl%vector(p)%p_data = d
   end subroutine refl_insert

   subroutine refl_del(vl,p)
      class(refl_vlist) :: vl
      integer, intent(in), optional :: p  ! position to delete
      type(reference),pointer :: dp
      integer i
      if (present(p)) then
         dp => vl%vector(p)%p_data
         call dp%destroy()
         deallocate(dp)
         do i=p,vl%cursize-1
            vl%vector(i)%p_data => vl%vector(i+1)%p_data
         end do
         vl%cursize = vl%cursize - 1
      else
         do i=1,vl%cursize
            dp => vl%vector(i)%p_data
            call dp%destroy()
            deallocate(dp)
         end do
         vl%cursize = 0
      end if
   end subroutine refl_del

   subroutine refl_delFirst(vl)
      class(refl_vlist) :: vl
      type(reference),pointer :: dp
      integer i
      dp => vl%vector(1)%p_data
      call dp%destroy()
      deallocate(dp)
      do i=1,vl%cursize-1
         vl%vector(i)%p_data => vl%vector(i+1)%p_data
      end do
      vl%cursize = vl%cursize - 1
   end subroutine refl_delFirst

   subroutine refl_delLast(vl)
      class(refl_vlist) :: vl
      type(reference),pointer :: dp
      dp => vl%vector(vl%cursize)%p_data
      call dp%destroy()
      deallocate(dp)
      vl%cursize = vl%cursize - 1
   end subroutine refl_delLast

   subroutine refl_sort(vl,cmp)
      class(refl_vlist) :: vl
      interface
         logical function cmp(i,j)
         import :: reference
         type(reference),intent(in) :: i,j
         end function cmp
      end interface
      call refl_quicksort(vl%vector(1:vl%cursize),cmp)
   end subroutine refl_sort

   recursive subroutine refl_quicksort(list,cmp)
      type(refl_data_pointer) :: list(:)
      interface
         logical function cmp(i,j)
         import :: reference
         type(reference),intent(in) :: i,j
         end function cmp
      end interface
      integer :: i,j,n
      type(reference), pointer :: chosen,temp  ! TYPE
      integer, parameter :: maxSimpleSortSize=3
      
      n = size(list)
      if (n <= maxSimpleSortSize) then
         call refl_simpleSort(list,cmp)
      else
         chosen => list(n/2)%p_data
         i = 0
         j = n+1
         do
            do
               i = i+1
               if (cmp(chosen,list(i)%p_data)) exit 
            enddo 
            do 
               j = j-1
               if (cmp(list(j)%p_data,chosen)) exit
            enddo
            if (i<j) then
               ! swap
               temp => list(i)%p_data
               list(i)%p_data => list(j)%p_data
               list(j)%p_data => temp
            else if (i==j) then
               i = i+1
               exit
            else
               exit
            end if
         end do
         if (j>1) call refl_quicksort(list(:j),cmp)
         if (i<n) call refl_quicksort(list(i:),cmp)
      end if
   end subroutine refl_quicksort
   
   subroutine refl_simpleSort(list,cmp)
      type(refl_data_pointer), intent(inout) :: list(:)
      interface
         logical function cmp(i,j)
         import :: reference
         type(reference),intent(in) :: i,j
         end function cmp
      end interface
      integer :: i,j
      type(reference),pointer  :: temp  
      do i=1,size(list)-1
         do j=i+1,size(list)
            if (cmp(list(j)%p_data,list(i)%p_data)) then
               ! swap
               temp => list(i)%p_data
               list(i)%p_data => list(j)%p_data
               list(j)%p_data => temp
            end if
         end do
      end do
   end subroutine refl_simpleSort
         
end module refVList_m

