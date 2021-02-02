! Copyright (C) 2012-2013 Arne Luechow
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
module reflistVList_m

   use refVList_m
   implicit none
   
   type refll_data_pointer
      private
      type(refl_vlist), pointer :: p_data => null()
   end type 

   type refll_vlist
      private
      type(refll_data_pointer), pointer :: vector(:) => null()
      integer          :: cursize=0
      integer          :: maxSize=0
   contains
      procedure :: create => refll_create
      procedure :: destroy => refll_destroy
      procedure :: resize => refll_resize
      procedure :: isEmpty => refll_isEmpty
      procedure :: size => refll_size
      procedure :: elem => refll_elem
      procedure :: append => refll_append
      procedure :: insert => refll_insert
      procedure :: del => refll_del
      procedure :: delFirst => refll_delFirst
      procedure :: delLast => refll_delLast
      procedure :: sort => refll_sort
   end type

   interface assignment(=)
      module procedure refll_copy
   end interface

contains

   subroutine refll_create(vl,maxsize)
      class(refll_vlist)    :: vl
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

   subroutine refll_copy(s,rhs)
      class(refll_vlist), intent(inout)  :: s
      type(refll_vlist), intent(in)  :: rhs
      type(refl_vlist), pointer :: rp
      integer i
      !DBG!call assert(associated(rhs%vector),"refll_copy: rhs not initialized")
      if (.not.associated(s%vector)) then 
         call refll_create(s,rhs%maxsize)
      else
         call refll_del(s)
      end if
      do i=1,rhs%cursize
         rp => rhs%vector(i)%p_data
         call refll_append(s,rp)
      end do
   end subroutine refll_copy

   subroutine refll_destroy(vl)
      class(refll_vlist)    :: vl
      call refll_del(vl)
      deallocate(vl%vector)
      vl%maxsize = 0
      vl%cursize = 0
   end subroutine

   subroutine refll_resize(vl,maxsize)
      class(refll_vlist)    :: vl
      integer, optional, intent(in) :: maxsize
      type(refll_data_pointer), pointer :: vp(:)
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
   end subroutine refll_resize

   logical function refll_isEmpty(vl)
      class(refll_vlist) :: vl
      refll_isEmpty = vl%cursize == 0
   end function refll_isEmpty

   pure integer function refll_size(vl)
      class(refll_vlist), intent(in) :: vl
      refll_size = vl%cursize
   end function refll_size

   function refll_elem(vl,i)
      class(refll_vlist) :: vl
      integer, intent(in) :: i
      type(refl_vlist), pointer :: refll_elem
      refll_elem => vl%vector(i)%p_data
   end function

   subroutine refll_append(vl,d)
      class(refll_vlist) :: vl
      type(refl_vlist), intent(in) :: d
      vl%cursize = vl%cursize + 1
      if (vl%cursize > vl%maxsize) call refll_resize(vl)
      allocate(vl%vector(vl%cursize)%p_data)
      vl%vector(vl%cursize)%p_data = d
   end subroutine refll_append

   subroutine refll_insert(vl,p,d)
      class(refll_vlist) :: vl
      integer, intent(in) :: p   ! position where to insert
      type(refl_vlist), intent(in) :: d   ! data
      integer i
      vl%cursize = vl%cursize + 1
      if (vl%cursize > vl%maxsize) call refll_resize(vl)
      do i=vl%cursize,p+1,-1
         vl%vector(i)%p_data => vl%vector(i-1)%p_data
      end do
      nullify(vl%vector(p)%p_data)
      allocate(vl%vector(p)%p_data)
      vl%vector(p)%p_data = d
   end subroutine refll_insert

   subroutine refll_del(vl,p)
      class(refll_vlist) :: vl
      integer, intent(in), optional :: p  ! position to delete
      type(refl_vlist),pointer :: dp
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
   end subroutine refll_del

   subroutine refll_delFirst(vl)
      class(refll_vlist) :: vl
      type(refl_vlist),pointer :: dp
      integer i
      dp => vl%vector(1)%p_data
      call dp%destroy()
      deallocate(dp)
      do i=1,vl%cursize-1
         vl%vector(i)%p_data => vl%vector(i+1)%p_data
      end do
      vl%cursize = vl%cursize - 1
   end subroutine refll_delFirst

   subroutine refll_delLast(vl)
      class(refll_vlist) :: vl
      type(refl_vlist),pointer :: dp
      dp => vl%vector(vl%cursize)%p_data
      call dp%destroy()
      deallocate(dp)
      vl%cursize = vl%cursize - 1
   end subroutine refll_delLast

   subroutine refll_sort(vl,cmp)
      class(refll_vlist) :: vl
      interface
         logical function cmp(i,j)
         import :: refl_vlist
         type(refl_vlist),intent(in) :: i,j
         end function cmp
      end interface
      call refll_quicksort(vl%vector(1:vl%cursize),cmp)
   end subroutine refll_sort

   recursive subroutine refll_quicksort(list,cmp)
      type(refll_data_pointer) :: list(:)
      interface
         logical function cmp(i,j)
         import :: refl_vlist
         type(refl_vlist),intent(in) :: i,j
         end function cmp
      end interface
      integer :: i,j,n
      type(refl_vlist), pointer :: chosen,temp  ! TYPE
      integer, parameter :: maxSimpleSortSize=3
      
      n = size(list)
      if (n <= maxSimpleSortSize) then
         call refll_simpleSort(list,cmp)
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
         if (j>1) call refll_quicksort(list(:j),cmp)
         if (i<n) call refll_quicksort(list(i:),cmp)
      end if
   end subroutine refll_quicksort
   
   subroutine refll_simpleSort(list,cmp)
      type(refll_data_pointer), intent(inout) :: list(:)
      interface
         logical function cmp(i,j)
         import :: refl_vlist
         type(refl_vlist),intent(in) :: i,j
         end function cmp
      end interface
      integer :: i,j
      type(refl_vlist),pointer  :: temp  
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
   end subroutine refll_simpleSort
         
end module reflistVList_m

