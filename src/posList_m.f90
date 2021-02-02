! Copyright (C) 2019 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later
!
!
! vlist template (vector list)
! template list based on array of pointers 
!
! insert/append insert *copies* of the data into the vlist
! therefore a copy constructor (and overloaded assignment) is required for the data type
! if deallocation of the data type is required set   appropriately to calling the destructor for pointer variable dp
!
! m4 -D_MODNAME_=IntVListModule -D_PRE_=int_ -D_TYPE_=integer -D_DA_= -D_USEMOD_= -D_IMP_= vlist_tmpl_m.m4 > int_vlist_m.f90
! m4 -D_MODNAME_=RefListModule -D_PRE_=refl_ -D_TYPE_=type(reference) -D_DA_="call dp%destroy()" -D_USEMOD_="use referenceADT" -D_IMP_="import reference"  vlist_tmpl_m.m4 > ref_vlist_m.f90
!
! careful: sort (with quick sort) requires in external "isSmaller" logical function 
! for increasing order (or "isGreater" for decreasing order)
! This function *must* have <= or >= semantics!!
!
module posList_m

   use kinds_m, only: r8   
   implicit none

   type posVal_t
      real(r8) :: pos(3) = 0.0_r8
      real(r8) :: value = 0.0_r8
      integer  :: count = 0
   end type posVal_t

   type dataPointer_t
      private
      type(posVal_t), pointer :: p_data
   end type dataPointer_t

   type pos_VList_t
      private
      type(dataPointer_t), pointer :: vector(:)
      integer          :: cursize
      integer          :: maxSize
   contains
      procedure :: create => pos_create
      procedure :: destroy => pos_destroy
      procedure :: resize => pos_resize
      procedure :: size => pos_size
      procedure :: elem => pos_elem
      procedure :: append => pos_append
      procedure :: insert => pos_insert
      procedure :: del => pos_del
      procedure :: sort => pos_sort
   end type pos_VList_t

contains

   subroutine pos_create(vl,maxsize)
      class(pos_VList_t)    :: vl
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

   subroutine pos_destroy(vl)
      class(pos_VList_t)    :: vl
      call pos_del(vl)
      deallocate(vl%vector)
      vl%maxsize = 0
      vl%cursize = 0
   end subroutine

   subroutine pos_resize(vl,maxsize)
      class(pos_VList_t)    :: vl
      integer, optional, intent(in) :: maxsize
      type(dataPointer_t), pointer :: vp(:)
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
   end subroutine

   integer function pos_size(vl)
      class(pos_VList_t) :: vl
      pos_size = vl%cursize
   end function pos_size

   function pos_elem(vl,i)
      class(pos_VList_t) :: vl
      integer, intent(in) :: i
      type(posVal_t), pointer :: pos_elem
      pos_elem => vl%vector(i)%p_data
   end function

   subroutine pos_append(vl,d)
      class(pos_VList_t) :: vl
      type(posVal_t), intent(in) :: d
      vl%cursize = vl%cursize + 1
      if (vl%cursize > vl%maxsize) call pos_resize(vl)
      allocate(vl%vector(vl%cursize)%p_data)
      vl%vector(vl%cursize)%p_data = d
   end subroutine pos_append

   subroutine pos_insert(vl,p,d)
      class(pos_VList_t) :: vl
      integer, intent(in) :: p   ! position where to insert
      type(posVal_t), intent(in) :: d   ! data
      integer i
      vl%cursize = vl%cursize + 1
      if (vl%cursize > vl%maxsize) call pos_resize(vl)
      do i=vl%cursize,p+1,-1
         vl%vector(i)%p_data => vl%vector(i-1)%p_data
      end do
      nullify(vl%vector(p)%p_data)
      allocate(vl%vector(p)%p_data)
      vl%vector(p)%p_data = d
   end subroutine pos_insert

   subroutine pos_del(vl,p)
      class(pos_VList_t) :: vl
      integer, intent(in), optional :: p  ! position to delete
      type(posVal_t),pointer :: dp
      integer i
      if (present(p)) then
         dp => vl%vector(p)%p_data
         
         deallocate(dp)
         do i=p,vl%cursize-1
            vl%vector(i)%p_data => vl%vector(i+1)%p_data
         end do
         vl%cursize = vl%cursize - 1
      else
         do i=1,vl%cursize
            dp => vl%vector(i)%p_data
            
            deallocate(dp)
         end do
         vl%cursize = 0
      end if
   end subroutine pos_del

   subroutine pos_sort(vl,cmp)
      class(pos_VList_t) :: vl
      interface
         logical function cmp(i,j)
            import posVal_t
            type(posVal_t),intent(in) :: i,j
         end function cmp
      end interface
      call pos_quicksort(vl%vector(1:vl%cursize),cmp)
   end subroutine pos_sort

   recursive subroutine pos_quicksort(list,cmp)
      type(dataPointer_t) :: list(:)
      interface
         logical function cmp(i,j)
            import posVal_t         
            type(posVal_t),intent(in) :: i,j
         end function cmp
      end interface
      integer :: i,j,n
      type(posVal_t), pointer :: chosen,temp  ! TYPE
      integer, parameter :: maxSimpleSortSize=3
      
      n = size(list)
      if (n <= maxSimpleSortSize) then
         call pos_simpleSort(list,cmp)
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
         if (j>1) call pos_quicksort(list(:j),cmp)
         if (i<n) call pos_quicksort(list(i:),cmp)
      end if
   end subroutine pos_quicksort
   
   subroutine pos_simpleSort(list,cmp)
      type(dataPointer_t), intent(inout) :: list(:)
      interface
         logical function cmp(i,j)
            import posVal_t         
            type(posVal_t),intent(in) :: i,j
         end function cmp
      end interface
      integer :: i,j
      type(posVal_t),pointer  :: temp  
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
   end subroutine pos_simpleSort
         
   logical function isSmaller(pvi, pvj) 
      type(posVal_t), intent(in) :: pvi, pvj
      isSmaller = pvi%count <= pvj%count
   end function isSmaller

   logical function isGreater(pvi, pvj) 
      type(posVal_t), intent(in) :: pvi, pvj
      isGreater = pvi%count >= pvj%count
   end function isGreater

end module posList_m

