! Copyright (C) 2012-2013, 2017 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

! this module contains various sorting routines

module sorting_m

   use kinds_m, only: r8
   use error_m
   implicit none
   private
   public :: sort, quickSortEntry,entry,entryToVector,entryFromVector, &
             quickSortIndex, inverse_permutation

   type entry
      integer :: idx=0
      integer :: tid=0
      real(r8)  :: data=0.d0
   end type

   interface sort
      module procedure double_sort
      module procedure int_sort
   end interface sort

   integer, parameter :: ENTRYDECIMAL=1000000

contains

   subroutine entryToVector(self,vec,n)
      type(entry), pointer :: self(:)
      real(r8), pointer      :: vec(:)
      integer, intent(in)  :: n      ! entry elements to convert
      integer i
      call assert(associated(self) .and. associated(vec)," entryToVector: not associtated")
      call assert(size(self)>=n .and. size(vec)>=2*n," entryToVector: sizes too small")
      do i=1,n
         vec(2*i-1) = ENTRYDECIMAL*self(i)%tid + self(i)%idx
         vec(2*i)   = self(i)%data
      end do
   end subroutine entryToVector

   subroutine entryFromVector(self,vec,n)
      type(entry), pointer :: self(:)
      real(r8), pointer      :: vec(:)
      integer, intent(in)  :: n      ! entry elements to convert
      integer i
      call assert(associated(self) .and. associated(vec)," entryToVector: not associtated")
      call assert(size(self)>=n .and. size(vec)>=2*n," entryToVector: sizes too small")
      do i=1,n
         self(i)%idx = mod(nint(vec(2*i-1)),ENTRYDECIMAL)
         self(i)%tid = nint(vec(2*i-1))/ENTRYDECIMAL
         self(i)%data = vec(2*i)
      end do
   end subroutine entryFromVector

   !----------------------------------------!
   recursive subroutine double_sort(list)
   !----------------------------------------!
   !
   !  quicksort (from Brainard F2003 book)
   !
      real(r8), intent(inout) :: list(:)     ! TYPE
      integer :: i,j,n
      real(r8)  :: chosen,temp  ! TYPE
      integer, parameter :: maxSimpleSortSize=5

      n = size(list)
      if (n <= maxSimpleSortSize) then
         call double_interchangeSort(list)
      else
         chosen = list(n/2)    ! TYPE
         i = 0
         j = n+1
         do
            do
               i = i+1
               if (list(i) >= chosen) exit   ! TYPE
            enddo
            do
               j = j-1
               if (list(j) <= chosen) exit   ! TYPE
            enddo
            if (i<j) then
               ! swap
               temp = list(i)                ! TYPE
               list(i) = list(j)
               list(j) = temp
            else if (i==j) then
               i = i+1
               exit
            else
               exit
            end if
         end do
         if (j>1) call double_sort(list(:j))
         if (i<n) call double_sort(list(i:))
      end if
   end subroutine double_sort

   !------------------------------------!
   subroutine double_interchangeSort(list)
   !------------------------------------!
      real(r8), intent(inout) :: list(:)         ! TYPE
      integer :: i,j
      real(r8)  :: temp                          ! TYPE
      do i=1,size(list)-1
         do j=i+1,size(list)
            if (list(i) > list(j)) then
               temp = list(i)                   ! TYPE
               list(i) = list(j)
               list(j) = temp
            end if
         end do
      end do
   end subroutine double_interchangeSort

   !--------------------------------!
   recursive subroutine int_sort(vec)
   !--------------------------------!
      ! same for integer vectors
      integer, intent(inout) :: vec(:)
      integer :: i,j,n
      integer :: chosen,temp 
      integer, parameter :: maxSimpleSortSize=5
      
      n = size(vec)
      if (n <= maxSimpleSortSize) then
         call int_simpleSort(vec)
      else
         chosen = vec(n/2)
         i = 0
         j = n+1
         do
            do
               i = i+1
               if (chosen <= vec(i)) exit 
            enddo 
            do 
               j = j-1
               if (vec(j) <= chosen) exit
            enddo
            if (i<j) then
               ! swap
               temp = vec(i)
               vec(i) = vec(j)
               vec(j) = temp
            else if (i==j) then
               i = i+1
               exit
            else
               exit
            end if
         end do
         if (j>1) call int_sort(vec(:j))
         if (i<n) call int_sort(vec(i:))
      end if
   end subroutine int_sort
   
   subroutine int_simpleSort(vec)
      integer, intent(inout) :: vec(:)
      integer :: i,j
      integer :: temp  
      do i=1,size(vec)-1
         do j=i+1,size(vec)
            if (vec(j) <= vec(i)) then
               ! swap
               temp = vec(i)
               vec(i) = vec(j)
               vec(j) = temp
            end if
         end do
      end do
   end subroutine int_simpleSort


   !
   !  quicksort (from Brainard book)
   !
   !---------------------------------------!
   recursive subroutine quickSortEntry(list)
   !---------------------------------------!
      type(entry), intent(inout) :: list(:)     ! TYPE
      integer :: i,j,n
      type(entry)  :: chosen,temp  ! TYPE
      integer, parameter :: maxSimpleSortSize=5

      n = size(list)
      if (n <= maxSimpleSortSize) then
         call interchangeSortEntry(list)
      else
         chosen = list(n/2)    ! TYPE
         i = 0
         j = n+1
         do
            do
               i = i+1
               if (list(i)%data >= chosen%data) exit   ! TYPE
            enddo
            do
               j = j-1
               if (list(j)%data <= chosen%data) exit   ! TYPE
            enddo
            if (i<j) then
               ! swap
               temp = list(i)                ! TYPE
               list(i) = list(j)
               list(j) = temp
            else if (i==j) then
               i = i+1
               exit
            else
               exit
            end if
         end do
         if (j>1) call quickSortEntry(list(:j))
         if (i<n) call quickSortEntry(list(i:))
      end if
   end subroutine quickSortEntry

   !-----------------------------------!
   subroutine interchangeSortEntry(list)
   !-----------------------------------!
      type(entry), intent(inout) :: list(:)         ! TYPE
      integer :: i,j
      type(entry)  :: temp                          ! TYPE
      do i=1,size(list)-1
         do j=i+1,size(list)
            if (list(i)%data > list(j)%data) then
               temp = list(i)                   ! TYPE
               list(i) = list(j)
               list(j) = temp
            end if
         end do
      end do
   end subroutine interchangeSortEntry

   !--------------------------------------------------------------------!
   recursive subroutine quickSortIndex(isGreaterEqual,isSmallerEqual,idx)
   !--------------------------------------------------------------------!
      integer, intent(inout) :: idx(:)
      integer            :: i,j,n
      integer            :: chosen,temp
      integer, parameter :: maxSimpleSortSize=5

      interface
         function isGreaterEqual(i,j) result(res)
         integer,intent(in) :: i,j
         logical            :: res
         end function isGreaterEqual
      end interface

      interface
         function isSmallerEqual(i,j) result(res)
         integer,intent(in) :: i,j
         logical            :: res
         end function isSmallerEqual
      end interface

      n = size(idx)
      if (n <= maxSimpleSortSize) then
         call interchangeSortIndex(isSmallerEqual,idx)
      else
         chosen = idx(n/2)    ! TYPE
         i = 0
         j = n+1
         do
            do
               i = i+1
               if (isGreaterEqual(idx(i),chosen)) exit
            enddo
            do
               j = j-1
               if (isSmallerEqual(idx(j),chosen)) exit
            enddo
            if (i<j) then
               ! swap
               temp = idx(i)
               idx(i) = idx(j)
               idx(j) = temp
            else if (i==j) then
               i = i+1
               exit
            else
               exit
            end if
         end do
         if (j>1) call quickSortIndex(isGreaterEqual,isSmallerEqual,idx(:j))
         if (i<n) call quickSortIndex(isGreaterEqual,isSmallerEqual,idx(i:))
      end if
   end subroutine quickSortIndex

   !-------------------------------------------------!
   subroutine interchangeSortIndex(isSmallerEqual,idx)
   !-------------------------------------------------!
      integer, intent(inout) :: idx(:)
      integer   :: i,j
      integer   :: temp

      interface
         function isSmallerEqual(i,j) result(res)
         integer,intent(in) :: i,j
         logical            :: res
         end function isSmallerEqual
      end interface

      do i=1,size(idx)-1
         do j=i+1,size(idx)
            if (.not.isSmallerEqual(idx(i),idx(j))) then
               temp = idx(i)
               idx(i) = idx(j)
               idx(j) = temp
            end if
         end do
      end do
   end subroutine interchangeSortIndex


   !-----------------------------------------------!
   recursive subroutine inverse_permutation(p, invp)
   !-----------------------------------------------!
      ! use quicksort to construct inverse permutation invp to permutation p
      ! on initial entry p must be a permuation, and invp=(1,2,...size(p))
      integer, intent(inout) :: p(:), invp(:)
      integer :: i,j,n
      integer :: chosen,temp 
      integer, parameter :: maxSimpleSortSize=4
      
      n = size(p)
      if (n <= maxSimpleSortSize) then
         call invp_simple(p, invp)
      else
         chosen = p(n/2)
         i = 0
         j = n+1
         do
            do
               i = i+1
               if (chosen <= p(i)) exit 
            enddo 
            do 
               j = j-1
               if (p(j) <= chosen) exit
            enddo
            if (i<j) then
               ! swap
               temp = p(i)
               p(i) = p(j)
               p(j) = temp
               temp    = invp(i)
               invp(i) = invp(j)
               invp(j) = temp
            else if (i==j) then
               i = i+1
               exit
            else
               exit
            end if
         end do
         if (j>1) call inverse_permutation(p(:j), invp(:j))
         if (i<n) call inverse_permutation(p(i:), invp(i:))
      end if
   end subroutine inverse_permutation
   
   subroutine invp_simple(p, invp)
      integer, intent(inout) :: p(:), invp(:)
      integer :: i,j
      integer :: temp  
      do i=1,size(p)-1
         do j=i+1,size(p)
            if (p(j) <= p(i)) then
               ! swap
               temp = p(i)
               p(i) = p(j)
               p(j) = temp
               temp    = invp(i)
               invp(i) = invp(j)
               invp(j) = temp
            end if
         end do
      end do
   end subroutine invp_simple





end module sorting_m
