! Copyright (C) 2013-2014 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refStrPerm_m

! references for "structure with electron permutations for maximal coincidence" 

   use kinds_m, only: r8
   use refBase_m
   use refUtils_m, only: calcRefDifference, calcRefDiffPerm

   implicit none

   integer, parameter :: SIMPLE=0, SPINPERM=1, ALLPERM=2     ! modes for determining reference differences

   public
   private :: SIMPLE,SPINPERM,ALLPERM

   type, extends(referenceContainer) :: refContainerStructPerm   ! like Struct but electrons are permuted for max likeness
      real(r8)            :: simThreshold = 0.1d0       ! identify structures as similar if diff is smaller
      real(r8)            :: sameThreshold = 0.01d0          ! sub list thresh
   contains
      procedure :: refcsp_create
      procedure :: insert => refcsp_insert
      procedure :: writeShortList => refcsp_writeShortList
   end type refContainerStructPerm


contains


   subroutine refcsp_create(this,listLength,doSortFreq,el,simThreshold,sameThreshold)
      class(refContainerStructPerm), intent(inout) :: this
      integer, intent(in)                            :: listLength
      logical, intent(in)                            :: doSortFreq
      integer, intent(in)                            :: el
      real(r8), optional, intent(in)                   :: simThreshold
      real(r8), optional, intent(in)                   :: sameThreshold
      call refc_create(this,listLength,doSortFreq,el)
      if (present(simThreshold)) this%simThreshold = simThreshold
      if (present(sameThreshold)) this%sameThreshold = sameThreshold
   end subroutine refcsp_create


   subroutine refcsp_insert(this,r)
      class(refContainerStructPerm), intent(inout) :: this
      type(reference), intent(in)              :: r
      type(reference)                          :: r1
      type(reference), pointer                 :: rp
      type(refl_vlist)                         :: rl
      type(refl_vlist), pointer                :: p
      real(r8)                                   :: dist,dist2thr,maxdist
      integer cnt,i,j
      logical firstChanged,foundAll,foundSpin,inserted

      firstChanged = .false.
      foundAll = .false.
      dist2thr = this%sameThreshold
      call r%setCount(1)
      call rl%create(this%getElemLength()+1)
      call r1%new(ne)
      !!write(999,*) '*** s insert value=',r%f
      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         call calcRefDifference(rp,r,ALLPERM,dist,maxdist)
         !!write(999,*) '   ** ALLPERM dist=',dist,this%distanceThreshold,ALLPERM
         if (maxdist < this%simThreshold) then
            !!write(999,*) '   ** allperm struct found'
            foundAll = .true.
            foundSpin = .false.
            do j=1,p%size()
               rp => p%elem(j)
               call calcRefDifference(rp,r,SPINPERM,dist,maxdist)
               !!write(999,*) '      * SPINPERM dist=',dist,dist2thr,SPINPERM
               if (maxdist < dist2thr) then
                  !!write(999,*) '      * found: incr'
                  foundSpin = .true.
                  call rp%incrementCount()
                  exit
               end if
            end do
            if (.not.foundSpin) then
               !!write(999,*) '   ** not found in sub list'
               ! new entry for sub list: sort in list according to value
               ! permute new entry (=r1) to match first entry as close as possible
               ! if necessary implement binary search
               inserted = .false.
               r1 = r 
               call calcRefDiffPerm(rp,r1,SPINPERM)
               do j=1,p%size()
                  rp => p%elem(j)
                  if (r1%f < rp%f) then
                     ! insert here
                     !!write(999,*) '   ** inserting in sub list'
                     if (j==1) firstChanged = .true.
                     call p%insert(j,r1)
                     inserted = .true.
                     exit
                  end if
               end do
               if (.not.inserted) call p%append(r1)
               if (p%size() == this%getElemLength() + 1) then
                  ! remove last sub list entry, but add count to new last entry
                  !!write(999,*) '   ** remove last sub list entry'
                  rp => p%elem(p%size())
                  cnt = rp%count
                  rp => p%elem(p%size()-1)
                  rp%count = rp%count + cnt
                  call p%delLast()
               end if
            end if
            exit
         end if
      end do
      if (.not.foundAll) then
         ! not found: append new entry
         !!write(999,*) '   ** append to list'
         call rl%del()
         if (this%getSize()==0) then ! very first entry to list
            call rl%append(r)
         else
            r1 = r
            call this%getElemPtr(1,p)
            rp => p%elem(1)
            call calcRefDiffPerm(rp,r1,ALLPERM)
            call rl%append(r1)
         end if 
         call this%appendList(rl)
         firstChanged = .true.
      end if
      if (firstChanged) call this%sortList(hasSmallerValue)
      if (this%getSize() == this%getMaxListLength() + 1) call this%delLastList()
      call rl%destroy()
      call r1%destroy()
      !!write(999,*) '*** done'
   end subroutine refcsp_insert

   subroutine refcsp_writeShortList(this,iu,str)
      class(refContainerStructPerm)             :: this
      integer, intent(in)                   :: iu     ! unit to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist), pointer             :: p
      integer i,j,totalCount
      real(r8) f
      integer, parameter :: MAXSHOW = 15

      write(iu,'(a)') ' list contains different structures (ignoring spin) sorted w.r.t function value -ln(psi**2):'
      write(iu,'(a/)') ' each structure contains different structures (with spin) sorted w.r.t function value -ln(psi**2):'

      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         ! sum counters of all sub list elements
         totalCount = rp%count
         f = rp%f
         do j=2,p%size()
            rp => p%elem(j)
            totalCount = totalCount + rp%count
         end do
         write(iu,'(i5,a,f13.6,a,i7)') i,' structure with best value=',f,'    # found:',totalCount
         do j=1,min(p%size(),MAXSHOW)
            rp => p%elem(j)
            write(iu,'(5x,i4,2x,f13.6,2x,i7)') j,rp%f,rp%count
         end do
      end do
   end subroutine refcsp_writeShortList
   

   logical function hasSmallerValue(li,lj)
      type(refl_vlist), intent(in)  :: li,lj
      type(reference), pointer      :: p,q
      p => li%elem(1)
      q => lj%elem(1)
      hasSmallerValue = p%f <= q%f
   end function hasSmallerValue

end module refStrPerm_m
