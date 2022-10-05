! Copyright (C) 2013-2014 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refStr_m

   ! "structure"-type references

   use kinds_m, only: r8
   use refBase_m
   use sorting_m, only: sort, quickSortIndex
   use utils_m, only: tokenize
   use atom_m, only: atoms_initPositionThresholds

   implicit none

   integer, parameter :: SIMPLE=0, SPINPERM=1, ALLPERM=2     ! modes for determining reference differences

   public
   private :: SIMPLE,SPINPERM,ALLPERM

   type, extends(referenceContainer) :: refContainerStruct
      real(r8)            :: simThreshold = 0.1d0             ! identify structures as similar if diff is smaller
      real(r8)            :: sameThreshold = 0.01d0           ! sub list thresh, i.e. for identical structures (w/ spin)
      integer, pointer  :: excl(:) => null()
      integer           :: exclMode = 0
   contains
      procedure :: refcs_create
      procedure :: destroy => refcs_destroy
      procedure :: insert => refcs_insert
      procedure :: calcDistances => refcs_calcDistances
      procedure :: writeShortList => refcs_writeShortList
   end type refContainerStruct


contains


   subroutine refcs_create(this,listLength,doSortFreq,el,exclFile,simThreshold,sameThreshold)
      class(refContainerStruct), intent(inout) :: this
      integer, intent(in)                            :: listLength
      logical, intent(in)                            :: doSortFreq
      integer, intent(in)                            :: el
      character(len=40), optional, intent(in)        :: exclFile
      real(r8), optional, intent(in)                   :: simThreshold
      real(r8), optional, intent(in)                   :: sameThreshold
      integer, parameter                             :: iu = 12
      character(len=120)                             :: line
      character(len=20)                              :: words(20)
      character(len=1)                               :: emode
      integer io,i,n,nWords
      real(r8) nucThresh,coreThresh,bondDiff

      call assert(MASTER,"refcs_create: must be called by MASTER only")
      call refc_create(this,listLength,doSortFreq,el)
      if (present(simThreshold)) this%simThreshold = simThreshold
      if (present(sameThreshold)) this%sameThreshold = sameThreshold
      if (present(exclFile)) then
         ! read file containing coded list of excluded position (as determined by atoms_m.f90:atoms_whatPosition)
         open(iu,file=exclFile,status='old',iostat=io)
         call assert(io==0,' exclusion file '//trim(exclFile)//' does not exist')
         read(iu,'(a)') line
         call tokenize(line,words,nWords)
         if (nWords==0) call abortp("refcs_create: illegal format in exclusion file")
         read(words(1),*) n
         if (nWords>=2) then
            read(words(2),'(a1)') emode
            if (emode=='i') this%exclMode = 1
         end if
         if (nWords >= 5) then
            read(words(3),*) nucThresh
            read(words(4),*) coreThresh
            read(words(5),*) bondDiff
            call atoms_initPositionThresholds(nucThresh,coreThresh,bondDiff)
         end if
         allocate(this%excl(n))
         do i=1,n
            read(iu,*) this%excl(i)
         end do
         close(iu)
         call sort(this%excl)
         write(line,'(a,a,a,i1)') 'Exclusion File ',trim(exclFile),' read. Exclusion Mode: ',this%exclMode
         call this%setInfoMessage(line)
      end if
   end subroutine refcs_create


   subroutine refcs_destroy(this)
      class(refContainerStruct), intent(inout) :: this
      if (associated(this%excl)) deallocate(this%excl)
      call refc_destroy(this)
   end subroutine refcs_destroy


   subroutine refcs_insert(this,r)
      class(refContainerStruct), intent(inout) :: this
      type(reference), intent(in)              :: r
      type(reference), pointer                 :: rp
      type(refl_vlist)                         :: rl
      type(refl_vlist), pointer                :: p
      real(r8)                                   :: meandist,maxdist
      integer cnt,i,j
      logical firstChanged,foundAll,foundSpin,inserted

      firstChanged = .false.
      foundAll = .false.
      call r%setCount(1)
      call rl%create(this%getElemLength()+1)
      !!write(999,*) '*** s insert value=',r%f
      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         if (associated(this%excl)) then
            call calcRefDifference(rp,r,ALLPERM,meandist,maxdist,exclist=this%excl,exclmode=this%exclMode)
         else
            call calcRefDifference(rp,r,ALLPERM,meandist,maxdist)
         end if
         !!write(999,*) '   ** ALLPERM maxdist=',maxdist,this%simThreshold,ALLPERM
         if (maxdist < this%simThreshold) then
            !!write(999,*) '   ** allperm struct found'
            foundAll = .true.
            foundSpin = .false.
            do j=1,p%size()
               rp => p%elem(j)
               if (associated(this%excl)) then
                  call calcRefDifference(rp,r,SPINPERM,meandist,maxdist,exclist=this%excl,exclmode=this%exclMode)
               else
                  call calcRefDifference(rp,r,SPINPERM,meandist,maxdist)
               end if
               !!write(999,*) '      * SPINPERM maxdist=',maxdist,this%sameThreshold,SPINPERM
               if (maxdist < this%sameThreshold) then
                  !!write(999,*) '      * found: incr'
                  foundSpin = .true.
                  call rp%incrementCount()
                  exit
               end if
            end do
            if (.not.foundSpin) then
               !!write(999,*) '   ** not found in sub list'
               ! new entry for sub list: sort in list according to value
               ! if necessary implement binary search
               inserted = .false.
               do j=1,p%size()
                  rp => p%elem(j)
                  if (r%f < rp%f) then
                     ! insert here
                     !!write(999,*) '   ** inserting in sub list'
                     if (j==1) firstChanged = .true.
                     call p%insert(j,r)
                     inserted = .true.
                     exit
                  end if
               end do
               if (.not.inserted) call p%append(r)
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
         !!write(999,*) '       ',rl%size()
         call rl%append(r)
         call this%appendList(rl)
         firstChanged = .true.
      end if
      if (firstChanged) call this%sortList(hasSmallerValue)
      if (this%getSize() == this%getMaxListLength() + 1) call this%delLastList()
      call rl%destroy()
      !!write(999,*) '*** done'
   end subroutine refcs_insert

   subroutine refcs_calcDistances(this,r,mode,maxIdx,outIdx,meanDist,maxDist,idxList)
      ! find smallest maxDist (and return meanDist as well) of r to any of the FIRST list entries (n,1)
      class(refContainerStruct), intent(inout)   :: this
      type(reference), intent(in)                :: r
      integer, intent(in)                        :: mode    !  0:SIMPLE, 1:SPINPERM, 2:ALLPERM
      integer, intent(in)                        :: maxIdx   ! 0: all, >0 compare up to maxIdx
      integer, intent(out)                       :: outIdx
      real(r8), intent(out)                        :: meanDist ! mean distance
      real(r8), intent(out)                        :: maxDist  ! maximum distance
      integer, intent(in), optional              :: idxList(:)
      integer smallestIdx,i,n
      integer, allocatable :: idx(:)
      real(r8) smallestDist,smallestMean,dmean,dmax
      type(reference), pointer                 :: rp
      type(refl_vlist), pointer                :: p

      smallestDist = 999999.d0
      smallestMean = 999999.d0
      smallestIdx = 1

      n = this%getSize()
      if (maxIdx > 0) n = min(n,maxIdx)

      allocate(idx(n))
      if (present(idxList)) then
         call assert(size(idxList) >= n,"refcs_calcDistances: idxList is not sufficently long")
         idx = idxList(1:n)
      else
         idx = (/ (i,i=1,n) /)
      end if

      do i=1,n
         call this%getElemPtr(idx(i),p)
         rp => p%elem(1)
         if (associated(this%excl)) then
            call calcRefDifference(rp,r,mode,dmean,dmax,exclist=this%excl,exclmode=this%exclMode)
         else
            call calcRefDifference(rp,r,mode,dmean,dmax)
         end if
         if (dmax < smallestDist) then
            smallestDist = dmax
            smallestMean = dmean
            smallestIdx = i
         end if
      end do
      outIdx = smallestIdx
      meanDist = smallestMean
      maxDist = smallestDist
      deallocate(idx)
   end subroutine refcs_calcDistances


   subroutine refcs_writeShortList(this,iu,str)
      class(refContainerStruct)             :: this
      integer, intent(in)                   :: iu     ! unit to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist), pointer             :: p
      integer, allocatable                  :: idx(:)
      integer i,j,totalCount,outIdx
      real(r8) f,meanDist,maxDist
      integer, parameter :: MAXSHOW = 15

      if (this%sortForFreqs()) then
         write(iu,'(a)') ' List contains different structures (ignoring spin) sorted w.r.t frequency:'
      else
         write(iu,'(a)') ' List contains different structures (ignoring spin) sorted w.r.t function value -ln(psi**2):'
      endif
      write(iu,'(a)') ' each structure contains a list of similar structures (with spin) sorted w.r.t function value -ln(psi**2):'
      totalCount = this%getCount()
      write(iu,'(a,i12/)') ' total number of maxima collected with these structure:',totalCount

      ! sort list according to totalCount of sublists
      allocate(idx(this%getSize()))
      idx = [ (i,i=1,size(idx)) ]
      if (this%sortForFreqs()) then
         call quickSortIndex(isSmallerEqual,isGreaterEqual,idx)
      endif

      do i=1,this%getSize()
         call this%getElemPtr(idx(i),p)
         rp => p%elem(1)
         ! sum counters of all sub list elements
         totalCount = rp%count
         f = rp%f
         do j=2,p%size()
            rp => p%elem(j)
            totalCount = totalCount + rp%count
         end do
         if (i>1) then
            call this%calcDistances(rp,ALLPERM,i-1,outIdx,meanDist,maxDist,idx)
            write(iu,'(i5,a,f13.6,a,i7,a,i5,a,2f10.3)') i,' structure with best value=',f,'    # found:',totalCount, &
                                          "   min dist to ",outIDx,"  max/mean dist:",maxDist*bohr2angs,meanDist*bohr2angs
         else
            write(iu,'(i5,a,f13.6,a,i7)') i,' structure with best value=',f,'    # found:',totalCount
         end if
         do j=1,min(p%size(),MAXSHOW)
            rp => p%elem(j)
            write(iu,'(5x,i4,2x,f13.6,2x,i7)') j,rp%f,rp%count
         end do
      end do

      deallocate(idx)

   contains

      logical function isGreaterEqual(i,j)
         integer, intent(in) :: i,j
         integer counti,countj,k
         ! sum counters of all sub list elements and compare
         call this%getElemPtr(i,p)
         counti = 0
         do k=1,p%size()
            rp => p%elem(k)
            counti = counti + rp%count
         end do
         call this%getElemPtr(j,p)
         countj = 0
         do k=1,p%size()
            rp => p%elem(k)
            countj = countj + rp%count
         end do
         isGreaterEqual = counti >= countj
      end function isGreaterEqual

      logical function isSmallerEqual(i,j)
         integer, intent(in) :: i,j
         integer counti,countj,k
         ! sum counters of all sub list elements and compare
         call this%getElemPtr(i,p)
         counti = 0
         do k=1,p%size()
            rp => p%elem(k)
            counti = counti + rp%count
         end do
         call this%getElemPtr(j,p)
         countj = 0
         do k=1,p%size()
            rp => p%elem(k)
            countj = countj + rp%count
         end do
         isSmallerEqual = counti <= countj
      end function isSmallerEqual

   end subroutine refcs_writeShortList


   logical function hasSmallerValue(li,lj)
      type(refl_vlist), intent(in)  :: li,lj
      type(reference), pointer      :: p,q
      p => li%elem(1)
      q => lj%elem(1)
      hasSmallerValue = p%f <= q%f
   end function hasSmallerValue


end module refStr_m
