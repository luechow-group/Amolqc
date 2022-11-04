! Copyright (C) 2013-2014, 2017 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refCtr_m

   ! "centroid"-type references: calc centroids for positions (indep of spin)

   use kinds_m, only: r8
   use refBase_m
   use refUtils_m, only: calcRefDifference, calcRefDiffPerm
   use sorting_m, only: sort,quickSortIndex


   implicit none

   integer, parameter :: SIMPLE=0, SPINPERM=1, ALLPERM=2     ! modes for determining reference differences

   public
   private :: SIMPLE,SPINPERM,ALLPERM

!!! sollte eigentlich refContainerStr als Basisklasse haben. Teste  !!

   type, extends(referenceContainer) :: refContainerCtr      ! like Struct but electrons are permuted for max likeness and averaged
                                                             ! the first element in each sublist contains the centroid reference (if available)
      real(r8)                :: simThreshold = 0.1d0          ! identify structures as similar if max dist is smaller
      real(r8)                :: sameThreshold = 0.01d0           ! sub list thresh, i.e. for identical structures (w/ spin)
      integer, pointer      :: excl(:) => null()             ! list of excluded positions, see "str" mode
      integer               :: exclMode = 0
      type(matrixstat),allocatable  :: ctrPos(:)             ! statistics calculating the centroid from maxima for each reference
   contains
      procedure :: refcc_create
      procedure :: destroy => refcc_destroy
      procedure :: getFirstF => refcc_getFirstF
      procedure :: insert => refcc_insert
      procedure :: calcDistances => refcc_calcDistances
      procedure :: writeShortList => refcc_writeShortList
      procedure :: writeFullList => refcc_writeFullList
   end type refContainerCtr


contains


   subroutine refcc_create(this,listLength,doSortFreq,el,refFile,simThreshold,sameThreshold,exclFile)
      class(refContainerCtr), intent(inout)        :: this
      integer, intent(in)                          :: listLength
      logical, intent(in)                          :: doSortFreq
      integer, intent(in)                          :: el
      character(len=40), optional, intent(in)      :: refFile
      real(r8), optional, intent(in)                 :: simThreshold
      real(r8), optional, intent(in)                 :: sameThreshold
      character(len=40), optional, intent(in)      :: exclFile
      real(r8) x(ne),y(ne),z(ne),xx(ne,3)
      integer l,i,j,n,nWords,io,nlen
      character(len=120) line
      character(len=20) words(12)
      character(len=1) emode
      character(len=4) suffix
      type(reference) :: r,r0
      real(r8) nucThresh,coreThresh,bondDiff
      integer, parameter :: iu = 12

      call assert(MASTER,"refcc_create: must be called by MASTER only")
      call refc_create(this,listLength,doSortFreq,el+1)        ! additional element for centroid data
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

      allocate(this%ctrPos(listLength))
      do l=1,listLength
         call this%ctrPos(l)%create(ne,3)
      end do

      if (present(refFile)) then
         nlen = len(trim(refFile))
         suffix = refFile(nlen-3:nlen)
         if (suffix == ".ref") then
            call read_refFile()            ! inner subroutine
         else if (suffix == ".ctr") then
            call read_ctrFile()            ! inner subroutine
         else
            call abortp("max_mode=ctr requires ref_file with .ref or .ctr ending (1)")
         end if
      else
         call abortp("max_mode=ctr requires ref_file with .ref or .ctr ending (2)")
      end if

   contains

      subroutine read_refFile()
         integer ref,ref2,count,nomax,noel,io,outIdx,tCount,elemLength,idx
         real(r8) wgt,f,meanDist,maxDist
         type(reference), pointer                 :: rp
         type(refl_vlist), pointer                :: p
         character(len=20) str
         integer, parameter :: iu=12, iu1=13

         open(iu,file=refFile,status='old',iostat=io)
         if (io /= 0) call abortp("refcc_create: refFile does not exist")
         read(iu,*) n   ! ncenter
         do i=1,n
            read(iu,*) j  ! index,elem,x,y,z
         end do

         ! accumulate refs (n,j=1,jmax) using ALLPERM permutation to 1st into ctrPos matrix stat
         do i=1,this%getMaxListLength()
            call this%ctrPos(i)%reset()
         end do
         read(iu,*) nomax,str,tCount,elemLength
         ref = 0

         call r%new(ne)
         call r0%new(ne)
         do
            read(iu,'(A)',iostat=io) line
         if (io /= 0) exit
            call tokenize(line,words,n)
            read(words(2),*) ref
            read(words(3),*) ref2
            read(words(5),*) f
            read(words(7),*) count
         if (ref > this%getMaxListLength()) exit
            read(iu,*) noel
            call assert(noel==ne,"refcc_create: wrong number of electrons in ref file")
            do i=1,noel
               read(iu,*) xx(i,1),xx(i,2),xx(i,3)
            enddo
            wgt = count
            if (ref2==1) then
               call r0%set(xx(:,1),xx(:,2),xx(:,3),f)
               call this%ctrPos(ref)%add(xx,wgt)
            else if (ref2 < elemLength) then ! skip last because of count accumulation
               call r%set(xx(:,1),xx(:,2),xx(:,3),f)
               call calcRefDiffPerm(r0,r,ALLPERM)
               call r%get(xx(:,1),xx(:,2),xx(:,3),f)
               call this%ctrPos(ref)%add(xx,wgt)
            end if
         enddo
         close(iu)
         ref = min(ref,this%getMaxListLength())

         ! construct initial references from accumulated read references, if distances are large enough
         xx = this%ctrPos(1)%localmean()
         call r%set(xx(:,1),xx(:,2),xx(:,3),0.d0)
         call this%append(r)
         idx = 1
         do i=2,ref
            xx = this%ctrPos(i)%localmean()
            call r%set(xx(:,1),xx(:,2),xx(:,3),0.d0)
            call this%calcDistances(r,ALLPERM,idx,outIdx,meanDist,maxDist)
            if (maxDist >= this%simThreshold) then
               idx = idx + 1
               call this%append(r)
            end if
         end do

         do i=1,ref
            call this%ctrPos(i)%reset()
         end do

         open(iu1,file="init.ctr")
         write(iu1,*) this%getSize()
         do i=1,this%getSize()
            call this%getElemPtr(i,p)
            rp => p%elem(1)
            call rp%get(xx(:,1),xx(:,2),xx(:,3),f)
            write(iu1,*) i
            do j=1,ne
               write(iu1,'(3g15.5)') xx(j,1), xx(j,2), xx(j,3)
            end do
         end do
         close(iu1)
         call r%destroy()
         call r0%destroy()

         write(line,'(i5,2a)') this%getSize()," centroid references created from ref file ",trim(refFile)
         call this%setInfoMessage(line)

      end subroutine read_refFile


      subroutine read_ctrFile()
         integer ref,ref2,count,nomax,noel,iref,outIdx,idx
         real(r8) f,meanDist,maxDist
         open(iu,file=refFile,status='old',iostat=io)
         if (io /= 0) call abortp("refcc_create: ctrFile does not exist")
         read(iu,*) n   ! ncenter
         do i=1,n
            read(iu,*) j  ! index,elem,x,y,z
         end do

         call r%new(ne)
         read(iu,*) nomax

         iref = 1
         idx = 0
         do
            read(iu,'(A)',iostat=io) line
         if (io /= 0) exit
            call tokenize(line,words,n)
            read(words(2),*) ref
            read(words(3),*) ref2
            read(words(5),*) f
            read(words(7),*) count
            read(iu,*) noel
            call assert(noel==ne,"refcc_create: wrong number of electrons in ref file")
            do i=1,noel
               read(iu,*) x(i),y(i),z(i)
            enddo
            if (ref == iref .and. ref2==1) then
               call r%set(x,y,z,f)
               call this%calcDistances(r,ALLPERM,idx,outIdx,meanDist,maxDist)
               if (maxDist >= this%simThreshold) then
                  idx = idx + 1
                  call this%append(r)
               end if
               iref = iref + 1
            endif
         if (this%getSize()==listLength .or. iref > nomax) exit
         enddo
         close(iu)
         call r%destroy()

         write(line,'(i5,2a)') this%getSize()," centroid references read from ctr file ",trim(refFile)
         call this%setInfoMessage(line)

      end subroutine read_ctrFile

   end subroutine refcc_create


   subroutine refcc_destroy(this)
      class(refContainerCtr), intent(inout) :: this
      integer l
      do l=1,this%getSize()
         call this%ctrPos(l)%destroy()
      end do
      deallocate(this%ctrPos)
      call refc_destroy(this)
   end subroutine refcc_destroy


   real(r8) function refcc_getFirstF(this)
      class(refContainerCtr), intent(inout) :: this
      type(reference), pointer                 :: rp
      type(refl_vlist), pointer                :: p
      call this%getElemPtr(1,p)
      if (p%size() > 1) then
         rp => p%elem(2)
      else
         rp => p%elem(1)
      end if
      refcc_getFirstF = rp%f
   end function refcc_getFirstF



   subroutine refcc_insert(this,r)
      class(refContainerCtr), intent(inout) :: this
      type(reference), intent(in)              :: r
      type(reference)                          :: r1
      type(reference), pointer                 :: rp
      type(refl_vlist)                         :: rl
      type(refl_vlist), pointer                :: p
      real(r8)                                   :: x(ne,3)
      real(r8)                                   :: f,maxdist,meandist,smallestDist
      logical foundSpin, inserted
      integer i,j,smallestIdx,cnt

      !!write(999,*) '*** s insert value=',r%f
      call r%setCount(1)


      smallestDist = 999999.d0
      smallestIdx = 1
      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         if (associated(this%excl)) then
            call calcRefDifference(rp,r,ALLPERM,meandist,maxdist,exclist=this%excl,exclmode=this%exclMode)
         else
            call calcRefDifference(rp,r,ALLPERM,meandist,maxdist)
         end if
         !!write(999,*) '   ** ALLPERM maxdist:',i,maxdist,this%simThreshold,ALLPERM
         if (maxdist < smallestDist) then
            smallestDist = maxdist
            smallestIdx = i
         end if
      end do
      !!write(999,*) '           smallest:',smallestIdx,smallestDist
      if (smallestDist < this%simThreshold) then
         !!write(999,*) '   ** allperm struct found'

         ! first add electron positions after permutation to centroid statistics
         call r1%new(ne)
         r1 = r      ! permute only a copy
         call this%getElemPtr(smallestIdx,p)
         rp => p%elem(1)
         call calcRefDiffPerm(rp,r1,ALLPERM)
         call r1%get(x(:,1),x(:,2),x(:,3),f)
         call this%ctrPos(smallestIdx)%add(x,1.d0)
         !!write(999,*) '          adding to ctrPos:',smallestIdx

         call r1%destroy()

         ! secondly, add electron positions into "str" reference sublist (code like in str)
         ! remember: 1st entry in sublist is the (old) centroid position
         foundSpin = .false.
         do j=2,p%size()
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
            do j=2,p%size()
               rp => p%elem(j)
               if (r%f < rp%f) then
                  ! insert here
                  !!write(999,*) '   ** inserting in sub list'
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

      else
         ! not found: append new entry: to "str" list, to "ctrRef". Add to ctrPos statistics
         !!write(999,*) '   ** append to list'
         if (this%getSize() < this%getMaxListLength()) then
            call rl%create(this%getElemLength()+1)
            !!write(999,*) '       ',rl%size()
            call rl%append(r)      ! 1st head entry (instead of centroid)
            call rl%append(r)      ! 2nd real entry
            call this%appendList(rl)
            call r%get(x(:,1),x(:,2),x(:,3),f)
            i = this%getSize()
            call this%ctrPos(i)%add(x,1.d0)
            call rl%destroy()
         else
            !!write(999,*), '       not appended max list length reached: size=',this%getSize()
         end if
      end if
      !!write(999,*) '*** done'
   end subroutine refcc_insert

   subroutine refcc_calcDistances(this,r,mode,maxIdx,outIdx,meanDist,maxDist,idxList)
      ! find smallest maxDist (and return meanDist as well) of r to any of the FIRST list entries (n,1)
      class(refContainerCtr), intent(inout)   :: this
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
         call assert(size(idxList) >= n,"refcc_calcDistances: idxList is not sufficently long")
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

   end subroutine refcc_calcDistances


   subroutine refcc_writeShortList(this,iu,str)
      class(refContainerCtr)             :: this
      integer, intent(in)                   :: iu     ! unit to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist), pointer             :: p
      integer i,j,totalCount,outIdx
      integer, allocatable                  :: idx(:)
      real(r8) f,meanDist,maxDist
      integer, parameter :: MAXSHOW = 15

      write(iu,'(a)') ' List contains different structures (when ignoring spin) determined with centroid reference.'
      write(iu,'(a)') ' Each structure contains list of similar structures (with spin) sorted w.r.t function value -ln(psi**2):'
      totalCount = this%getCount()
      write(iu,'(a,i12/)') ' total number of maxima collected with these structure:',totalCount

      ! sort list according to totalCount of sublists
      allocate(idx(this%getSize()))
      idx = [ (i,i=1,size(idx)) ]
      call quickSortIndex(isSmallerEqual,isGreaterEqual,idx)

      do i=1,this%getSize()
         call this%getElemPtr(idx(i),p)
         ! sum counters of all sub list elements
         f = 0
         totalCount = 0
         do j=2,p%size()
            rp => p%elem(j)
            if (j==2) f = rp%f
            totalCount = totalCount + rp%count
         end do
         if (i>1) then
            rp => p%elem(1)    ! compare with head element == reference structure
            call this%calcDistances(rp,ALLPERM,i-1,outIdx,meanDist,maxDist,idx)
            write(iu,'(i5,a,f13.6,a,i7,a,i5,a,2f10.3)') i,' structure with best value=',f,'    # found:',totalCount, &
                                          "   min dist to ",outIDx,"  max/mean dist:",maxDist*bohr2ang,meanDist*bohr2ang
         else
            write(iu,'(i5,a,f13.6,a,i7)') i,' structure with best value=',f,'    # found:',totalCount
         end if
         do j=2,min(p%size(),MAXSHOW)
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
         do k=2,p%size()
            rp => p%elem(k)
            counti = counti + rp%count
         end do
         call this%getElemPtr(j,p)
         countj = 0
         do k=2,p%size()
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
         do k=2,p%size()
            rp => p%elem(k)
            counti = counti + rp%count
         end do
         call this%getElemPtr(j,p)
         countj = 0
         do k=2,p%size()
            rp => p%elem(k)
            countj = countj + rp%count
         end do
         isSmallerEqual = counti <= countj
      end function isSmallerEqual


   end subroutine refcc_writeShortList


   subroutine refcc_writeFullList(this,fname,str)
      class(refContainerCtr), intent(inout) :: this
      character(len=*), intent(in)          :: fname  ! file to write to (without .ref!!!)
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist), pointer             :: p
      integer                               :: naNuc, nbNuc    ! alpha, beta elecs at nucleus (not varied)
      integer i,n,m,io
      integer(i8) totalcnt
      integer, allocatable                  :: idx(:)
      real(r8) xx(ne,3),sg(ne,3)
      integer, parameter :: iu1 = 12
      integer, parameter :: iu2 = 13
      real(r8), parameter  :: thresh = 0.002d0

      open(iu1,file=trim(fname)//'.ref',iostat=io)
      if (io/=0) call abortp('(references_write): opening file '//trim(fname)//'.ref failed')
      open(iu2,file=trim(fname)//'.ctr',iostat=io)
      if (io/=0) call abortp('(references_write): opening file '//trim(fname)//'.ctr failed')
      write(iu1,'(i5,a)') ncenter, " nuclei:"
      write(iu2,'(i5,a)') ncenter, " nuclei:"
      do i=1,getNNuc()
         write(iu1,'(i4,1x,a2,1x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
         write(iu2,'(i4,1x,a2,1x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
      enddo

      write(iu1,'(i5,1x,a,1x,i9,i6)') this%getSize(),trim(str),this%getCount(),this%getElemLength()-1
      write(iu2,'(i5,1x,a,1x,i9,i6)') this%getSize(),"centroids",this%getCount(),1

      ! sort list according to totalCount of sublists
      allocate(idx(this%getSize()))
      idx = [ (i,i=1,size(idx)) ]
      call quickSortIndex(isSmallerEqual,isGreaterEqual,idx)

      do n=1,this%getSize()
         call this%getElemPtr(idx(n),p)
         do m=2,p%size()
            rp => p%elem(m)
            call findNucElecs(thresh, rp%x, rp%y, rp%z, naNuc, nbNuc)
            write(iu1,'(a,2i5,a,f15.6,a,i6,4i4)') trim(str)//":",n,m-1," F("//trim(str)//"):", &
               rp%f, "  found:", rp%count, naNuc, getNAlpha(), nbNuc, getNelec()
            write(iu1,'(i5)') size(rp%x)
            do i=1,size(rp%x)
               write(iu1,'(3f13.6)') rp%x(i),rp%y(i),rp%z(i)
            end do
         end do
         totalcnt = this%ctrPos(idx(n))%localcount()
         write(iu2,'(a,i5,a,i6,a)') "CTR:",n," 1  F(CTR):  0.0 found: ",totalcnt," 0 0 0 0"
         write(iu2,'(i5)') getNelec()
         if (totalcnt > 0) then
            xx = this%ctrPos(idx(n))%localmean()
         else
            xx = 0
         end if
         if (totalcnt > 1) then
            sg = this%ctrPos(idx(n))%localstddev()
         else
            sg = 0
         end if
         do i=1,getNelec()
            write(iu2,'(3f13.6,2x,3f13.6)') xx(i,1),xx(i,2),xx(i,3),sg(i,1),sg(i,2),sg(i,3)
         end do
      end do
      close(iu1)
      close(iu2)

      deallocate(idx)
   contains

      logical function isGreaterEqual(i,j)
         integer, intent(in) :: i,j
         integer counti,countj,k
         ! sum counters of all sub list elements and compare
         call this%getElemPtr(i,p)
         counti = 0
         do k=2,p%size()
            rp => p%elem(k)
            counti = counti + rp%count
         end do
         call this%getElemPtr(j,p)
         countj = 0
         do k=2,p%size()
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
         do k=2,p%size()
            rp => p%elem(k)
            counti = counti + rp%count
         end do
         call this%getElemPtr(j,p)
         countj = 0
         do k=2,p%size()
            rp => p%elem(k)
            countj = countj + rp%count
         end do
         isSmallerEqual = counti <= countj
      end function isSmallerEqual



   end subroutine refcc_writeFullList

end module refCtr_m
