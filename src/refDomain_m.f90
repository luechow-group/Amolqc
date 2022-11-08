! Copyright (C) 2014, 2017 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refDomain_m

   ! Like the centroid type collecting maxima by comparison with
   ! given references (centroid .ctr file or in situ averaged from .ref file)
   ! Additionally a new function returns the  permutation index of the
   ! last inserted reference for domain sampling.
   ! No averaging for reference centroids.

   use kinds_m, only: r8
   use refBase_m
   use refUtils_m, only: calcRefDifference, calcRefDiffPerm
   use sorting_m, only: sort

   implicit none

   integer, parameter :: SIMPLE=0, SPINPERM=1, ALLPERM=2     ! modes for determining reference differences

   public
   private :: SIMPLE,SPINPERM,ALLPERM

   !!! could be derived from refContainerCtr. Only getLastMaxPermIdx is added, but create and insert are modified

   type, extends(referenceContainer) :: refContainerDom      ! like Struct but electrons are permuted for max likeness and averaged
                                                             ! the first element in each sublist contains the centroid reference (if available)
      real(r8)                :: simThreshold = 0.1d0          ! identify structures as similar if max dist is smaller
      real(r8)                :: sameThreshold = 0.01d0           ! sub list thresh, i.e. for identical structures (w/ spin)
      integer, pointer      :: excl(:) => null()             ! list of excluded positions, see "str" mode
      integer               :: exclMode = 0
      type(matrixstat),allocatable  :: ctrPos(:)             ! statistics calculating the centroid from maxima for each reference
      integer               :: lastm = 0                     ! save the number of the last maximum found at insert
      integer, pointer      :: lastpermidx(:) => null()      ! save the corresponding permutation index
   contains
      procedure :: refcd_create
      procedure :: destroy => refcd_destroy
      procedure :: getFirstF => refcd_getFirstF
      procedure :: insert => refcd_insert
      procedure :: getLastMaxPermIdx => refcd_getLastMaxPermIdx
      procedure :: calcDistances => refcd_calcDistances
      procedure :: writeShortList => refcd_writeShortList
      procedure :: writeFullList => refcd_writeFullList
   end type refContainerDom


contains


   subroutine refcd_create(this,listLength,doSortFreq,el,refFile,simThreshold,sameThreshold,exclFile)
      class(refContainerDom), intent(inout)        :: this
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

      call assert(MASTER,"refcd_create: must be called by MASTER only")
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

      allocate(this%lastpermidx(ne))
      this%lastpermidx(1:ne) = (/ (i,i=1,ne) /)

      if (present(refFile)) then
         nlen = len(trim(refFile))
         suffix = refFile(nlen-3:nlen)
         if (suffix == ".ref") then
            call read_refFile()            ! inner subroutine
         else if (suffix == ".ctr") then
            call read_ctrFile()            ! inner subroutine
         else
            call abortp("max_mode=dom requires ref_file with .ref or .ctr ending (1)")
         end if
      else
         call abortp("max_mode=dom requires ref_file with .ref or .ctr ending (2)")
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
         if (io /= 0) call abortp("refcd_create: refFile does not exist")
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
            call assert(noel==ne,"refcd_create: wrong number of electrons in ref file")
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
         if (io /= 0) call abortp("refcd_create: ctrFile does not exist")
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
            call assert(noel==ne,"refcd_create: wrong number of electrons in ref file")
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

   end subroutine refcd_create


   subroutine refcd_destroy(this)
      class(refContainerDom), intent(inout) :: this
      integer l
      do l=1,this%getSize()
         call this%ctrPos(l)%destroy()
      end do
      deallocate(this%ctrPos)
      call refc_destroy(this)
      deallocate(this%lastpermidx)
   end subroutine refcd_destroy


   real(r8) function refcd_getFirstF(this)
      class(refContainerDom), intent(inout) :: this
      type(reference), pointer                 :: rp
      type(refl_vlist), pointer                :: p
      call this%getElemPtr(1,p)
      if (p%size() > 1) then
         rp => p%elem(2)
      else
         rp => p%elem(1)
      end if
      refcd_getFirstF = rp%f
   end function refcd_getFirstF



   subroutine refcd_insert(this,r)
      class(refContainerDom), intent(inout) :: this
      type(reference), intent(in)              :: r
      type(reference), pointer                 :: rp
      type(refl_vlist), pointer                :: p
      real(r8)                                   :: maxdist,meandist,smallestDist
      integer i,smallestIdx

      ! refDomain insert is not really inserting the reference
      ! Only the closest of the stored references is identified, the index and
      ! the permutation is stored for later retrieval
      ! The reference is counted, though

      !!!write(999,*) '*** s insert value=',r%f
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
      !!!write(999,*) '           smallest:',smallestIdx,smallestDist
      if (smallestDist < this%simThreshold) then
         !!!write(999,*) '   ** allperm struct found'

         ! save the idx
         this%lastm = smallestIdx

         ! get permutation index of the structure
         call this%getElemPtr(smallestIdx,p)
         rp => p%elem(1)
         call calcRefDifference(rp,r,ALLPERM,meandist,maxdist,permidx=this%lastpermidx)
         !!!write(999,*) '          found smallestIdx:',smallestIdx
         call rp%incrementCount()
      else
         this%lastm = 0
      end if
      !!!write(999,*) '*** done'
   end subroutine refcd_insert


   subroutine refcd_getLastMaxPermIdx(this,m,pidx)
      ! returns the maximum number m and the permutation index pidx saved
      ! at the last insert call
      class(refContainerDom), intent(inout) :: this
      integer, intent(inout)                :: m
      integer, intent(inout)                :: pidx(:)
      call assert(size(pidx)==ne,"refcd_getLastPermIdx: illegal size of pidx")
      m = this%lastm
      pidx = this%lastpermidx
   end subroutine refcd_getLastMaxPermIdx


   subroutine refcd_calcDistances(this,r,mode,maxIdx,outIdx,meanDist,maxDist,idxList)
      ! find smallest maxDist (and return meanDist as well) of r to any of the FIRST list entries (n,1)
      class(refContainerDom), intent(inout)   :: this
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
         call assert(size(idxList) >= n,"refcd_calcDistances: idxList is not sufficently long")
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

   end subroutine refcd_calcDistances


   subroutine refcd_writeShortList(this,iu,str)
      class(refContainerDom)             :: this
      integer, intent(in)                   :: iu     ! unit to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist), pointer             :: p
      integer i,totalCount

      write(iu,'(a)') ' List contains different structures (when ignoring spin) determined with centroid reference.'
      totalCount = this%getCount()
      write(iu,'(a,i12/)') ' total number of maxima collected with these structure:',totalCount
      write(iu,'(a)') '   distribution into the reference structures:'
      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         write(iu,'(5x,i4,i10)') i,rp%count
      end do

   end subroutine refcd_writeShortList


   subroutine refcd_writeFullList(this,fname,str)
      class(refContainerDom), intent(inout) :: this
      character(len=*), intent(in)          :: fname  ! file to write to (without .ref!!!)
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist), pointer             :: p
      integer                               :: naNuc, nbNuc    ! alpha, beta elecs at nucleus (not varied)
      integer i,n,io
      integer, parameter :: iu2 = 13
      real(r8), parameter  :: thresh = 0.002d0

      open(iu2,file=trim(fname)//'.dom',iostat=io)
      if (io/=0) call abortp('(references_write): opening file '//trim(fname)//'.dom failed')
      write(iu2,'(i5,a)') ncenter, " nuclei:"
      do i=1,getNNuc()
         write(iu2,'(i4,1x,a2,1x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
      enddo

      write(iu2,'(i5,1x,a,1x,i10,i6)') this%getSize(),"ref/ctr",this%getCount(),1

      do n=1,this%getSize()
         call this%getElemPtr(n,p)
         rp => p%elem(1)
         call findNucElecs(thresh, rp%x, rp%y, rp%z, naNuc, nbNuc)
         write(iu2,'(a,i5,a,f15.6,a,i10,4i4)') trim(str)//":",n," 1   F("//trim(str)//"):", &
         rp%f, "  found:", rp%count, naNuc, getNAlpha(), nbNuc, getNelec()
         write(iu2,'(i5)') size(rp%x)
         do i=1,size(rp%x)
            write(iu2,'(3f13.6)') rp%x(i),rp%y(i),rp%z(i)
         end do
      end do
      close(iu2)
   end subroutine refcd_writeFullList

end module refDomain_m
