! Copyright (C) 2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refVector_m

   ! this module defines a class for vectors of references (refADT_m)
   ! purpose: compare given electron arrangements with the current list
   ! and return permutation index and vector index.
   !
   ! parallelism: this container is locally duplicated, no communication needed

   use kinds_m, only: r8
   use global_m, only: getNElec,iul,MASTER,ne,logmode
   use wfData_m
   use utils_m, only: tokenize
   use refADT_m
   use refUtils_m, only: calcRefDifference, calcRefDiffPerm
   use sorting_m, only: sort
   use newStatistics_m
   use mpiInterface_m, only: myMPIBcastInteger, myMPIBcastDouble

   implicit none

   integer, parameter :: SIMPLE=0, SPINPERM=1, ALLPERM=2     ! modes for determining reference differences

   public
   private :: SIMPLE,SPINPERM,ALLPERM

   type referenceVector
      ! wrapper type to encapsulate vector of references data structure
      private
      type(reference), pointer :: rv(:)
      integer                  :: nElems
      integer, pointer         :: excl(:) => null()             ! list of excluded positions, see "str" mode
      integer                  :: exclMode = 0
      real(r8)                   :: simThreshold = 0.1d0          ! identify structures as similar if max dist is smaller
      type(matrixstat), allocatable  :: ctrPos(:)            ! statistics calculating the centroid from maxima for each reference
      character(len=120)       :: infoMessage1 = "    no exclude list used"
      character(len=120)       :: infoMessage2 = "    no initial references provided"
      type(atom),pointer       :: mol(:) => null()   ! molecule geometry if different from wf
      logical                  :: doSpaceWarpTransform = .true.  ! if references are read from file
   contains
      procedure :: create => refv_create
      procedure :: destroy => refv_destroy
      procedure :: getSize => refv_getSize
      procedure :: getInfoMessage => refv_getInfoMessage
      procedure :: setInfoMessage => refv_setInfoMessage
      procedure :: getRefPtr => refv_getRefPtr
      procedure :: setRef => refv_setRef
      procedure :: calcDistances => refv_calcDistances
      procedure :: getIdxAndPermutation => refv_getIdxAndPermutation
   end type referenceVector

contains

   subroutine refv_create(this,listLength,refFile,simThreshold,exclFile)
      class(referenceVector), intent(inout)        :: this
      integer, intent(in)                          :: listLength
      character(len=40), intent(in)                :: refFile
      real(r8), optional, intent(in)                 :: simThreshold
      character(len=40), optional, intent(in)      :: exclFile
      real(r8) x(ne),y(ne),z(ne),xx(ne,3)
      real(r8), allocatable :: vec(:)
      integer l,i,j,n,nWords,io,nlen
      character(len=120) line
      character(len=20) words(12)
      character(len=1) emode
      character(len=4) suffix
      type(reference) :: r,r0
      real(r8) nucThresh,coreThresh,bondDiff
      integer, parameter :: iu = 12

      allocate(this%rv(listLength))
      do i=1,listLength
         call this%rv(i)%new(ne)
      enddo
      this%nElems = 0

      if (present(simThreshold)) this%simThreshold = simThreshold

      ! read file containing coded list of excluded position (as determined by atoms_m.f90:atoms_whatPosition)
      if (present(exclFile)) then
         if (MASTER) then
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
            call this%setInfoMessage(1,line)
         endif
         call myMPIBcastInteger(n,1)
         if (.not.MASTER) allocate(this%excl(n))
         call myMPIBcastInteger(this%exclMode,1)
         call myMPIBcastInteger(this%excl,n)
      endif

      allocate(this%ctrPos(listLength))
      do l=1,listLength
         call this%ctrPos(l)%create(ne,3)
      end do

      ! read references for assignment from refFile and broadcast to all nodes
      nlen = len(trim(refFile))
      suffix = refFile(nlen-3:nlen)
      if (suffix == ".ref") then
         if (MASTER) then
            call read_refFile()            ! inner subroutine
         endif
         call broadcast_rv()               ! inner subroutine
      else if (suffix == ".ctr") then
         if (MASTER) then
            call read_ctrFile()            ! inner subroutine
         endif
         call broadcast_rv()               ! inner subroutine
      else
         call abortp("references or centroid file (.ref or .ctr) required")
      end if

      do l=1,listLength
         call this%ctrPos(l)%destroy()
      end do
      deallocate(this%ctrPos)

   contains

      subroutine read_refFile()
         integer ref,ref2,count,nomax,noel,io,outIdx,tCount,elemLength,idx
         real(r8) wgt,f,meanDist,maxDist
         character(len=20) str
         integer, parameter :: iu=12, iu1=13, iu3=14


         open(iu,file=refFile,status='old',iostat=io)
         if (io /= 0) call abortp("refv_create: refFile does not exist")
         read(iu,*) n   ! ncenter
         do i=1,n
            read(iu,*) j  ! index,elem,x,y,z
         end do

         ! accumulate refs (n,j=1,jmax) using ALLPERM permutation to 1st into ctrPos matrix stat
         do i=1,listLength
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
         if (ref > listLength) exit
            read(iu,*) noel
            call assert(noel==ne,"refv_create: wrong number of electrons in ref file")
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
         ref = min(ref,listLength)

         ! construct initial references from accumulated read references, if distances are large enough
         xx = this%ctrPos(1)%localmean()
         call r%set(xx(:,1),xx(:,2),xx(:,3),0.d0)


         open(iu3,file='CTR-assignment.txt')
         write(iu3,'(a17)') 'Ref #   |  CTR # '

         this%rv(1) = r
         idx = 1
         write(iu3,'(i8,a,i8)') 1,'|',1
         do i=2,ref
            xx = this%ctrPos(i)%localmean()
            call r%set(xx(:,1),xx(:,2),xx(:,3),0.d0)
            call this%calcDistances(r,ALLPERM,idx,outIdx,meanDist,maxDist)
            if (maxDist >= this%simThreshold) then
               idx = idx + 1
               this%rv(idx) = r
            end if
            write(iu3,'(i8,a,i8)') i,'|',idx
         end do

         close(iu3)

         this%nElems = idx

         do i=1,ref
            call this%ctrPos(i)%reset()
         end do

         call writeRefFile(this)
         call r%destroy()
         call r0%destroy()

         write(line,'(i5,2a)') this%nElems," centroid references created from ref file ",trim(refFile)
         call this%setInfoMessage(2,line)

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
                  this%rv(idx) = r
                  idx = idx + 1
               end if
               iref = iref + 1
            endif
         if (this%getSize()==listLength .or. iref > nomax) exit
         enddo
         close(iu)
         call r%destroy()
         this%nElems = idx

         write(line,'(i5,2a)') this%nElems," centroid references read from ctr file ",trim(refFile)
         call this%setInfoMessage(2,line)

      end subroutine read_ctrFile

      subroutine broadcast_rv()
         integer idx,l
         allocate(vec(3*ne*listlength+1))
         if (MASTER) then
            idx = 0
            do l=1,listlength
               vec(idx+1:idx+ne)        = this%rv(l)%x(1:ne)
               vec(idx+ne+1:idx+2*ne)   = this%rv(l)%y(1:ne)
               vec(idx+2*ne+1:idx+3*ne) = this%rv(l)%z(1:ne)
               idx = idx + 3*ne
            enddo
            vec(idx+1) = this%nElems
         endif
         call myMPIBcastDouble(vec,3*ne*listLength+1)
         if (.not.MASTER) then
            idx = 0
            do l=1,listlength
               this%rv(l)%x(1:ne)     = vec(idx+1:idx+ne)
               this%rv(l)%y(1:ne)     = vec(idx+ne+1:idx+2*ne)
               this%rv(l)%z(1:ne)     = vec(idx+2*ne+1:idx+3*ne)
               idx = idx + 3*ne
            enddo
            this%nElems = nint(vec(idx+1))
         endif
         deallocate(vec)
      end subroutine broadcast_rv

   end subroutine refv_create


   subroutine refv_destroy(this)
      class(referenceVector), intent(inout) :: this
      integer i
      do i=1,size(this%rv)
         call this%rv(i)%destroy()
      enddo
      deallocate(this%rv)
      if (associated(this%mol)) deallocate(this%mol)
   end subroutine refv_destroy

   function refv_getSize(this) result(vecsize)
      class(referenceVector), intent(in) :: this
      integer                            :: vecsize
      vecsize = size(this%rv)
   end function refv_getSize

   subroutine refv_getInfoMessage(this,n,s)
      class(referenceVector), intent(in) :: this
      integer, intent(in)                :: n
      character(len=120), intent(out)    :: s
      if (n==1) then
         s = this%infoMessage1
      else if (n==2) then
         s = this%infoMessage2
      endif
   end subroutine refv_getInfoMessage

   subroutine refv_setInfoMessage(this,n,s)
      class(referenceVector), intent(inout) :: this
      integer, intent(in)                   :: n
      character(len=120), intent(in)        :: s
      if (n==1) then
         this%infoMessage1 = s
      else if (n==2) then
         this%infoMessage2 = s
      endif
   end subroutine refv_setInfoMessage

   subroutine refv_getRefPtr(this,i,rptr)
      class(referenceVector), intent(inout)  :: this
      integer, intent(in)                    :: i
      type(reference), pointer, intent(out)  :: rptr
      rptr => this%rv(i)
   end subroutine refv_getRefPtr

   subroutine refv_setRef(this,i,rptr)
      class(referenceVector), intent(inout)  :: this
      integer, intent(in)                    :: i
      type(reference), pointer, intent(out)  :: rptr
      this%rv(i) = rptr
   end subroutine refv_setRef


   subroutine refv_calcDistances(this,r,mode,maxIdx,outIdx,meanDist,maxDist,idxList)
      ! find smallest maxDist (and return meanDist as well) of r to any of the FIRST list entries (n,1)
      class(referenceVector), intent(inout)      :: this
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
         call this%getRefPtr(idx(i),rp)
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

   end subroutine refv_calcDistances


   subroutine refv_getIdxAndPermutation(this,r,n,permIdx)
      class(referenceVector), intent(inout) :: this
      type(reference), intent(in)           :: r          ! reference to check
      integer, intent(inout)                :: n          ! return idx of closest reference from this%rv
      integer, intent(inout)                :: permIdx(:) ! return permutation of closest reference
      type(reference), pointer              :: rp
      real(r8)                                :: maxdist,meandist,smallestDist
      integer i,smallestIdx

      !!write(900+mytid,*) '***   getIdxAndPermutation:'
      call r%setCount(1)

      smallestDist = 999999.d0
      smallestIdx = 1
      do i=1,this%nElems
         call this%getRefPtr(i,rp)
         if (associated(this%excl)) then
            call calcRefDifference(rp,r,ALLPERM,meandist,maxdist,this%excl,exclmode=this%exclMode)
         else
            call calcRefDifference(rp,r,ALLPERM,meandist,maxdist)
         end if
         !!write(900+mytid,*) '   ** ALLPERM maxdist:',i,maxdist,this%simThreshold,ALLPERM
         if (maxdist < smallestDist) then
            smallestDist = maxdist
            smallestIdx = i
         end if
      end do
      !!write(900+mytid,*) '           smallest:',smallestIdx,smallestDist
      if (smallestDist < this%simThreshold) then
         !!write(900+mytid,*) '   ** allperm struct found'
         n = smallestIdx
         ! get permutation of the structure
         call this%getRefPtr(smallestIdx,rp)
         call calcRefDifference(rp,r,ALLPERM,meandist,maxdist,permidx=permIdx)
         !!write(900+mytid,*) '          found smallestIdx:',smallestIdx
      else
         n = 0
         permIdx = (/ (i,i=1,n) /)
      end if
      !!write(900+mytid,*) '*** done'

   end subroutine refv_getIdxAndPermutation

   subroutine writeRefFile(this)
      class(referenceVector), intent(inout)        :: this
      type(reference), pointer              :: rp
      integer i,n,io
      integer, parameter :: iu2 = 13
      real(r8), parameter  :: thresh = 0.002d0

      open(iu2,file=trim(basename)//'.ref',iostat=io)
      if (io/=0) call abortp('(rwriteRefFile): opening file '//trim(basename)//'.ref failed')
      write(iu2,'(i5,a)') ncenter, " nuclei:"
      do i=1,getNNuc()
         write(iu2,'(i4,1x,a2,1x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
      enddo

      write(iu2,'(i5,1x,a,1x,i10,i6)') this%nElems,"ref/ctr",0,0

      do n=1,this%nElems
         call this%getRefPtr(n,rp)
         write(iu2,'(a,i5,a,f15.6,a,i10,4i4)') "CTR:",n," 1   F(CTR):", &
         rp%f,"  found:",rp%count,0,0,0,0
         write(iu2,'(i5)') size(rp%x)
         do i=1,size(rp%x)
            write(iu2,'(3f13.6)') rp%x(i),rp%y(i),rp%z(i)
         end do
      end do
      close(iu2)

   end subroutine writeRefFile

end module refVector_m
