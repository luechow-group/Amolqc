! Copyright (C) 2013-2014, 2017 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refBase_m

   ! this module creates an abstract base class (referenceContainer) and a "create" method
   ! that creates the actual class. Only the base class methods need to be used in applications
   ! (polymorphic objects)
   ! The basis data structure is a list of list of reference electron configurations, requiring
   ! two indices for access
   !
   ! parallelism: this container should be kept on the MASTER (only)
   !   collect data (GATHER) on the MASTER and use insert/append for all gathered data
   !   no combining/merging is done in this module! use local statistics only!

   use kinds_m, only: r8
   use global_m, only: getNElec,iul,MASTER,ne,logmode
   use wfData_m
   use findNucElecs_m, only: findNucElecs
   use utils_m, only: tokenize
   use refADT_m
   use refVList_m
   use reflistVList_m
   use refUtils_m, only: calcRefDifference
   use sorting_m, only: quickSortIndex
   use newStatistics_m
   use atom_m, only: atom

   implicit none

   public

   type, abstract :: referenceContainer
      ! wrapper type to encapsulate List of List data structure
      private
      type(refll_vlist)  :: rll
      integer            :: maxListLength = 0
      integer            :: maxElemLength = 1        ! one element in list element only
      character(len=240) :: infoMessage = "no initial references provided"
      type(atom),pointer :: mol(:) => null()   ! molecule geometry if different from wf
      logical            :: doSpaceWarpTransform = .true.  ! if references are read from file
      logical            :: sortFreq = .false.               ! sort wrt frequency of maxima or function value (-ln(psi^2))?
   contains
      procedure :: create => refc_create
      procedure :: destroy => refc_destroy
      procedure :: getSize => refc_getSize
      procedure :: getCount => refc_getCount
      procedure :: getMaxListLength => refc_getMaxListLength
      procedure :: getInfoMessage => refc_getInfoMessage
      procedure :: setInfoMessage => refc_setInfoMessage
      procedure :: getFirstF => refc_getFirstF
      procedure :: getElemLength => refc_getElemLength
      procedure :: getElemPtr => refc_getElemPtr
      procedure :: insertList => refc_insertList
      procedure :: appendList => refc_appendList
      procedure :: delLastList => refc_delLastList
      procedure :: sortList => refc_sortList
      procedure :: sortForFreqs => refc_sortForFreqs
      procedure :: append => refc_append
      procedure (refc_insert), deferred :: insert
      procedure :: calcDistances => refc_calcDistances
      procedure :: readFromFile => refc_readFromFile
      procedure :: appendRefFromFile => refc_appendRefFromFile
      procedure (refc_writeShortList), deferred :: writeShortList
      procedure :: writeSimpleList => refc_writeSimpleList
      procedure :: writeFullList => refc_writeFullList
   end type referenceContainer

   abstract interface
      subroutine refc_insert(this,r)
         import :: referenceContainer,reference
         class(referenceContainer), intent(inout) :: this
         type(reference), intent(in)              :: r
      end subroutine refc_insert
   end interface

   abstract interface
      subroutine refc_writeShortList(this,iu,str)
         import :: referenceContainer,reference
         class(referenceContainer)             :: this
         integer, intent(in)                   :: iu     ! unit to write to
         character(len=*), intent(in)          :: str    ! output string e.g. "Ref" or "Max"
      end subroutine refc_writeShortList
   end interface

contains

   subroutine refc_create(this,listLength,doSortFreq,elemLength)
      class(referenceContainer), intent(inout) :: this
      integer, intent(in)                      :: listLength
      logical, intent(in)                      :: doSortFreq
      integer, optional, intent(in)            :: elemLength
      call this%rll%create(listLength+1)
      this%maxListLength = listLength
      this%sortFreq = doSortFreq
      if (present(elemLength)) this%maxElemLength = elemLength
   end subroutine refc_create

   subroutine refc_destroy(this)
      class(referenceContainer), intent(inout) :: this
      call this%rll%destroy()
      if (associated(this%mol)) deallocate(this%mol)
   end subroutine refc_destroy

   integer function refc_getSize(this,i)
      class(referenceContainer), intent(in) :: this
      integer, optional, intent(in)         :: i       ! if given, return size of i-th sublist
      type(refl_vlist), pointer                :: p
      if (present(i)) then
         if (i <= this%rll%size()) then
            p => this%rll%elem(i)
            refc_getSize = p%size()
         else
            refc_getSize = 0
         end if
      else
         refc_getSize = this%rll%size()
      end if
   end function refc_getSize

   integer function refc_getCount(this,i,j)
      class(referenceContainer), intent(in) :: this
      integer, optional, intent(in)         :: i       ! if given, return all counts of i-th sublist
      integer, optional, intent(in)         :: j       ! if given, return count of j-th element in i-th sublist
                                                       ! if neither given, return total count
      type(refl_vlist), pointer             :: p
      type(reference), pointer              :: rp
      integer totalCount,k,l

      totalCount = 0
      if (present(i)) then
         if (i <= this%rll%size()) then
            p => this%rll%elem(i)
            if (present(j)) then
               if (j <= p%size()) then
                  rp => p%elem(j)
                  totalCount = rp%getCount()
               end if
            else
               do k=1,p%size()
                  rp => p%elem(k)
                  totalCount = totalCount + rp%getCount()
               end do
            end if
         end if
      else
         do l=1,this%rll%size()
            p => this%rll%elem(l)
            do k=1,p%size()
               rp => p%elem(k)
               totalCount = totalCount + rp%getCount()
            end do
         end do
      end if
      refc_getCount = totalCount
   end function refc_getCount


   pure integer function refc_getMaxListLength(this)
      class(referenceContainer), intent(in) :: this
      refc_getMaxListLength = this%maxListLength
   end function refc_getMaxListLength


   subroutine refc_getInfoMessage(this,s)
      class(referenceContainer), intent(in) :: this
      character(len=120), intent(out)        :: s
      s = this%infoMessage
   end subroutine refc_getInfoMessage

   subroutine refc_setInfoMessage(this,s)
      class(referenceContainer), intent(inout) :: this
      character(len=120), intent(in)            :: s
      if (this%infoMessage == "no initial references provided") then
         this%infoMessage = s
      else
         write(this%infoMessage,'(A)') trim(this%infoMessage)//NEW_LINE('A')//trim(s)
      endif
   end subroutine refc_setInfoMessage


   real(r8) function refc_getFirstF(this)
      class(referenceContainer), intent(inout) :: this
      type(reference), pointer                 :: rp
      type(refl_vlist), pointer                :: p
      if (this%rll%size() == 0) then
         !if reference list is empty, return NaN as first f value
         refc_getFirstF = 99999.99999d0
      else
        p => this%rll%elem(1)
        rp => p%elem(1)
        refc_getFirstF = rp%f
      endif
   end function refc_getFirstF

   pure integer function refc_getElemLength(this)
      class(referenceContainer), intent(in) :: this
      refc_getElemLength = this%maxElemLength
   end function refc_getElemLength

   subroutine refc_getElemPtr(this,i,ptr)
      class(referenceContainer), intent(inout) :: this
      integer, intent(in)                      :: i
      type(refl_vlist), pointer                :: ptr
      ptr => this%rll%elem(i)
   end subroutine refc_getElemPtr

   subroutine refc_insertList(this,i,rl)
      class(referenceContainer), intent(inout) :: this
      integer, intent(in)                      :: i
      type(refl_vlist),intent(in)              :: rl
      call this%rll%insert(i,rl)
   end subroutine refc_insertList

   subroutine refc_appendList(this,rl)
      class(referenceContainer), intent(inout) :: this
      type(refl_vlist),intent(in)              :: rl
      call this%rll%append(rl)
   end subroutine refc_appendList

   subroutine refc_delLastList(this)
      class(referenceContainer), intent(inout) :: this
      call this%rll%delLast()
   end subroutine refc_delLastList

   subroutine refc_sortList(this,cmp)
      class(referenceContainer), intent(inout) :: this
      interface
         logical function cmp(i,j)
         import :: refl_vlist
         type(refl_vlist),intent(in) :: i,j
         end function cmp
      end interface
      call this%rll%sort(cmp)
   end subroutine refc_sortList

   logical function refc_sortForFreqs(this)
      class(referenceContainer), intent(inout) :: this
      refc_sortForFreqs = this%sortFreq
   end function

   subroutine refc_append(this,r)
      class(referenceContainer), intent(inout) :: this
      type(reference), intent(in)              :: r
      type(refl_vlist)                         :: rl
      call rl%create(this%maxElemLength+1)
      call rl%append(r)
      call this%rll%append(rl)
      call rl%destroy()
   end subroutine refc_append

   subroutine refc_calcDistances(this,r,mode,maxIdx,outIdx,meanDist,maxDist,idxList)
      ! find smallest maxDist (and return meanDist as well) of r to any of the FIRST list entries (n,1)
      class(referenceContainer), intent(inout)   :: this
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
         call assert(size(idxList) >= n,"refc_calcDistances: idxList is not sufficently long")
         idx = idxList(1:n)
      else
         idx = (/ (i,i=1,n) /)
      end if

      do i=1,n
         call this%getElemPtr(idx(i),p)
         rp => p%elem(1)
         call calcRefDifference(rp,r,mode,dmean,dmax)
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
   end subroutine refc_calcDistances

   subroutine refc_appendRefFromFile(this,f,count,lines,nl,idx)
      class(referenceContainer), intent(inout) :: this
      real(r8), intent(in)                  :: f
      integer, intent(in)                 :: count
      character(len=*), intent(in)        :: lines(:)
      integer, intent(in)                 :: nl        ! true # of lines
      integer, intent(inout)              :: idx       ! start reading line idx
                                                       ! return next line
      type(refl_vlist)                    :: rl
      type(reference)                     :: r
      real(r8)                              :: x(getNElec()),y(getNElec()),z(getNElec())
      integer                             :: nElec,i,n
      character(len=20)                   :: words(10)

      read(lines(idx),*) nElec
      if (nElec/=getNElec()) call abortp('reading ref file: incorrect # elecs ')
      do i=1,nElec
         read(lines(idx+i),*) x(i),y(i),z(i)
      end do
      if (this%doSpaceWarpTransform) then
         call atoms_spaceWarpTransform(atoms,this%mol,x,y,z)
      end if
      if (.not. all(abs(x)<huge(1.d0)) .or. .not. all(abs(y)<huge(1.d0)) .or. .not. all(abs(z)<huge(1.d0))) &
         call abortp("reading ref file: not all data are numbers")
      idx = idx + nElec + 1
      call r%new(getNElec())
      call r%set(x,y,z,f)
      call r%setCount(count)

      call rl%create(this%maxElemLength+1)
      call rl%append(r)
      if (this%getSize() < this%maxListLength) call this%rll%append(rl)
      call rl%destroy()
      call r%destroy()
   end subroutine refc_appendRefFromFile

   subroutine refc_readFromFile(this,refs,lines,nl,doSpaceWarpTransform)
      class(referenceContainer), intent(inout) :: this
      integer, intent(in)                      :: refs(2)    ! 0: read first refs only; n: read all refs n m=1,...refs(2)
      character(len=60), intent(in)            :: lines(:)
      integer, intent(in)                      :: nl
      logical, intent(in)                      :: doSpaceWarpTransform
      integer i,j,idx,nRef,nElec,n,nw,nn,mm,iRef,maxRef,count
      real(r8) f
      character(len=120)                       :: line
      character(len=40)                        :: words(10)
      logical startReading,stopReading

      read(lines(1),*) n   ! ncenter
      if (associated(this%mol)) deallocate(this%mol)
      allocate(this%mol(n))
      do i=1,n
         read(lines(i+1),*) j,this%mol(i)%elem, this%mol(j)%cx , this%mol(j)%cy , this%mol(j)%cz   ! index,elem,x,y,z
      end do
      this%doSpaceWarpTransform = doSpaceWarpTransform
      read(lines(n+2),*) nRef
      nElec = getNElec()
      call assert(nElec>0,'reading ref file: # of elecs not set ($wf required!)')
      idx = n + 3
      select case (refs(1))
      case (0)
         maxRef = min(refs(2),nRef)
         call this%create(maxRef,.false.)
         iRef = 0
         do
            call tokenize(lines(idx),words,nw)
            if (nw<7) call abortp('reading ref file: ref format incorrect')
            read(words(2),*) nn
            read(words(3),*) mm
            read(words(5),*) f
            read(words(7),*) count
            if (mm == 1) then
               idx = idx + 1
               call refc_appendRefFromFile(this,f,count,lines,nl,idx)
               iRef = iRef + 1
            else
               idx = idx + nElec + 2
            end if
            if (iRef==nRef) exit
         end do
      case (1:)
         call this%create(refs(2),.false.)
         iRef = 0
         startReading = .false.
         stopReading = .false.
         do
            call tokenize(lines(idx),words,nw)
            if (nw==1) exit
            if (nw<7) call abortp('reading ref file: ref format incorrect')
            read(words(2),*) nn
            read(words(3),*) mm
            read(words(5),*) f
            read(words(7),*) count
            if (nn == refs(1)) then
               startReading = .true.
               idx = idx + 1
               call refc_appendRefFromFile(this,f,count,lines,nl,idx)
               iRef = iRef + 1
            else
               if (startReading) stopReading = .true.
               idx = idx + nElec + 2
            end if
            if (iRef==refs(2) .or. stopReading) exit
         end do
      case default
         call abortp("references-readFromFile: illegal refs")
      end select
   end subroutine refc_readFromFile


   subroutine refc_writeSimpleList(this,fname,str)
      class(referenceContainer)             :: this
      character(len=*), intent(in)          :: fname  ! file to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist)                      :: rl
      type(refl_vlist), pointer             :: p
      integer                               :: naNuc, nbNuc    ! alpha, beta elecs at nucleus (not varied)
      integer, allocatable                  :: idx(:)
      integer n,i,j,io,totalCount
      integer, parameter            :: iu = 13
      real(r8) f,thresh

      !! careful: usually only from from MASTER

      thresh = 0.002d0     ! could be input parameter

      ! sort list according to totalCount of sublists
      allocate(idx(this%getSize()))
      idx = [ (i,i=1,size(idx)) ]
      if (this%sortFreq) call quickSortIndex(isSmallerEqual,isGreaterEqual,idx)

      open(iu,file=trim(fname)//'.ref',iostat=io)
      if (io/=0) call abortp('(references_write): opening file '//trim(fname)//' failed')
      write(iu,'(i5,a)') ncenter, " nuclei:"
      do i=1,getNNuc()
         write(iu,'(i4,1X,a2,1X,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
      enddo

      write(iu,'(i5,1x,a)') this%getSize(),trim(str)

      do n=1,this%rll%size()
         p => this%rll%elem(idx(n))
         rp => p%elem(1)
         ! sum counters of all sub list elements
         totalCount = rp%count
         f = rp%f
         call findNucElecs(thresh, rp%x, rp%y, rp%z, naNuc, nbNuc)
         do j=2,p%size()
            rp => p%elem(j)
            totalCount = totalCount + rp%count
         end do
         write(iu,'(a,i5,a,f15.5,a,i6,4i4)') trim(str)//":",n,"  1 F("//trim(str)//"):", &
            f, "  found:", totalCount, naNuc, getNAlpha(), nbNuc, getNelec()
         write(iu,'(i5)') size(rp%x)
         do i=1,size(rp%x)
            write(iu,'(3f13.6)') rp%x(i),rp%y(i),rp%z(i)
         end do
      end do
      close(iu)
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

   end subroutine refc_writeSimpleList

   subroutine refc_writeFullList(this,fname,str)
      class(referenceContainer), intent(inout) :: this
      character(len=*), intent(in)          :: fname  ! file to write to (without .ref!!!)
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist)                      :: rl
      type(refl_vlist), pointer             :: p
      integer                               :: naNuc, nbNuc    ! alpha, beta elecs at nucleus (not varied)
      integer                               :: aNucElec(getNNuc()), bNucElec(getNNuc()) ! nuc a -> elec i (at nuc)
      integer, allocatable                  :: idx(:)
      integer i,n,m,io
      real(r8) thresh
      character(len=3) nstr
      integer, parameter            :: iu = 13

      thresh = 0.002d0     ! could be input parameter

      ! sort list according to totalCount of sublists
      allocate(idx(this%getSize()))
      idx = (/ (i,i=1,size(idx)) /)
      if (this%sortFreq) call quickSortIndex(isSmallerEqual,isGreaterEqual,idx)

      open(iu,file=trim(fname)//'.ref',iostat=io)
      if (io/=0) call abortp('(references_write): opening file '//trim(fname)//' failed')
      write(iu,'(i5,a)') ncenter, " nuclei:"
      do i=1,getNNuc()
         write(iu,'(i4,1x,a2,1x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
      enddo

      write(iu,'(i5,1x,a,1x,i9,i6)') this%getSize(),trim(str),this%getCount(),this%getElemLength()

      do n=1,this%rll%size()
         p => this%rll%elem(idx(n))
         do m=1,p%size()
            rp => p%elem(m)
            call findNucElecs(thresh, rp%x, rp%y, rp%z, naNuc, nbNuc, aNucElec, bNucElec)
            write(iu,'(a,2i5,a,f15.6,a,i6,4i4)') trim(str)//":",n,m," F("//trim(str)//"):", &
               rp%f, "  found:", rp%count, naNuc, getNAlpha(), nbNuc, getNelec()
            write(iu,'(i5)') size(rp%x)
            do i=1,size(rp%x)
               write(iu,'(3f13.6)') rp%x(i),rp%y(i),rp%z(i)
            end do
         end do
      end do
      close(iu)
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

   end subroutine refc_writeFullList

end module refBase_m
