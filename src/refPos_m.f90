! Copyright (C) 2013-2014 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refPos_m

   ! "pos"-type references: identify w.r.t position only: with AND without spin

   use kinds_m, only: r8, i8
   use refBase_m
   use refUtils_m, only: calcRefDifference, calcRefDiffPerm
   use sorting_m, only: sort


   implicit none

   integer, parameter :: SIMPLE=0, SPINPERM=1, ALLPERM=2     ! modes for determining reference differences

   public
   private :: SIMPLE,SPINPERM,ALLPERM

   type, extends(referenceContainer) :: refContainerPos      ! like Struct but electrons are permuted for max likeness and averaged
      real(r8)            :: simThreshold = 0.1d0            ! identify structures as similar if max dist is smaller
      real(r8)            :: simmaxThreshold = 0.1d0         ! for "max positions": similar elec position if distance smaller
      real(r8), allocatable  :: maxrefpos(:,:,:)               ! collect "max positions" of alpha/beta electrons for each reference
      integer, allocatable :: nref(:)                        ! how many actual max positions in maxrefpos
      integer, pointer     :: irep(:,:) => null()            ! pointer to allocated array: list of reference electrons to ignore
      integer, allocatable          :: excl(:)               ! list of "excluded positions"
      integer                       :: exclMode = 0
      type(matrixstat),allocatable  :: pairpositions(:)      ! statistics for alpha/beta position of electron in e pair, for each reference
      type(vectorstat),allocatable  :: maxpositions(:,:)     ! for each reference for each elec pos in reference
      integer           :: totalFoundCount=0
      integer           :: notFoundCount=0
      real(r8)            :: totalWeight=0.d0
   contains
      procedure :: refcp_create
      procedure :: destroy => refcp_destroy
      procedure :: insert => refcp_insert
      procedure :: writeShortList => refcp_writeShortList
      procedure :: writeSimpleList => refcp_writeSimpleList
   end type refContainerPos


contains


   subroutine refcp_create(this,listLength,doSortFreq,refFile,simThreshold,simmaxThreshold,ire,exclFile)
      class(refContainerPos), intent(inout)        :: this
      integer, intent(in)                          :: listLength
      logical, intent(in)                          :: doSortFreq
      character(len=40), intent(in)                :: refFile
      real(r8), optional, intent(in)                 :: simThreshold
      real(r8), optional, intent(in)                 :: simmaxThreshold
      integer, optional, pointer, intent(in)       :: ire(:,:)
      character(len=40), optional, intent(in)      :: exclFile
      real(r8) x(ne),y(ne),z(ne),ri(3),d,f
      integer l,i,io,j,k,n,ref,noel,refNo,refNo2,nomax,nWords
      character(len=120) line
      character(len=20) words(12)
      character(len=1) emode
      type(reference) :: r
      logical found
      integer, parameter :: iu=12
      real(r8) nucThresh,coreThresh,bondDiff

      call assert(.not. (present(ire) .and. present(exclFile)),"refcp_create: ire and exclFile exclude each other")
      call refc_create(this,listLength,doSortFreq)
      if (present(simThreshold)) this%simThreshold = simThreshold
      if (present(simmaxThreshold)) this%simmaxThreshold = simmaxThreshold
      if (present(ire)) this%irep => ire
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

      allocate(this%pairpositions(listLength),this%maxpositions(10*ne,listLength))
      do l=1,listLength
         call this%pairpositions(l)%create(ne,3)
      end do
      do l=1,listLength
         do i=1,10*ne
            call this%maxpositions(i,l)%create(3)
         end do
      end do
      allocate(this%maxrefpos(3,10*ne,listLength),this%nref(listlength))
      this%nref = 0
      call r%new(ne)
      !! read ref file: this is (should be) called only for MASTER
      open(iu,file=refFile,status='old',iostat=io)
      call assert(io==0,' ref file '//trim(refFile)//' does not exist')
      read(iu,*) n   ! ncenter
      do i=1,n
         read(iu,*) j  ! index,elem,x,y,z
      end do

      ref=1
      read(iu,*) nomax
      do
         read(iu,'(A)') line
         call tokenize(line,words,n)
         read(words(2),*) refNo
         read(words(3),*) refNo2
         read(words(5),*) f
         read(iu,*) noel
         call assert(noel==ne,"refcp_create: wrong number of electrons in ref file")
         do i=1,noel
            read(iu,*) x(i),y(i),z(i)
         enddo
         if (refNo == ref .and. refNo2 == 1) then      ! read first of each list
            call r%set(x,y,z,f)
            call this%append(r)
            do i=1,ne
               ri = (/ x(i),y(i),z(i) /)
               found = .false.
               do k=1,this%nref(ref)
                  d = sqrt(dot_product(ri(:)-this%maxrefpos(:,k,ref),ri(:)-this%maxrefpos(:,k,ref)))
                  if (d < this%simmaxThreshold) then
                     found = .true.
                     exit
                  end if
               end do
               if (.not.found) then
                  this%nref(ref) = this%nref(ref) + 1
                  this%maxrefpos(:,this%nref(ref),ref) = ri(:)
               end if
            end do
            ref = ref + 1
            if (ref > listLength) exit
         endif
      enddo
      close(iu)
      call r%destroy()

   end subroutine refcp_create


   subroutine refcp_destroy(this)
      class(refContainerPos), intent(inout) :: this
      integer l,i
      do l=1,size(this%pairpositions)
         call this%pairpositions(l)%destroy()
      end do
      do l=1,size(this%maxpositions,2)
         do i=1,2*ne
            call this%maxpositions(i,l)%destroy()
         end do
      end do
      deallocate(this%maxrefpos,this%pairpositions,this%maxpositions)
      if (associated(this%irep)) deallocate(this%irep)
      call refc_destroy(this)
   end subroutine refcp_destroy


   subroutine refcp_insert(this,r)
      class(refContainerPos), intent(inout) :: this
      type(reference), intent(in)              :: r
      type(reference)                          :: r1
      type(reference), pointer                 :: rp
      type(refl_vlist)                         :: rl
      type(refl_vlist), pointer                :: p
      real(r8)                                   :: x(ne,3),xx(3),rr(3)
      real(r8)                                   :: f,wgt,dist,d,maxdist
      logical found, sfound
      integer i,j,k,ns
      integer, pointer :: irev(:)
      real(r8), parameter :: EPS=1.d-4

      call r%setCount(1)
      call rl%create(this%getElemLength()+1)
      call r1%new(ne)
      !!write(999,*) '*** s insert value=',r%f
      sfound = .false.
      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         r1 = r
         ns = 0
         dist = 0
         if (allocated(this%excl)) then
            call calcRefDifference(rp,r1,ALLPERM,dist,maxdist,exclist=this%excl,exclmode=this%exclMode)
         else if (associated(this%irep)) then
            irev => this%irep(:,i)
            call calcRefDifference(rp,r1,ALLPERM,dist,maxdist,irev=irev)
         else
            call calcRefDifference(rp,r1,ALLPERM,dist,maxdist)
         endif
         if (maxdist < this%simThreshold) sfound = .true.
         !!write(999,*) '   ** ALLPERM dist,ns=',dist,this%distanceThreshold,ns,this%nSameMin
         if (sfound) then
            this%totalFoundCount = this%totalFoundCount + 1
            !!write(999,*) '   ** allperm struct found:',i
            r1 = r
            call calcRefDiffPerm(rp,r1,SPINPERM)
            call r1%get(x(:,1),x(:,2),x(:,3),f)
            wgt = exp(-r%f+rp%f)
            !!write(999,*) '    * adding with weight ',wgt,r%f,rp%f
            !!do j=1,size(r1%x)
            !!   write(999,'(3f13.6)') x(j,1),x(j,2),x(j,3)
            !!end do
            if (wgt > EPS) then
               this%totalWeight = this%totalWeight + wgt
               call this%pairpositions(i)%add(x,wgt)
               do j=1,ne
                  found = .false.
                  rr(1) = r%x(j); rr(2) = r%y(j); rr(3) = r%z(j)
                  do k=1,this%nref(i)
                     xx = this%maxrefpos(:,k,i)
                     d = sqrt(dot_product((xx-rr),(xx-rr)))
                     if (d < this%simmaxThreshold) then
                        found = .true.
                        !!write(999,*) ' ** found ',k,this%nref(i),d,this%posThreshold
                        call this%maxpositions(k,i)%add(rr,wgt)
                        exit
                     end if
                  end do
                  if (.not.found) then
                     if (this%nref(i)<10*ne) then
                        this%nref(i) = this%nref(i) + 1
                        !!write(999,*) ' ** not found ',this%nref(i),this%posThreshold
                        this%maxrefpos(:,this%nref(i),i) = rr
                        call this%maxpositions(this%nref(i),i)%add(rr,wgt)
                     else
                        !!write(999,*) ' ** not found: ref pos overflow at ',i
                     end if
                  end if
               end do
            end if
            exit
         end if
      end do
      if (.not.sfound) this%notFoundCount = this%notFoundCount + 1
      call rl%destroy()
      call r1%destroy()
      !!write(999,*) '*** done'
   end subroutine refcp_insert


   subroutine refcp_writeShortList(this,iu,str)
      class(refContainerPos)                :: this
      integer, intent(in)                   :: iu     ! unit to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      integer i,j
      integer(i8) c
      real(r8) pos(3),sd(3),var(3),w
      integer, parameter :: MAXSHOW = 15
      integer, parameter :: iuf = 19

      call assert(MASTER,'refcp_writeShortList: shall be called only for MASTER')

      write(iu,'(a)') ' all maximum positions with similarity to reference configurations'
      write(iu,'(a,i10)')  ' total sample size  : ',this%totalFoundCount
      write(iu,'(a,f6.3)') ' not matching ratio : ',dble(this%notFoundCount)/(this%totalFoundCount+this%notFoundCount)
      write(iu,'(2a)') '  ref   count weightratio            pos                                ', &
                       'sigma                              stderr'
      write(iu,'(2a)') '------------------------------------------------------------------------', &
                       '-----------------------------------------'
      open(iuf,file=trim(basename)//'.max')
      write(iuf,'(a)') 'all max positions with similarity to reference configurations'
      write(iuf,'(a,i10,a,f6.3)') 'total sample size:',this%totalFoundCount,' no-match-ratio:',  &
         dble(this%notFoundCount)/(this%totalFoundCount+this%notFoundCount)
      write(iuf,*) this%getSize()
      do i=1,this%getSize()
         write(iu,'(a,i4)') ' ref=',i
         write(iuf,'(a,i4,a,i4)') ' ref=',i,' #pos=',this%nref(i)
         do j=1,this%nref(i)
            c   = this%maxpositions(j,i)%localcount()
            if (c > 0) then
               pos = this%maxpositions(j,i)%localmean()
               w   = this%maxpositions(j,i)%localweight()
               var = this%maxpositions(j,i)%localvar()
               sd  = this%maxpositions(j,i)%localstddev()
            else
               pos = 0
               w = 0
               var = 0
               sd  = 0
            end if
            write(iu,'(i4,i12,f8.3,3f11.4,2(3x,3f11.4))') j,c,w/this%totalWeight,pos(:),sqrt(var(:)),sd(:)
            write(iuf,'(i4,i12,f8.3,3f11.4,2(3x,3f11.4))') j,c,w/this%totalWeight,pos(:),sqrt(var(:))
         end do
      end do
      close(iuf)
   end subroutine refcp_writeShortList

   subroutine refcp_writeSimpleList(this,fname,str)
      class(refContainerPos)                :: this
      character(len=*), intent(in)          :: fname  ! file to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      integer n,i,io
      integer(i8) cnt
      integer, parameter            :: iu = 13
      real(r8) thresh
      real(r8) :: mean(ne,3), stddev(ne,3)

      !! careful: usually only from from MASTER

      thresh = 0.002d0     ! could be input parameter

      open(iu,file=trim(fname)//'.ref',iostat=io)
      if (io/=0) call abortp('(references_write): opening file '//trim(fname)//' failed')
      write(iu,'(i5,a)') ncenter, " nuclei:"
      do i=1,getNNuc()
         write(iu,'(i4,1x,a2,1x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
      enddo

      write(iu,'(i5,1x,a)') this%getSize(), "PREF"

      do n=1,this%getSize()
         !!write(999,*) " *** writeSimpleList: ",n,mytid
         cnt = this%pairpositions(n)%localcount()
         !!write(999,*) "   * cnt: ",n,mytid,cnt
         if (cnt > 0) then
            mean = this%pairpositions(n)%localmean()
            !! write(999,'(a,2i4,100g10.3)') "   * mean: ",mytid,n,mean(:,1),mean(:,2),mean(:,3)
            stddev = this%pairpositions(n)%localstddev()
            !! write(999,'(a,2i4,100g10.3)') "   * stddev: ",mytid,n,stddev(:,1),stddev(:,2),stddev(:,3)
         else
            mean = 0
            stddev = 0
         end if
         if (MASTER) then
            write(iu,'(a,i5,a,i10,a)') "PAIR:",n,"  1 F(ref): 0.0 found: ",cnt,"  0 0 0 0"
            write(iu,'(i5)') ne
            do i=1,ne
               write(iu,'(3f13.6,3x,3f13.6)') mean(i,:),stddev(i,:)
            end do
         end if
      end do
      if (MASTER) close(iu)
   end subroutine refcp_writeSimpleList

end module refPos_m
