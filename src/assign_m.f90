! Copyright (C) 2013 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module assign_m

   use kinds_m, only: r8
   use global_m, only: ne,iul
   use wfData_m, only: ncenter,nalpha,nbeta,atoms,calcNucElecDists
   use error_m
   use parsing_m
   use utils_m, only: tokenize
   use randomWalker_m
   use hungarian_m, only: munkres
   implicit none

   private
   public :: assign_init, assign_destroy, assign_readRef, assign_getRef, assign_setRef, assign_get

   integer, parameter   :: SQUARED=1, SIMPLE=2

   integer              :: mCoreEl = 0      ! # of alpha electrons (= #beta elecs) assigned greedy as core
   integer, allocatable :: mCoreList(:)     ! list of alpha core core electrons at nuclei
   integer              :: mCostMatrix = SIMPLE   ! cost matrix type
   integer              :: mRef = 0         ! reference
   integer              :: mRef2 = 0        ! 2nd reference
   real(r8),allocatable   :: mX(:), mY(:), mZ(:)
   real(r8),allocatable   :: mXRef(:), mYRef(:), mZRef(:)
   real(r8),allocatable     :: maData(:,:), mbData(:,:)

contains

   !-------------------------------
   subroutine assign_init(lines,nl)
   !-------------------------------
      integer, intent(in)           :: nl
      character(len=120), intent(in) :: lines(:)
      integer                       :: iflag
      integer                       :: j,a
      character(len=7)              :: s

      mCoreEl = 0
      allocate(mCoreList(ncenter))
      mCoreList = 0
      j = 0
      if (.not.finda(lines,nl,'nocoresep')) then
         do a=1,ncenter
            if (atoms(a)%za > 2) then
               j = j + 1
               mCoreList(j) = a
            endif
         enddo
         mCoreEl = j
      endif

      call getstra(lines,nl,'assign=',s,iflag)
      if (iflag == 0) then
         if (s=='squared') then
            mCostMatrix = SQUARED
         else if (s=='simple') then
            mCostMatrix = SIMPLE
         else
            call abortp("illegal value given for assign=")
         end if
      else
         mCostMatrix = SIMPLE
      end if

      ! index of reference to read from file
      mRef = 1
      call getinta(lines,nl,'ref_nr=',mRef,iflag)
      mRef2 = 1
      call getinta(lines,nl,'ref_nr2=',mRef2,iflag)

      ! allocate arrays
      allocate(mX(ne),mY(ne),mZ(ne),mXRef(ne),mYRef(ne),mZRef(ne))
      allocate(maData(nalpha-mCoreEl,nalpha-mCoreEl),mbData(nbeta-mCoreEl,nbeta-mCoreEl))
   end subroutine assign_init


   !-------------------------!
   subroutine assign_destroy()
   !-------------------------!
      deallocate(mX,mY,mZ,mXRef,mYRef,mZRef)
      deallocate(maData,mbData)
      deallocate(mCoreList)
   end subroutine assign_destroy


   !-----------------------------------------------!
   subroutine assign_readRef(lines,nl,hasbeenSorted)
   !-----------------------------------------------!
      ! read reference from ref-File
      integer, intent(in)            :: nl
      character(len=120), intent(in) :: lines(:)
      logical, intent(out)           :: hasbeenSorted

      integer, parameter            :: iu=45
      integer                       :: i,j,iflag,io,n
      real(r8)                        :: x_tmp2(ne),y_tmp2(ne),z_tmp2(ne)
      real(r8)                        :: F,Ftmp
      integer                       :: nomax,noel
      character(len=120)            :: line
      character(len=20)             :: words(12)
      integer                       :: refNo,refNo2,foundtmp,found
      character(len=40)             :: fname

      call getstra(lines,nl,'ref_file=',fname,iflag)
      if (iflag /= 0) then
         fname = trim(baseName)//'.ref'
      end if
      open(iu,file=fname,status='old',iostat=io)
      call assert(io==0,' ref file '//trim(fname)//' does not exist')
      read(iu,*) n   ! ncenter
      do i=1,n
         read(iu,*) j  ! index,elem,x,y,z
      end do

      read(iu,*) nomax
      do j=1,nomax
         read(iu,'(A)') line
         call tokenize(line,words,n)
         read(words(2),*) refNo
         read(words(3),*) refNo2
         read(words(5),*) Ftmp
         read(words(7),*) foundtmp
         read(iu,*) noel
         do i=1,noel
            read(iu,*) x_tmp2(i),y_tmp2(i),z_tmp2(i)
         enddo
         if (refNo == mRef .and. refNo2 == mRef2) then
            mXRef = x_tmp2
            mYRef = y_tmp2
            mZRef = z_tmp2
            found = foundtmp
            F = Ftmp
            exit
         endif
      enddo

      if (MASTER .and. logmode >=2) then
         write(iul,'(3A)') "reference file ",trim(fname)," read!"
         write(iul,'(A,2I3,A,F13.5,A,I8)') "reference indices: ",mRef,mRef2," with function value: ",F," found:",found
         do i=1,ne
             write(iul,'(3F12.5)') mXRef(i),mYRef(i),mZRef(i)
         enddo
         write(iul,*)
      endif

      close(iu)

      if (mCoreEl>0) call findCoreElecsSwap(mXRef,mYRef,mZRef,hasbeenSorted)

   end subroutine assign_readRef

   !-----------------------------!
   subroutine assign_getRef(x,y,z)
   !-----------------------------!
      real(r8),intent(out) :: x(ne),y(ne),z(ne)
      x = mXRef
      y = mYRef
      z = mZRef
   end subroutine assign_getRef

   !-------------------------------------------!
   subroutine assign_setRef(x,y,z,hasbeenSorted)
   !-------------------------------------------!
      real(r8),intent(in)  :: x(:),y(:),z(:)
      logical, intent(out) :: hasbeenSorted

      mXRef = x
      mYRef = y
      mZRef = z

      if (mCoreEl>0) call findCoreElecsSwap(mXRef,mYRef,mZRef,hasbeenSorted)
   end subroutine assign_setRef

   !----------------------------------!
   subroutine assign_get(rwp,asgn,dist)
   !----------------------------------!
     
      ! this routine identifies the "assignment" using the assignment algorithm (Hungarian algorithm, munkres)
      ! meaning that asgn(:,w) contains the permutation of the w-th walker in walker block rwp
      ! The "permutation" is represented by a list of indices:
      ! asgn(1,w) is the index of the 1st electron after permutation
      ! asgn(i,w)=1 means: the 1st electron after permutation is electron i
      ! note that asgn does not permute the electrons. This can be done with asgn if necessary
      ! if mCoreEl is set in assign_init, the core electrons are assigned "greedy", then all
      ! other electrons are assigned with the Hungarian algorithm using the definition of
      ! the cost matrix as given in assign_init.
      ! On algorithm:
      ! The assignment is done separately for alpha and beta electrons.
      ! findCoreElecsIdx constructs the permutation corresponding to the greedy core assignment: "idxc"
      ! the munkres code works on the core-permuted electrons using idxc
      ! The final permutation index asgn contains the final permutation wrt to the original electrons
      ! this is therefore the concatenation of the munkres permutation and the core permutation
      !
      type(RandomWalker),pointer,intent(in) :: rwp(:)
      integer,intent(inout)                 :: asgn(:,:)
      real(r8), optional                      :: dist(:)
      real, allocatable                     :: dmata(:,:),dmatb(:,:)
      real(r8) dista,distb
      integer i,j,beta1,w
      integer idxa(nalpha-mCoreEl),idxb(nbeta-mCoreEl),idxc(ne)

      call assert(associated(rwp),"assign_get:: rw pointer not associated")
      call assert(allocated(mX).and.allocated(mXRef),"assign_get:: not initialized")

      do w = 1,size(rwp)
         call pos(rwp(w),mX,mY,mZ)

         idxc = (/ (i,i=1,ne) /)      

         if (mCoreEl>0) call findCoreElecsIdx(mX,mY,mZ,idxc)

         select case (mCostMatrix)
         case (SIMPLE)
            forall (i=mCoreEl+1:nalpha,j=mCoreEl+1:nalpha)
               maData(i-mCoreEl,j-mCoreEl) = sqrt( (mXRef(i)-mX(idxc(j)))**2 + & 
                                                   (mYRef(i)-mY(idxc(j)))**2 + (mZRef(i)-mZ(idxc(j)))**2 )
            end forall
            beta1 = nalpha+mCoreEl
            forall (i=beta1+1:ne,j=beta1+1:ne)
                mbData(i-beta1,j-beta1) = sqrt( (mXRef(i)-mX(idxc(j)))**2 + &
                                                (mYRef(i)-mY(idxc(j)))**2 + (mZRef(i)-mZ(idxc(j)))**2 )
            end forall
         case (SQUARED)
            forall (i=mCoreEl+1:nalpha,j=mCoreEl+1:nalpha)
               maData(i-mCoreEl,j-mCoreEl) = (mXRef(i)-mX(idxc(j)))**2 + & 
                                             (mYRef(i)-mY(idxc(j)))**2 + (mZRef(i)-mZ(idxc(j)))**2 
            end forall
            beta1 = nalpha+mCoreEl
            forall (i=beta1+1:ne,j=beta1+1:ne)
               mbData(i-beta1,j-beta1) = (mXRef(i)-mX(idxc(j)))**2 + &
                                         (mYRef(i)-mY(idxc(j)))**2 + (mZRef(i)-mZ(idxc(j)))**2
            end forall
         case default
            call abortp("assign_get: illegal mCostMatrix")
         end select
         call munkres(1, maData, nalpha-mCoreEl, nalpha-mCoreEl, idxa, dista)
         call munkres(1, mbData, nbeta-mCoreEl, nbeta-mCoreEl, idxb, distb)

         asgn(1:mCoreEl,w)               = idxc(1:mCoreEl)
         asgn(mCoreEl+1:nalpha,w)        = idxc(idxa(:) + mCoreEl)
         asgn(nalpha+1:nalpha+mCoreEl,w) = idxc(nalpha+1:nalpha+mCoreEl)
         asgn(nalpha+mCoreEl+1:ne,w)     = idxc(idxb(:) + mCoreEl+nalpha)

         if (present(dist)) then
            dist(w) = dista + distb
         end if
      end do

   end subroutine assign_get

   subroutine findCoreElecsSwap(x,y,z,hasbeenSwapped)
      ! identifies K shell core electrons by comparing distances
      ! swap electrons such that the core electrons are the first in alpha and
      ! beta list
      real(r8), intent(inout) :: x(:),y(:),z(:)  ! coords
      logical, intent(inout) :: hasbeenSwapped
      integer ia,a,iamin,ibmin,ica,icb
      real(r8) tmp
      real(r8) rai(ncenter,ne)

      hasbeenSwapped = .false.
      call calcNucElecDists(x,y,z,rai)

      do ia=1,mCoreEl
         a = mCoreList(ia)
         iamin = minloc(rai(a,1:nalpha),dim=1)
         ibmin = minloc(rai(a,nalpha+1:ne),dim=1) + nalpha
         ica = ia
         icb = nalpha + ia
         ! swap electrons
         if (ica /= iamin) then
            tmp = x(iamin); x(iamin) = x(ica); x(ica) = tmp
            tmp = y(iamin); y(iamin) = y(ica); y(ica) = tmp
            tmp = z(iamin); z(iamin) = z(ica); z(ica) = tmp
            hasbeenSwapped = .true.
         end if
         if (icb /= ibmin) then
            tmp = x(ibmin); x(ibmin) = x(icb); x(icb) = tmp
            tmp = y(ibmin); y(ibmin) = y(icb); y(icb) = tmp
            tmp = z(ibmin); z(ibmin) = z(icb); z(icb) = tmp
            hasbeenSwapped = .true.
         end if
         call calcNucElecDists(x,y,z,rai)
      end do
   end subroutine findCoreElecsSwap

   subroutine findCoreElecsIdx(x,y,z,idx)
      ! identifies K shell core electrons by comparing distances
      ! idx is permutation as idx list, isChanged signals deviation from no permuation
      ! permutation swaps core electrons in the order of the nuclei to the
      ! first positions in alpha and beta part
      ! idx is defined such that idx(i) is the real index of the i-th electron after permutation
      ! the alpha electron closest to the 1st nucleus is idx(1)
      !   - " -                           2nd            idx(2) ...
      ! the beta electron closest to the 1st nucleus has idx(nalpha+1) ...

      real(r8), intent(inout) :: x(:),y(:),z(:)  ! coords
      integer, intent(inout) :: idx(:)         ! MUST be on input 1,2,3, ..,n 
      integer ia,a,iamin,ibmin,ica,icb,i,j,itmp,iamin2,ibmin2
      integer cidx(nalpha)
      real(r8) rai(ncenter,ne),rmin2
      logical changed

      changed = .false.
      call calcNucElecDists(x,y,z,rai)

      do ia=1,mCoreEl
         a = mCoreList(ia)
         iamin = minloc(rai(a,1:nalpha),dim=1)
         if (idx(iamin)==0) then
            write(iul,*) "Warning: attempt to doubly assign core electron"
            ! find 2nd best
            iamin2 = 0
            rmin2  = 1000.d0
            do i=1,nalpha
               if (rai(a,i) < rmin2 .and. i/=iamin) then
                  iamin2 = i
                  rmin2 = rai(a,i)
               end if
            end do
            iamin = iamin2
         end if
         cidx(ia) = iamin
         idx(iamin) = 0
      end do
      j=0
      do i=mCoreEl+1,nalpha
         do
            j = j+1
            if (idx(j) /= 0) then
               cidx(i) = idx(j)
               exit
            end if
         end do
      end do
      idx(1:nalpha) = cidx(1:nalpha)

      do ia=1,mCoreEl
         a = mCoreList(ia)
         ibmin = minloc(rai(a,nalpha+1:ne),dim=1) + nalpha
         if (idx(ibmin)==0) then
            write(iul,*) "Warning: attempt to doubly assign core electron"
            ! find 2nd best
            ibmin2 = 0
            rmin2  = 1000.d0
            do i=nalpha+1,ne
               if (rai(a,i) < rmin2 .and. i/=ibmin) then
                  ibmin2 = i
                  rmin2 = rai(a,i)
               end if
            end do
            ibmin = ibmin2
         end if
         cidx(ia) = ibmin
         idx(ibmin) = 0
      end do
      j=nalpha
      do i=mCoreEl+1,ne-nalpha
         do
            j = j+1
            if (idx(j) /= 0) then
               cidx(i) = idx(j)
               exit
            end if
         end do
      end do
      idx(nalpha+1:ne) = cidx(1:ne-nalpha)
   end subroutine findCoreElecsIdx

end module assign_m
