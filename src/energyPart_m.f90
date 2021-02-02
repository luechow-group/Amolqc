! Copyright (C) 2013, 2015, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module energyPart_m

   ! energy partitioning based on local energy partitioning into
   ! T(i),Vee(i,j),Ven(i,j)
   ! simple block statistic with add: walker and electron permutation
   ! walker must contain the current T(i),Vee(i,j),Ven(i,j) contributions
   ! This version allows energy partitioning for several "domains"
   ! Total energy (including T,Vee,Ven) is exactly the sum over all contributions

   use kinds_m, only: r8, i8
   use global_m
   use randomWalker_m
   use newStatistics_m
   use assign_m, only: assign_getRef
   use wfData_m, only: Vnni,atoms,do_epart
   use utils_m, only: tokenize
   implicit none

   private
   public :: epart_init, epart_destroy, epart_reset, epart_add, epart_write

   type :: chemGroup
      integer, pointer :: elecs(:) => null()
      integer, pointer :: nucs(:) => null()
      integer, pointer :: attached(:) => null()
      character(len=4) :: typ = ''
   contains
      procedure :: init => chemGroup_init
      procedure :: destroy => chemGroup_destroy
      procedure :: write => chemGroup_write
      procedure :: intraEnergy => chemGroup_intraEnergy
      procedure :: interEnergy => chemGroup_interEnergy
   end type

   type(blockvectorstat), allocatable :: mEkin(:)
   type(blockmatrixstat), allocatable :: mVne(:)
   type(blockmatrixstat), allocatable :: mVee(:)
   type(blockvectorstat), allocatable :: mCGintra(:)
   type(blockmatrixstat), allocatable :: mCGinter(:)
   type(chemGroup), pointer :: mCGlist(:) => null()

   integer :: mKMax = 0           ! # of domains
   logical :: mDomain = .false.   ! domain or sed based partitioning

contains

   subroutine epart_init(lines,nl)
   !------------------------------
      integer, intent(in)            :: nl
      character(len=120), intent(in) :: lines(:)
      integer, allocatable           :: elecidx(:),nucidx(:),atidx(:)
      character(len=10) :: token(50)
      character(len=4)  :: typ
      integer ncg,blockSize, iflag, i,j,k,nToken, idx,nelec,nnuc,nat, kmax

      call assert(do_epart,"epart_init: epart not set in wf file")
      mDomain = .not.finda(lines,nl,"$sed")

      !!!blockSize = 100
      call getinta(lines,nl,'blocklen=',blockSize,iflag)
      if (iflag /= 0) call abortp("epart_init: blocklen required")
      kmax = 1
      call getinta(lines,nl,'kmax=',kmax,iflag)
      mKMax = kmax
      allocate(mEkin(kmax),mVne(kmax),mVee(kmax),mCGintra(kmax),mCGinter(kmax))
      do i=1,kmax
         call mEkin(i)%create(ne,blockSize)
         call mVne(i)%create(ncenter,ne,blockSize)
         call mVee(i)%create(ne,ne,blockSize)
      enddo
      !!! note: currently chem groups only for first domain n==1 
      ncg = 0
      call getinta(lines,nl,'chem_groups=',ncg,iflag)
      if (iflag==0 .and. ncg>0) then
         do i=1,nl
            if (index(lines(i),'chem_groups')>0) exit
         enddo
         call assert(nl-i >= ncg,'epart_init:: not enough chem_group lines')
         allocate(mCGlist(ncg))
         do j=i+1,i+ncg
            call tokenize(lines(j),token,nToken)
            read(token(1),*) idx
            call assert(idx==j-i,'epart_init:: chem_group list has incorrect format(1)')
            read(token(2),'(a)') typ
            read(token(3),*) nelec
            read(token(4),*) nnuc
            read(token(5),*) nat
            call assert(nToken==nelec+nnuc+nat+5,'epart_init:: chem_group list has incorrect format(2)')
            allocate(elecidx(nelec),nucidx(nnuc),atidx(nat))
            do k=1,nelec
               read(token(5+k),*) elecidx(k)
            enddo
            do k=1,nnuc
               read(token(5+nelec+k),*) nucidx(k)
            enddo
            do k=1,nat
               read(token(5+nelec+nnuc+k),*) atidx(k)
            enddo
            call mCGlist(idx)%init(elecidx,nucidx,atidx,typ)
            deallocate(elecidx,nucidx,atidx)
         enddo
         call mCGintra(1)%create(ncg,blockSize)
         call mCGinter(1)%create(ncg,ncg,blockSize)
      endif
   end subroutine epart_init


   subroutine epart_destroy()
   !------------------------!
      integer i
      do i=1,mKMax
         call mEkin(i)%destroy()
         call mVne(i)%destroy()
         call mVee(i)%destroy()
         if (mCGintra(i)%isCreated()) call mCGintra(i)%destroy()
         if (mCGinter(i)%isCreated()) call mCGinter(i)%destroy()
      enddo
      if (associated(mCGlist)) then
         do i=1,size(mCGlist)
            call mCGlist(i)%destroy()
         enddo
         deallocate(mCGlist)
         mCGlist => null()
      endif
      deallocate(mEkin,mVne,mVee,mCGintra,mCGinter)
   end subroutine



   subroutine epart_reset()
   !----------------------!
      integer i
      do i=1,mKMax
         call mEkin(i)%reset()
         call mVne(i)%reset()
         call mVee(i)%reset()
      enddo
      if (associated(mCGlist)) then
         call mCGintra(1)%reset()
         call mCGinter(1)%reset()
      endif
   end subroutine



   subroutine epart_add(rw,asgn,dIdx)
   !----------------------------!
      ! use assignment "asgn" to collect energy contributions
      type(RandomWalker)            :: rw
      integer, intent(in)           :: asgn(:)   ! assign permutation for each walker
      integer, intent(in), optional :: dIdx     ! domain idx    
      real(r8)            :: Ekin(ne),Vee(ne,ne),Vne(ncenter,ne)
      real(r8)            :: ee_tmp(ne,ne),ne_tmp(ncenter,ne),kin_tmp(ne)
      real(r8)            :: weight
      real(r8), allocatable :: Eintra(:),Einter(:,:)
      integer           :: i,j,ii,jj,a,n

      if (present(dIdx)) then
         n = dIdx
      else
         n = 1
      endif

      if (associated(mCGlist)) then
         allocate(Eintra(size(mCGlist)),Einter(size(mCGlist),size(mCGlist)))
         Eintra = 0; Einter = 0
      endif

      call Ekini(rw,kin_tmp)
      call EVee(rw,ee_tmp)
      call EVne(rw,ne_tmp)
      weight = wgt(rw)

      ! filling up array required for assignment
      forall (i=1:ne)
         forall (j=i+1:ne)
            ee_tmp(j,i) = ee_tmp(i,j)
         end forall
      end forall

      forall (i=1:ne)
         Ekin(i) = kin_tmp(asgn(i))
         forall (j=i+1:ne)
            Vee(j,i) = ee_tmp(asgn(i),asgn(j))
            Vee(i,j) = Vee(j,i)
         end forall
         forall (a=1:ncenter)
            Vne(a,i) = ne_tmp(a,asgn(i))
         end forall
      end forall

      call mEkin(n)%add(Ekin,weight)
      call mVne(n)%add(Vne,weight)
      call mVee(n)%add(Vee,weight)

      if (associated(mCGlist) .and. n==1) then
         do i=1,size(mCGlist)
            Eintra(i) = mCGlist(i)%intraEnergy(Ekin,Vne,Vee)
         end do
         Einter = 0
         do i=1,size(mCGlist)
            do j=i+1,size(mCGlist)
               Einter(i,j) = mCGlist(i)%interEnergy(mCGlist(j),Ekin,Vne,Vee)
               Einter(j,i) = Einter(i,j)
            end do
         end do
         do i=1,size(mCGlist)
            do j=1,size(mCGlist(i)%attached)
               jj = mCGlist(i)%attached(j)
               Eintra(i) = Eintra(i) + Einter(i,jj)
               Einter(i,jj) = 0
               Einter(jj,i) = 0
            end do
         end do
         call mCGintra(1)%add(Eintra,weight)
         call mCGinter(1)%add(Einter,weight)
      end if
      if (associated(mCGlist)) then
         deallocate(Einter,Eintra)
      end if

   end subroutine epart_add


   subroutine epart_write(iu,mode)
   !-----------------------------!

      integer, intent(in) :: iu    ! file unit to write to
      integer, intent(in) :: mode  ! 1: for log file, 2: for data file
      real(r8)  :: etot
      real(r8)  :: Ekin(ne,mKMax),Vee(ne,ne,mKMax),Vne(ncenter,ne,mKMax)
      real(r8)  :: Ekinsd(ne,mKMax),Veesd(ne,ne,mKMax),Vnesd(ncenter,ne,mKMax)
      real(r8)  :: x(ne),y(ne),z(ne)
      real(r8), allocatable :: Eintra(:), Einter(:,:), Eintrasd(:), Eintersd(:,:)
      integer i,j,n
      integer(i8) nd

      do i=1,mKMax
         nd = mEkin(i)%count()
         if (nd > 1) then
            Ekin(:,i)     = mEkin(i)%mean()
            Vee(:,:,i)    = mVee(i)%mean()
            Vne(:,:,i)    = mVne(i)%mean()
            Ekinsd(:,i)   = mEkin(i)%stddev()
            Veesd(:,:,i)  = mVee(i)%stddev()
            Vnesd(:,:,i)  = mVne(i)%stddev()
         else
            Ekin(:,i)     = 0.d0
            Vee(:,:,i)    = 0.d0
            Vne(:,:,i)    = 0.d0
            Ekinsd(:,i)   = 0.d0
            Veesd(:,:,i)  = 0.d0
            Vnesd(:,:,i)  = 0.d0
         endif
      enddo

      if (MASTER) then
         if (mode==1) then
            write(iu,*)
            write(iu,*) "Energy partitioning results:"
            write(iu,*) "(see full results in .see file)"
            write(iu,*)
         else if (mode==2) then
            write(iu,*) ne," ",ncenter
            write(iu,*) "Geometry"
            do i=1,ncenter
                write(iu,'(i4,1x,a2,1x,3f12.5)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
            enddo
            write(iu,*)
            if (.not.mDomain) then
               call assign_getRef(x,y,z)
               write(iu,*) "Reference:"
               do i=1,ne
                   write(iu,'(3f12.5)') x(i),y(i),z(i)
               enddo
            endif
            do n=1,mKMax
               if (mKMax > 1) then
                  write(iu,'(//a,i5)') "Domain # ",n
               endif
               etot = 0
               write(iu,*)
               write(iu,*) "Ekin(i):"
               do i=1,ne
                   write(iu,'(i4,2f15.5)') i, Ekin(i,n),Ekinsd(i,n)
                   etot = etot + Ekin(i,n)
               enddo
               write(iu,*)
               write(iu,*) "Vee(i,j):"
               do i=1,ne
                   do j=i+1,ne
                       write(iu,'(2i4,2f15.5)') i,j,Vee(i,j,n),Veesd(i,j,n)
                       etot = etot + Vee(i,j,n)
                   enddo
               enddo
               write(iu,*)
               write(iu,*) "Ven(i,j)"
               do i=1,ne
                   do j=1,ncenter
                       write(iu,'(2i4,2f15.5)') i,j,Vne(j,i,n),Vnesd(j,i,n)
                       etot = etot + Vne(j,i,n)
                   enddo
               enddo
               write(iu,*)
               write(iu,*) "Vnn(i,j)"
               do i=1,ncenter
                   do j=i+1,ncenter
                       write(iu,'(2i4,f15.5)') i,j,Vnni(i,j)
                       etot = etot + Vnni(i,j)
                   enddo
               enddo
               write(iu,*)
               write(iu,'(a,f15.5)') "sum of energies: Etot=",etot
            enddo
         end if
      end if

      if (associated(mCGlist)) then
         allocate(Eintra(size(mCGlist)),Einter(size(mCGlist),size(mCGlist)))
         allocate(Eintrasd(size(mCGlist)),Eintersd(size(mCGlist),size(mCGlist)))
         Eintra = mCGintra(1)%mean()
         Eintrasd = mCGintra(1)%stddev()
         Einter = mCGinter(1)%mean()
         Eintersd = mCGinter(1)%stddev()
         if (MASTER) then
            write(iu,*)
            write(iu,*) " chemical group definition:"
            do i=1,size(mCGlist)
               write(iu,'(i4,1x,a4,100i3)') i,mCGlist(i)%typ,size(mCGlist(i)%elecs),size(mCGlist(i)%nucs), &
                 size(mCGlist(i)%attached),mCGlist(i)%elecs,mCGlist(i)%nucs,mCGlist(i)%attached
            end do
            write(iu,*)
            write(iu,*) " chemical group energies:"
            etot = 0
            do i=1,size(mCGlist)
               write(iu,'(i4,1x,a4,2f15.5)') i,mCGlist(i)%typ,Eintra(i),Eintrasd(i)
               etot = etot + Eintra(i)
            end do
            write(iu,*)
            write(iu,*) " inter group energies:"
            do i=1,size(mCGlist)
               do j=i+1,size(mCGlist)
                  write(iu,'(2i4,1x,a4,a1,a4,2f15.5)') i,j,mCGlist(i)%typ,'-',mCGlist(j)%typ,Einter(i,j),Eintersd(i,j)
                  etot = etot + Einter(i,j)
               end do
            end do
            write(iu,*)
            write(iu,'(a,f15.5)') "sum of chem group energies: Etot=",etot
         end if
         deallocate(Eintra,Einter,Eintrasd,Eintersd)
      end if
   end subroutine epart_write

   subroutine chemGroup_init(self,elidx,nucidx,atidx,typ)
      class(chemGroup), intent(inout) :: self
      integer, intent(in)             :: elidx(:),nucidx(:),atidx(:)
      character(len=4), intent(in)    :: typ
      call assert(.not.(associated(self%nucs).or.associated(self%elecs)),"chemGroup_init:: already associated")
      allocate(self%nucs(size(nucidx)),self%elecs(size(elidx)),self%attached(size(atidx)))
      self%elecs = elidx
      self%nucs = nucidx
      self%attached = atidx
      self%typ = typ
   end subroutine

   subroutine chemGroup_destroy(self)
      class(chemGroup), intent(inout) :: self
      call assert(associated(self%nucs).and.associated(self%elecs),"chemGroup_destroy:: not associated")
      deallocate(self%nucs,self%elecs,self%attached)
   end subroutine

   subroutine chemGroup_write(self,iu)
      class(chemGroup), intent(in) :: self
      integer, intent(in)          :: iu
      integer i
      write(iu,'(a,100i3)') "chemGroup "//self%typ//" - nucs:",self%nucs
      write(iu,'(a,100i3)') "               - elec:",self%elecs
   end subroutine

   real(r8) function chemGroup_intraEnergy(self,ekin,vne,vee)
      class(chemGroup), intent(in) :: self
      real(r8)                      :: ekin(:),vne(:,:),Vee(:,:)
      real(r8) sum
      integer i,j,a
      sum = 0
      do i=1,size(self%elecs)
         sum = sum + ekin(self%elecs(i))
      end do
      do a=1,size(self%nucs)
         do i=1,size(self%elecs)
            sum = sum + vne(self%nucs(a),self%elecs(i))
         end do
      end do
      do i=1,size(self%elecs)
         do j=i+1,size(self%elecs)
            sum = sum + vee(self%elecs(i),self%elecs(j))
         end do
      end do
      ! note Vnni could be calculated once and stored
      do i=1,size(self%nucs)
         do j=i+1,size(self%nucs)
            sum = sum + Vnni(self%nucs(i),self%nucs(j))
         end do
      end do
      chemGroup_intraEnergy = sum
   end function

   real(r8) function chemGroup_interEnergy(self,other,ekin,vne,vee)
      class(chemGroup), intent(in) :: self
      type(chemGroup), intent(in) :: other
      real(r8)                      :: ekin(:),vne(:,:),Vee(:,:)
      real(r8) sum
      integer i,j,a,b
      sum = 0
      do a=1,size(self%nucs)
         do i=1,size(other%elecs)
            sum = sum + vne(self%nucs(a),other%elecs(i))
         end do
      end do
      do a=1,size(other%nucs)
         do i=1,size(self%elecs)
            sum = sum + vne(other%nucs(a),self%elecs(i))
         end do
      end do
      do i=1,size(other%elecs)
         do j=1,size(self%elecs)
            sum = sum + vee(other%elecs(i),self%elecs(j))
         end do
      end do
      ! note Vnni could be calculated once and stored
      do i=1,size(self%nucs)
         do j=1,size(other%nucs)
            sum = sum + Vnni(self%nucs(i),other%nucs(j))
         end do
      end do
      chemGroup_interEnergy = sum
   end function


end module energyPart_m
