! Copyright (C) 2014-2015, 2017-2018 Arne Luechow
! Copyright (C) 2015-2016 Christoph Schulte
! Copyright (C) 2015 Marko Hermsen
!
! SPDX-License-Identifier: GPL-3.0-or-later


! maxdomain_m.f90 contains maxdomainModule
! which collects the data for the probability maximum domains
! domain is determined with refVector of references (and the hungarian algo)
! analysis based on domain assignment:
!  - cube file generation
!  - electron pair spin correlation analysis
!  - export xyz files e.g. for correlation of positions in electron pairs
!  - energy partitioning
!
! initial version: AL 8/2014
!

module maxBasins_m
   use kinds_m, only: r8, i8
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use error_m
   use global_m
   use wfData_m, only: ncenter,atoms,nalpha
   use randomWalker_m
   use energyPart_m, only: epart_init, epart_add, epart_destroy, epart_write
   use parsing_m, only: finda, getdbla, getinta, getstra
   use gcube_m
   use refVector_m
   use statistics
   use mpiInterface_m, only: myMPIGatherDouble, myMPIGatherInteger, myMPIBarrier

   implicit none
   private
   public :: maxbas_init, maxbas_destroy, maxbas_IsInitialized, maxbas_writeParams, maxbas_add, maxbas_writeResults

   type(gcube), allocatable        :: mgcube(:,:)                               ! cube data structure for electron and several maxima
   real(r8)                          :: mGridWidth=0.d0                           ! *half* the grid width
   integer                         :: mNBin=0                                   ! # bins in each dimension
   integer                         :: mKMax=0                                   ! # of maxima for which to collect domains d
   integer                         :: mNCore=0                                  ! # of alpha core electrons (=# beta core electrons)
   integer                         :: mVerbose=0
   real(r8)                          :: lraconv=0.001
   logical                         :: mWriteXYZ = .false.
   logical                         :: mWriteMax = .false.
   logical                         :: mWriteCube = .false.
   logical                         :: mWriteIdx = .false.
   logical                         :: mDoEPart = .false.
   logical                         :: mPsi2Ratio = .false.
   logical                         :: mDoLRAnalysis = .false.                   ! trigger left-right analysis of electron positions
   integer                         :: mAnalyseEPairs=0                          ! >0 trigger electron pair analysis
   integer, parameter              :: mIU=90
   integer, parameter              :: mIU1=91
   integer, parameter              :: mIU2=92
   integer, parameter              :: mIU3=93
   integer, parameter              :: mIU4=94
   integer                         :: mLRPairCount=0
   integer, allocatable            :: mCount(:),mIterations(:),mLRPairs(:,:)
   real(r8), allocatable             :: lrx(:,:),lry(:,:),lrz(:,:)
   integer                         :: mLCount=0                                 ! counter for walkers in optional xyz output
   integer                         :: mIdxCount = 0                             ! counter for idx (= max permutation) output
   type(simpleStat), allocatable   :: mEPairStat(:,:)
   type(simpleStat), allocatable   :: mPsiRatio(:)                              ! statistic for f/f0
   type(vectorStat), allocatable   :: mLRStat(:,:)
   type(vectorStat), allocatable   :: mLRVecStat(:,:,:)
   type(matrixStat), allocatable   :: mLRMatStat(:,:)
   real(r8), allocatable             :: mTensor(:,:,:,:)
   real(r8), allocatable             :: mEigenValues(:,:,:)
   real(r8), allocatable             :: mLRmean(:,:,:)
   real(r8), allocatable             :: mLRResult(:,:)                            ! 1) total count 2) ll/rr count, 3) lr/rl count
   type(matrixStat), allocatable   :: mCovMat(:,:,:)

   type(referenceVector), pointer  :: mRV_p => null()                         ! vector of references for domain assignment

contains

   subroutine maxbas_init(lines,nl)
      integer, intent(in)             :: nl
      character(len=120), intent(in)  :: lines(:)
      integer                         :: i,iflag,iflag1,iflag2,iflag3,j,a,sbin
      real(r8)                          :: ax,bx,ay,by,az,bz
      real(r8)                          :: cx(ncenter),cy(ncenter),cz(ncenter)
      real(r8)                          :: tolSim
      integer                         :: elemIdx(ncenter)
      logical                         :: found
      character(len=3) ic1,ic2
      character(len=80) line1,line2
      character(len=40)               :: reffile,exclfile,lralist

      mVerbose = logmode
      call getinta(lines,nl,'verbose=',mVerbose,iflag)
      mWriteXYZ = finda(lines,nl,'xyz')
      mWriteIdx = finda(lines,nl,'idx')
      mWriteMax = finda(lines,nl,'write_max')
      call getinta(lines,nl,'elec_pair_analysis=',mAnalyseEPairs,iflag)
      mWriteCube = finda(lines,nl,'cubes')
      mDoEPart = finda(lines,nl,'energies')
      mPsi2Ratio = finda(lines,nl,'psi2ratio')
      mDoLRAnalysis = finda(lines,nl,'lra')
      call getdbla(lines,nl,'lraconv=',lraconv,iflag)

      call getinta(lines,nl,'kmax=',mKMax,iflag)  ! # of maxima for which to collect max domains
      if (iflag /= 0) call error("$init_basin_analysis: kmax required")

      tolSim = 1.d-1
      call getdbla(lines,nl,'tol_sim=',tolSim,iflag)      ! similar structure if max distance is smaller
      tolSim = tolSim / bohr2angs

      call getstra(lines,nl,'ref_file=',reffile,iflag1)
      call getstra(lines,nl,'excl_file=',exclfile,iflag2)
      call getstra(lines,nl,'lrafile=',lralist,iflag3)

      if (iflag1 /= 0) call error("$init_basin_analysis: ref_file argument required")
      allocate(mRV_p)
      if (iflag2 == 0) then
         call refv_create(mRV_p,mKMax,refFile=reffile,simThreshold=tolSim,exclFile=exclfile)
      else
         call refv_create(mRV_p,mKMax,refFile=reffile,simThreshold=tolSim)
      endif
      if (iflag3 /= 0 .and. mDoLRAnalysis) call error("$init_basin_analysis: lralist argument required")

      if (mDoEPart) call epart_init(lines,nl)

      if (mAnalyseEPairs > 0) call internal_epair_init()

      if (mPsi2Ratio) allocate(mPsiRatio(mKMax))


      if (MASTER) then

         if (mDoLRAnalysis) call internal_lr_init()

         if (mWriteCube) call internal_cube_init()

         if (mWriteXYZ) call internal_xyz_init()

         if (mWriteIdx) open(mIU1,file=trim(basename)//'.idx')

         if (mVerbose >= 2) call maxbas_writeParams(iul)

      endif

   contains

      subroutine internal_epair_init()
         integer a,j,n,k,np
         np = (ne*(ne-1))/2            ! # pairs
         allocate(mEPairStat(np,mKMax))
         do k=1,mKMax
            do n=1,np
               call reset(mEPairStat(n,k))
            enddo
         enddo
         j = 0
         do a=1,ncenter
            if (atoms(a)%za > 2 .and. (.not.atoms(a)%ecp)) then
               j = j + 1
            endif
         enddo
         mNCore = j
      end subroutine internal_epair_init

      subroutine internal_cube_init()
         mNBin = 80
         call getinta(lines,nl,'nbin=',mNBin,iflag)    ! # grid points in each dimension
         sbin=0
         call getinta(lines,nl,'sbin=',sbin,iflag)
         call getdbla(lines,nl,'grid=',mGridWidth,iflag)
         if (iflag /= 0) then
            call getdbla(lines,nl,'ax=',ax,iflag)
            if (iflag /= 0) call error("(maxbas_init): either grid or ax,bx,.. required for max domain calculation")
            call getdbla(lines,nl,'bx=',bx,iflag)
            call getdbla(lines,nl,'ay=',ay,iflag)
            call getdbla(lines,nl,'by=',by,iflag)
            call getdbla(lines,nl,'az=',az,iflag)
            call getdbla(lines,nl,'bz=',bz,iflag)
            if (ax > bx .or. ay > by .or. az > bz) call error('(maxbas_init): cube definition requires a < b')
         else
            ax = -1d0*mGridWidth; ay = -1d0*mGridWidth; az = -1d0*mGridWidth
            bx = mGridWidth; by = mGridWidth; bz = mGridWidth
         endif
         !!! NOTE: atoms_getBox is available to calculate optimal box

         allocate(mgcube(ne,mKMax))

         ! molecule data for cube files
         do i=1,ncenter
            cx(i) = atoms(i)%cx; cy(i) = atoms(i)%cy; cz(i) = atoms(i)%cz
            ElemIdx(i)=atoms(i)%elemIdx
         enddo

         ! initialize cube data structures
         line1 = 'cube file generated from max domain module'
         do j=1,mKMax
            do i=1,ne
               call mgcube(i,j)%initBins(mNBin,ax,bx,ay,by,az,bz,sbin)
               call mgcube(i,j)%setMolecule(cx,cy,cz,elemIdx)
               write(ic1,'(I3)') i
               ic1 = adjustl(ic1)
               write(ic2,'(I3)') j
               ic2 = adjustl(ic2)
               line2 = ' max='//trim(ic2)//'  electron pos='//trim(ic1)
               call mgcube(i,j)%setLines(line1,line2)
            end do
         end do
      end subroutine internal_cube_init

      subroutine internal_xyz_init()
         open(mIU,file=trim(basename)//'.xyz')
         write(mIU,'(i5,a)') ncenter,' nuclei:'
         do i=1,ncenter
             write(mIU,'(i4,1x,a2,1x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
         enddo
         write(mIU,'(a)') ' UNKNOWN_SIZE  xyz'
         mLCount = 0
      end subroutine internal_xyz_init

      subroutine internal_lr_init()
      integer                 :: i,j,io,nWords
      character(len=20) words(12)
      character(len=120) line

         allocate(mLRPairs(ne,2))
         mLRPairs = 0
         ! read lralist and write pair array
         open(mIU4,file=lralist,status='old',iostat=io)
            call assert(io==0,' lralist file '//trim(lralist)//' does not exist')
            do
               read(mIU4,'(a)',iostat=io) line
               if (io /= 0) exit
               mLRPairCount = mLRPairCount + 1
               call tokenize(line,words,nWords)
               read(words(1),*) mLRPairs(mLRPairCount,1)
               read(words(2),*) mLRPairs(mLRPairCount,2)
            enddo
         close(mIU4)

         allocate(mCount(mKMax))
         mCount = 0
         allocate(mIterations(mKMax))
         mIterations = 0
         allocate(lrx(mKMax,ne))
         allocate(lry(mKMax,ne))
         allocate(lrz(mKMax,ne))
         lrx=0d0
         lry=0d0
         lrz=0d0
         allocate(mLRmean(mKMax,mLRPairCount,3))
         allocate(mLRStat(mKMax,mLRPairCount))
         allocate(mLRVecStat(mKMax,mLRPairCount,2))
         allocate(mLRMatStat(mKMax,mLRPairCount))
         allocate(mTensor(mKMax,mLRPairCount,3,3))
         allocate(mEigenValues(mKMax,mLRPairCount,3))
         allocate(mCovMat(mKMax,mLRPairCount,3))
         allocate(mLRResult(mLRPairCount,5))
         mLRResult = 0d0
         do i=1,mKMax
            do j=1,mLRPairCount
               call mLRMatStat(i,j)%create(3,3)
               call mLRStat(i,j)%create(3)
               call mLRVecStat(i,j,1)%create(3)
               call mLRVecStat(i,j,2)%create(3)
               call mCovMat(i,j,1)%create(2,2)
               call mCovMat(i,j,2)%create(2,2)
               call mCovMat(i,j,3)%create(2,2)
            enddo
         enddo
      end subroutine internal_lr_init

   end subroutine maxbas_init


   subroutine maxbas_destroy()
      integer i,j
      if (mAnalyseEPairs > 0) then
         deallocate(mEPairStat)
      endif
      if (mDoEPart) call epart_destroy()
      if (.not.MASTER) return
      if (mWriteCube) then
         do j=1,mKMax
            do i=1,ne
               call mgcube(i,j)%destroy()
            enddo
         enddo
         deallocate(mgcube)
      endif
      mNBin = 0
      mKMax = 0
      if (mPsi2Ratio) deallocate(mPsiRatio)
      if (mWriteXYZ) close(mIU)
      if (mWriteIdx) close(mIU1)
      if (associated(mRV_p)) then
         call mRV_p%destroy()
         deallocate(mRV_p)
      endif
   end subroutine maxbas_destroy


   function maxbas_isInitialized() result(res)
      logical res
      res = associated(mRV_p)
   end function maxbas_isInitialized


   subroutine maxbas_writeParams(iu)
      integer, intent(in) :: iu
      character(len=120)    :: s
      if (.not.MASTER) return
      write(iu,'(A/)') "    maximum domain sampling parameters:"
      if (mWriteCube) then
         write(iu,'(1X,A21,F10.4,5X,A21,I4)') "grid width =",2*mGridWidth,"# bins =",mNbin
         write(iu,'(1X,A21,I4)') "# of maxima =",mKMax
      endif
      if (mWriteXYZ) then
         write(iu,'(A)') "    writing .xyz file of assigned walker coordinates"
      endif
      if (mWriteIdx) then
         write(iu,'(A)') "    writing .idx file of assignment array for each walker"
      endif
      if (mAnalyseEPairs>0) then
         write(iu,'(A,i2)') "    writing -epairs.txt file of electron pair statistic with mode=",mAnalyseEPairs
      endif
      if (mDoEPart) then
         write(iu,'(A,i3)') "    energy partitioning with # references = ",mKMax
      endif
      write(iu,*)
      if (associated(mRV_p)) then
         call mRV_p%getInfoMessage(1,s)
         write(iu,'(A)') trim(s)
         call mRV_p%getInfoMessage(2,s)
         write(iu,'(A)') trim(s)
      endif
      write(iu,*)
   end subroutine maxbas_writeParams


   subroutine maxbas_add(x,y,z,f,rw,idx1,iflag)
      real(r8), intent(in)             :: x(:),y(:),z(:)
      real(r8), intent(in)             :: f
      type(Randomwalker)             :: rw
      integer, intent(in)            :: idx1(:), iflag
      integer                        :: idx2(ne),asgn(ne)
      real(r8)                         :: x0(ne),y0(ne),z0(ne)
      integer l,n
      real(r8) f0
      type(reference)     :: r

      call assert(associated(mRV_p),"maxbas_add: mRV_p not associated")

      ! parallel part

      if (((mDoEPart .or. mAnalyseEPairs==2) .and. .not. ieee_is_nan(f)) .and. iflag==0) then

         call r%new(ne)
         r%x = x
         r%y = y
         r%z = z
         r%f = 0.d0
         call mRV_p%getIdxAndPermutation(r,n,idx2)

         if (n>0) then    ! ignore maxima not assigned to a reference

            ! concatenate core and non-core permutation
            if (mDoEPart) then
               do l=1,ne
                  asgn(l) = idx1(idx2(l))
               end do
               call epart_add(rw,asgn,n)
            endif

            ! collect electron pairs for maximum n
            if (mAnalyseEPairs==2) call maxbas_collectEPairs(n,idx2)

         endif

         call r%destroy()

      endif

      ! serial part

      if (mWriteCube .or. mWriteIdx .or. mPsi2Ratio &
          .or. mAnalyseEPairs==1 .or. mWriteMax .or. mWriteXYZ .or. mDoLRAnalysis) then

          call pos(rw,x0,y0,z0)
          f0 = -lnpsi2(rw)
          if (mWriteIdx) then
            call maxbas_addIdxSerially(x,y,z,x0,y0,z0,f,f0,idx1,iflag)
          else 
            call maxbas_addSerially(x,y,z,x0,y0,z0,f,f0,idx1,iflag)
          end if 
      endif

   end subroutine maxbas_add



   subroutine maxbas_addSerially(x,y,z,x0,y0,z0,f,f0,idx1,iflag)
      ! collect maxima on MASTER, determine reference and permutation
      ! add to cube, or write to max, idx files, or ...
      real(r8),intent(in)   :: x(:),y(:),z(:)      ! max coords
      real(r8),intent(in)   :: x0(:),y0(:),z0(:)   ! original coords
      real(r8), intent(in)  :: f,f0                ! function values (-ln(psi^2))
      integer, intent(in) :: idx1(:), iflag
      integer             :: allidx1(nproc*ne)
      real(r8)              :: vec(6*ne+2)
      real(r8)              :: vec_rcv((6*ne+2)*nproc)
      integer             :: ierr,i,n,m,l,k,j
      integer             :: idx2(ne),idx(ne),int_rcv(nproc)
      type(reference)     :: r
      real(r8), parameter   :: EPS = 1.d-4
      real(r8)              :: temp1(3),temp2(3),tempMat(2,2)

      vec(1:ne) = x(1:ne)
      vec(ne+1:2*ne) = y(1:ne)
      vec(2*ne+1:3*ne) = z(1:ne)
      vec(3*ne+1:4*ne) = x0(1:ne)
      vec(4*ne+1:5*ne) = y0(1:ne)
      vec(5*ne+1:6*ne) = z0(1:ne)
      vec(6*ne+1) = f
      vec(6*ne+2) = f0
      vec_rcv = 0d0
      n = 6*ne+2

      ! gather all coordinates to MASTER
      call myMPIGatherDouble(vec,n,vec_rcv,ierr)
      call myMPIGatherInteger(iflag,1,int_rcv,ierr)

      ! gather all permutation indices to MASTER
      call myMPIGatherInteger(idx1,ne,allidx1,ierr)

      if (MASTER) then
         call r%new(ne)
         do i=0,nproc-1
            r%x = vec_rcv(i*n+1:i*n+ne)
            r%y = vec_rcv(i*n+ne+1:i*n+2*ne)
            r%z = vec_rcv(i*n+2*ne+1:i*n+3*ne)
            r%f = vec_rcv(i*n+n-1)
            if (int_rcv(i+1) == 0) then
               call mRV_p%getIdxAndPermutation(r,m,idx2)

               ! discard maxima not belonging to any of the reference max
               if (m > 0) then

                  ! concatenation of core assignment permutation (idx1) and assignment to ref permutation (idx2)
                  ! select correct idx1 vector in allidx1
                  do l=1,ne
                     idx(l) = allidx1(i*ne+idx2(l))
                  end do

                  if (mVerbose>=4) then
                     write(999,*) ' psimax_add_to_domain_list: proc=',i,'  max=',m
                     do l=1,ne
                        write(999,'(i5,3f15.6)') l,vec_rcv(i*n+3*ne+l),vec_rcv(i*n+4*ne+l),vec_rcv(i*n+5*ne+l)
                     end do
                     write(999,'(50i3)') (/ (k,k=1,ne) /)
                     write(999,'(50i3)') allidx1(i*ne+1:i*ne+ne)
                     write(999,'(50i3)') idx2(1:ne)
                     write(999,'(50i3)') idx(1:ne)
                     do l=1,ne
                        write(999,'(i5,3f15.6)') l,vec_rcv(i*n+l),vec_rcv(i*n+ne+l),vec_rcv(i*n+2*ne+l)
                     end do
                     call internal_printDistToPermRef()
                  end if

                  if (mWriteMax) then
                     write(999,*) ' Max:',m,'  0'
                     write(999,*) ne
                     do l=1,ne
                        write(999,'(3f15.6)') vec_rcv(i*n+l),vec_rcv(i*n+ne+l),vec_rcv(i*n+2*ne+l)
                     end do
                     write(999,'(50i4)') idx(1:ne)
                  endif

                  ! add permuted walker positions to domain sampling

                  if (mWriteCube) then
                     call internal_addCube(m,vec_rcv(i*n+3*ne+1:i*n+4*ne),vec_rcv(i*n+4*ne+1:i*n+5*ne), &
                                           vec_rcv(i*n+5*ne+1:i*n+6*ne),idx,1.d0)
                  endif

                  if (mWriteXYZ) then
                     call internal_addXYZ(m,vec_rcv(i*n+3*ne+1:i*n+4*ne),vec_rcv(i*n+4*ne+1:i*n+5*ne), &
                                          vec_rcv(i*n+5*ne+1:i*n+6*ne),idx)
                  endif

                  ! perform left-right analysis
                  if (mDoLRAnalysis) then
                     if (mCount(m) >= 0) then
                        do k=1,ne
                           lrx(m,k) = vec_rcv(i*n+3*ne+idx(k))
                           lry(m,k) = vec_rcv(i*n+4*ne+idx(k))
                           lrz(m,k) = vec_rcv(i*n+5*ne+idx(k))
                        end do
                           call maxbas_collectLR(lrx(m,:),lry(m,:),lrz(m,:),m,mTensor(m,:,:,:),mEigenValues(m,:,:))
                     else
                        do k=1,mLRPairCount
                           mLRResult(k,1) = mLRResult(k,1) + 1

                           temp1 = (mLRVecStat(m,k,1)%localmean() + mLRVecStat(m,k,2)%localmean())/2
                           mLRmean(m,k,:) = temp1

                           do l=1,3
                              temp1(l) = vec_rcv(i*n+(2+l)*ne+idx(mLRPairs(k,1))) - mLRmean(m,k,l)
                              temp2(l) = vec_rcv(i*n+(2+l)*ne+idx(mLRPairs(k,2))) - mLRmean(m,k,l)
                           enddo

                           temp1(:) = matmul(mTensor(m,k,:,:),temp1(:))
                           temp2(:) = matmul(mTensor(m,k,:,:),temp2(:))

                           do l=1,3
                              tempMat(1,1) = temp1(l)*temp1(l)
                              tempMat(1,2) = temp1(l)*temp2(l)
                              tempMat(2,2) = temp2(l)*temp2(l)
                              tempMat(2,1) = tempMat(1,2)
                              call mCovMat(m,k,l)%add(tempMat)
                           enddo

                           call mLRStat(m,k)%add(temp1)
                           call mLRStat(m,k)%add(temp2)

                           if (((temp1(3) .gt. 0) .and. (temp2(3) .gt. 0)) .or. ((temp1(3) .lt. 0) .and. (temp2(3) .lt. 0))) then
                              mLRResult(k,2) = mLRResult(k,2) + 1
                           else
                              mLRResult(k,3) = mLRResult(k,3) + 1
                           endif
                        enddo
                     endif
                  endif

                  ! add ratio psi**2/psi0**2 to stat
                  if (mPsi2Ratio) call addData(mPsiRatio(m),exp(vec_rcv(i*n+n-1) - vec_rcv(i*n+n)))

                  ! collect electron pairs for maximum m
                  if (mAnalyseEPairs==1) call maxbas_collectEPairs(m,idx2)
               endif

            endif
         end do
         call r%destroy()
      end if

      call myMPIBarrier(ierr)

   contains

      subroutine internal_addCube(m,x,y,z,idx,val)
         integer    :: m                ! maximum
         real(r8)     :: x(:),y(:),z(:)   ! electron coords
         integer    :: idx(:)           ! permutation
         real(r8)     :: val              ! walker weight
         integer i,ii

         do i=1,ne
            ii = idx(i)
            call mgcube(i,m)%addData(x(ii),y(ii),z(ii),val)
            !!!write(999,'(2i5,3f15.6)') i,ii,x(ii),y(ii),z(ii)
         end do
      end subroutine internal_addCube

      subroutine internal_addXYZ(m,x,y,z,idx)
         integer    :: m                ! maximum
         real(r8)     :: x(:),y(:),z(:)   ! electron coords
         integer    :: idx(:)           ! permutation
         integer i,ii

         mLCount = mLCount + 1
         write(mIU,'(a,i5,i9)') 'xyz:  ',m,mLCount
         write(mIU,'(i5)') ne
         do i=1,ne
            ii=idx(i)
            write(mIU,'(3f15.6,i6)') x(ii),y(ii),z(ii),ii
         enddo
      end subroutine internal_addXYZ

      subroutine internal_printDistToPermRef()
         type(reference), pointer                 :: rp
         integer j,jj
         real(r8) xx,yy,zz,xx1,yy1,zz1,d,dsum,d1sum

         ! access reference positions
         call mRV_p%getRefPtr(m,rp)
         write(999,*) 'proc ',i,' diff to ref:',m
         write(999,*) ' i idx   distance           x             y              z'
         dsum = 0.d0
         d1sum = 0.d0
         do j=1,ne
            jj = idx(j)
            xx = vec_rcv(i*n+3*ne+jj)
            yy = vec_rcv(i*n+4*ne+jj)
            zz = vec_rcv(i*n+5*ne+jj)
            d = sqrt( (rp%x(j)-xx)**2 + (rp%y(j)-yy)**2 + (rp%z(j)-zz)**2 )
            dsum = dsum + d
            ! walker: orig idx, distance to ref, walker elec coords
            write(999,'(2i3,4f15.6)') j,jj,d,xx,yy,zz
            jj = idx2(j)
            xx1 = vec_rcv(i*n+jj)
            yy1 = vec_rcv(i*n+ne+jj)
            zz1 = vec_rcv(i*n+2*ne+jj)
            d = sqrt( (xx1-xx)**2 + (yy1-yy)**2 + (zz1-zz)**2 )
            d1sum = d1sum + d
            ! maximized walker: orig idx (from nextpsimax), distance to ref, maximum elec coords
            write(999,'(3x,i3,4f15.6)') jj,d,xx1,yy1,zz1
            ! reference electron coords
            write(999,'(21x,3f15.6)') rp%x(j), rp%y(j), rp%z(j)
         end do
         write(999,'(2(a,f15.6))') 'dmean(walker-ref)=',dsum/ne, ' mean distance(walker-maximum)=',d1sum/ne
         write(999,*) '-------------------------------------------'
      end subroutine internal_printDistToPermRef

   end subroutine maxbas_addSerially


   subroutine maxbas_addIdxSerially(x, y, z, x0, y0, z0, f, f0, idx1, iflag)
      ! collect maxima on MASTER, determine reference and permutation
      ! add to cube, or write to max, idx files, or ...
      real(r8),intent(in)   :: x(:), y(:), z(:)      ! max coords
      real(r8),intent(in)   :: x0(:), y0(:), z0(:)   ! original coords
      real(r8), intent(in)  :: f, f0                 ! function values (-ln(psi^2))
      integer, intent(in) :: idx1(:), iflag
      integer             :: allidx1(nproc*ne)
      real(r8)              :: vec(6*ne+2)
      real(r8)              :: vec_rcv((6*ne+2)*nproc)
      integer             :: ierr,i,n,m,l,k,j
      integer             :: idx2(ne), idx(ne), int_rcv(nproc)
      type(reference)     :: r
      real(r8), parameter   :: EPS = 1.d-4
      real(r8)              :: temp1(3),temp2(3),tempMat(2,2)

      vec(1:ne) = x(1:ne)
      vec(ne+1:2*ne) = y(1:ne)
      vec(2*ne+1:3*ne) = z(1:ne)
      vec(3*ne+1:4*ne) = x0(1:ne)
      vec(4*ne+1:5*ne) = y0(1:ne)
      vec(5*ne+1:6*ne) = z0(1:ne)
      vec(6*ne+1) = f
      vec(6*ne+2) = f0
      vec_rcv = 0d0
      n = 6*ne+2

      ! gather all coordinates to MASTER
      call myMPIGatherDouble(vec,n,vec_rcv,ierr)
      call myMPIGatherInteger(iflag,1,int_rcv,ierr)

      ! gather all permutation indices to MASTER
      call myMPIGatherInteger(idx1,ne,allidx1,ierr)

      if (MASTER) then
         call r%new(ne)
         do i = 0, nproc - 1
            r%x = vec_rcv(i*n+1:i*n+ne)
            r%y = vec_rcv(i*n+ne+1:i*n+2*ne)
            r%z = vec_rcv(i*n+2*ne+1:i*n+3*ne)
            r%f = vec_rcv(i*n+n-1)
            if (int_rcv(i + 1) == 0) then
               call mRV_p%getIdxAndPermutation(r, m, idx2)
               ierr = 0
               do l = 1, ne
                  idx(l) = allidx1(i*ne + idx2(l))
               end do
            else 
               idx = [ (k, k = 1, ne) ]
               m = 0
               ierr = int_rcv(i + 1)
            end if

            mIdxCount = mIdxCount + 1

            write(mIU1,'(i7,2i5,f20.4)') mIdxCount, m, ierr, r%f
            write(mIU1,'(5x,50i4)') idx
         end do

      end if

   end subroutine maxbas_addIdxSerially



   subroutine maxbas_writeResults()
      character(len=3)                          :: ic1,ic2
      character(len=80)                         :: filename
      integer                                   :: i,n,m,j,np,k,l
      integer(i8)                                 :: nd,mynd
      real(r8)                                    :: pairRatio,stddev,mymean,mystddev
      real(r8)                                    :: tempMat(3,2,2)
      real(r8)                                    :: tempMat0(2,2)
      real(r8), parameter                         :: PAIRTHRESH = 0.95d0, WEAKPAIR=0.1d0, &
                                                   WEAKTRIP=-0.1d0, TRIPTHRESH=-0.95d0

      if (mWriteCube .and. MASTER) then
         do n=1,mKMax
            write(ic1,'(I3)') n
            ic1 = adjustl(ic1)
            do i=1,ne
               write(ic2,'(I3)') i
               ic2 = adjustl(ic2)
               filename = trim(basename)//'-'//trim(ic1)//'-'//trim(ic2)//'.cube'
               call mgcube(i,n)%writeToFile(filename)
            end do
         end do
      endif

      if (mAnalyseEPairs > 0) then
         if (MASTER) open(mIU2,file=trim(basename)//'-epa.txt')
         do m=1,mKMax
            if (mAnalyseEPairs==1) then
               nd = dataCount(mEPairStat(1,m))
            else if (mAnalyseEPairs==2) then
               nd = dataCountAllNodes(mEPairStat(1,m))
            endif
            if (MASTER) write(mIU2,'(a,i4,a,i8,a,i2)') ' Maximum ',m,' data Count=',nd,' mode=',mAnalyseEPairs
            np = 0
            do i=1,nalpha
               do j=i+1,nalpha
                  np = np + 1
                  call internal_writePair()
               enddo
               do j=nalpha+1,ne
                  np = np + 1
                  call internal_writePair()
               enddo
            enddo
            do i=nalpha+1,ne-1
               do j=i+1,ne
                  np = np + 1
                  call internal_writePair()
               enddo
            enddo
         enddo
         if (MASTER) close(mIU2)
      endif

      if (mDoEPart) then
         call epart_write(iul,1)
         open(42,file=trim(basename)//'.see')
         call epart_write(42,2)
         close(42)
      endif

      if (mPsi2Ratio .and. MASTER) then
         write(iul,'(a)') "    <psi**2/psi(max)**2> for each max structure:"
         do i=1,mKMax
            if (dataCount(mPsiRatio(i))>=2) then
               write(iul,'(i4,g11.2,a5,g11.2)') i,mean(mPsiRatio(i)),' +/- ',sqrt(variance(mPsiRatio(i)))
            endif
         enddo
      endif

      if(MASTER .and. mDoLRAnalysis) then
         open(mIU3,file=trim(basename)//'-lra.txt')
                  write(mIU3,'(a,f10.8)') 'LR Analysis with convergence criterion (lraconv) @ ',lraconv
         do i=1,mKMax
            if (mLRStat(i,1)%localcount() > 0) then
                  write(mIU3,'(a,i7,a)') 'Convergence reached after ',mIterations(i),' collected maxima.'
                  write(mIU3,'(a,i3)') 'Maximum ',i
                  write(mIU3,'(i7)') mLRPairCount
               do j=1,mLRPairCount
                  write(mIU3,'(4i3)') j,mLRPairs(j,1),mLRPairs(j,2),int(mLRResult(j,1))
                  write(mIU3,'(a12,3f10.6)') 'Center      ',mLRmean(i,j,1),mLRmean(i,j,2),mLRmean(i,j,3)
                  write(mIU3,'(a12,3f10.6)') 'Mean        ',mLRStat(i,j)%localmean()
                  write(mIU3,'(a12,3f10.6)') 'Var         ',mLRStat(i,j)%localvar()
                  write(mIU3,'(a12,3f10.6)') 'StdDev      ',mLRStat(i,j)%localstddev()
                  write(mIU3,'(a12,3f10.6)') 'EigenVal    ',mEigenValues(i,j,1),mEigenValues(i,j,2),mEigenValues(i,j,3)
                  write(mIU3,'(a12,3f10.6)') 'CovMatEVec  ',(mTensor(i,j,l,1), l = 1,3)
                  do k = 2,3
                     write(mIU3,'(a12,3f10.6)') ' ',(mTensor(i,j,l,k), l = 1,3)
                  enddo
                  do l=1,3
                     tempMat0 = mCovMat(i,j,l)%localmean()
                     tempMat(l,:,:) = tempMat0
                  enddo
                  write(mIU3,'(a12,3f10.6)') 'Rho X Y Z   ',tempMat(1,1,2)/(sqrt(tempMat(1,1,1))*sqrt(tempMat(1,2,2))), &
                            tempMat(2,1,2)/(sqrt(tempMat(2,1,1))*sqrt(tempMat(2,2,2))), &
                            tempMat(3,1,2)/(sqrt(tempMat(3,1,1))*sqrt(tempMat(3,2,2)))
                  write(mIU3,'(a12,f10.6,a5,f8.6,a2,f8.6,a1)') 'LR/RL Ratio ',mLRResult(j,3)/mLRResult(j,1), &
                            ' +/- ',getAgrestiCoullInterval(mLRResult(j,1),mLRResult(j,3),dble(2.5758293035489)), &
                            ' (',getAgrestiCoullInterval(dble(mIterations(i)),&
                             dble((mLRResult(j,3)/mLRResult(j,1))*mIterations(i)),dble(2.5758293035489)),')'
                  write(mIU3,'(a12,f10.6,a5,f8.6,a2,f8.6,a1)') 'LL/RR Ratio ',mLRResult(j,2)/mLRResult(j,1), &
                            ' +/- ',getAgrestiCoullInterval(mLRResult(j,1),mLRResult(j,2),dble(2.5758293035489)), &
                            ' (',getAgrestiCoullInterval(dble(mIterations(i)),&
                             dble((mLRResult(j,2)/mLRResult(j,1))*mIterations(i)),dble(2.5758293035489)),')'
               enddo
               write(mIU3,*) '----------------------------------------------------------------'
            endif
         enddo
         close(mIU3)
      endif

   contains

      subroutine internal_writePair()
        if (mAnalyseEPairs==1) then
          if (nd > 1) then
            pairRatio = mean(mEPairStat(np,m))
            stddev = stdDevMean(mEPairStat(np,m))
            write(mIU2,'(3i4,2f8.3,3x)',advance='no') i,j,np,pairRatio,stddev
            if (pairRatio > PAIRTHRESH) then
              write(mIU2,*) 'ELECTRON_PAIR!'
            else if (pairRatio > WEAKPAIR) then
              write(mIU2,*) 'weak_pair'
            else if (pairRatio > WEAKTRIP) then
              write(mIU2,*) 'no_correlation'
            else if (pairRatio > TRIPTHRESH) then
              write(mIU2,*) 'weak_triplet'
            else
              write(mIU2,*) 'triplet'
            endif
          else if (nd==1) then
            pairRatio = mean(mEPairStat(np,m))
            write(mIU2,'(3i4,f7.3)') i,j,np,pairRatio
          else
            write(mIU2,'(3i4,a)') i,j,np,' no data!'
          endif
        else if (mAnalyseEPairs==2) then
          if (nd > 1) then
            !!mynd = dataCount(mEPairStat(np,m))
            !!mymean = 0.d0
            !!if (mynd > 0) mymean = mean(mEPairStat(np,m))
            !!mystddev = 0.d0
            !!if (mynd > 1) mystddev = stdDevMean(mEPairStat(np,m))
            !!write(900+mytid,'(4i5,3f20.3)') mytid,i,j,np,mynd,mymean,mystddev
            !!flush(900+mytid)
            pairRatio = meanAllNodes(mEPairStat(np,m))
            stddev = stdDevMeanAllNodes(mEPairStat(np,m))
            !!write(900+mytid,'(4i5,3f20.3)') mytid,i,j,np,nd,mean,stddev
            !!flush(900+mytid)
            if (MASTER) then
              write(mIU2,'(3i4,2f8.3,3x)',advance='no') i,j,np,pairRatio,stddev
              if (pairRatio > PAIRTHRESH) then
                write(mIU2,*) 'ELECTRON_PAIR!'
              else if (pairRatio > WEAKPAIR) then
                write(mIU2,*) 'weak_pair'
              else if (pairRatio > WEAKTRIP) then
                write(mIU2,*) 'no_correlation'
              else if (pairRatio > TRIPTHRESH) then
                write(mIU2,*) 'weak_triplet'
              else
                write(mIU2,*) 'triplet'
              endif
            endif
          else if (nd==1) then
            pairRatio = meanAllNodes(mEPairStat(np,m))
            if (MASTER) write(mIU2,'(3i4,f7.3)') i,j,np,pairRatio
          else
            if (MASTER) write(mIU2,'(3i4,a)') i,j,np,' no data!'
          endif
        endif
      end subroutine
   end subroutine maxbas_writeResults


   subroutine maxbas_collectEPairs(m,idx)
       integer, intent(in)    :: m                ! maximum
       integer, intent(in)    :: idx(:)           ! permutation
       integer i,j,ii,jj,itmp,np
       logical spinPair

       !!!write(998,'(a,i4)') 'new max ',m
       np = 0
       do i=1,nalpha
         ii = idx(i)
         do j=i+1,nalpha
           jj = idx(j)
           spinPair = ( ii <= nalpha .eqv. jj > nalpha )
           np = np + 1
           !!!np = ((j-1)*(j-2))/2 + i          ! idx in pair list
           if (spinPair) then
             !!!write(998,'(5i4,a)') i,j,ii,jj,np,' pair'
             call addData(mEPairStat(np,m),1.d0)
           else
             !!!write(998,'(5i4,a)') i,j,ii,jj,np,' no pair'
             call addData(mEPairStat(np,m),-1.d0)
           endif
         enddo
         do j=nalpha+1,ne
           jj = idx(j)
           spinPair = ( ii <= nalpha .eqv. jj > nalpha )
           np = np + 1
           !!!np = ((j-1)*(j-2))/2 + i          ! idx in pair list
           if (spinPair) then
             !!!write(998,'(5i4,a)') i,j,ii,jj,np,' pair'
             call addData(mEPairStat(np,m),1.d0)
           else
             !!!write(998,'(5i4,a)') i,j,ii,jj,np,' no pair'
             call addData(mEPairStat(np,m),-1.d0)
           endif
         enddo
       enddo
       do i=nalpha+1,ne-1
         ii = idx(i)
         do j=i+1,ne
           jj = idx(j)
           spinPair = ( ii <= nalpha .eqv. jj > nalpha )
           np = np + 1
           !!!np = ((j-1)*(j-2))/2 + i          ! idx in pair list
           if (spinPair) then
             !!!write(998,'(5i4,a)') i,j,ii,jj,np,' pair'
             call addData(mEPairStat(np,m),1.d0)
           else
             !!!write(998,'(5i4,a)') i,j,ii,jj,np,' no pair'
             call addData(mEPairStat(np,m),-1.d0)
           endif
         enddo
       enddo
   end subroutine maxbas_collectEPairs

   subroutine maxbas_collectLR(x,y,z,m,pitensor,pev)
      real(r8), intent(in)                         :: x(:),y(:),z(:)
      integer, intent(in)                        :: m                ! maximum
      real(r8), allocatable                        :: itensor(:,:,:)
      real(r8), intent(out)                        :: pitensor(:,:,:)
      real(r8), allocatable                        :: ev(:,:)
      real(r8), intent(out)                        :: pev(:,:)
      real(r8)                                     :: work(8),diff
      real(r8)                                     :: tempMat(3,3)
      real(r8)                                     :: tempVec(3)
      integer                                    :: alstat,i,k,l,n

      ! calculate mean and inertia tensor

      ! check for convergence
      diff = 1d0
      mIterations(m) = mIterations(m) + 1
      do i=1,mLRPairCount
         do l=1,2
            call mLRVecStat(m,i,l)%add((/ x(mLRPairs(i,l)),y(mLRPairs(i,l)),z(mLRPairs(i,l)) /))
            tempMat(1,1) = x(mLRPairs(i,l))*x(mLRPairs(i,l))
            tempMat(1,2) = x(mLRPairs(i,l))*y(mLRPairs(i,l))
            tempMat(1,3) = x(mLRPairs(i,l))*z(mLRPairs(i,l))
            tempMat(2,2) = y(mLRPairs(i,l))*y(mLRPairs(i,l))
            tempMat(2,3) = y(mLRPairs(i,l))*z(mLRPairs(i,l))
            tempMat(3,3) = z(mLRPairs(i,l))*z(mLRPairs(i,l))
            tempMat(2,1) = tempMat(1,2)
            tempMat(3,1) = tempMat(1,3)
            tempMat(3,2) = tempMat(2,3)
            call mLRMatStat(m,i)%add(tempMat)
            tempVec = mLRVecStat(m,i,l)%localstddev()
            if (mIterations(m) >= 128) then
               diff = maxval(tempVec)
            endif
         enddo
      enddo
      if (diff .le. lraconv) then
         allocate(itensor(mLRPairCount,3,3),ev(ne,3),stat=alstat)
         itensor = 0d0
         ! construct covariance matrix / inertia tensor
         do i=1,mLRPairCount
            tempVec = (mLRVecStat(m,i,1)%localmean()+mLRVecStat(m,i,2)%localmean())/2
            tempMat = mLRMatStat(m,i)%localmean()
            do k=1,3
               do n=1,3
                  itensor(i,k,n) = tempMat(k,n) - tempVec(k) * tempVec(n)
               enddo
            enddo
         ! calculate inertia tensor eigenvalues and eigenvectors (principal axes)
            call dsyev('V','U',3,itensor(i,:,:),3,ev(i,:),work(8),8,alstat)
            itensor(i,:,:) = transpose(itensor(i,:,:))
         enddo
         pitensor = itensor
         pev = ev
         ! clear data
         mCount(m) = -1
      else
         mCount(m) = 0
      endif
   end subroutine maxbas_collectLR

   function getAgrestiCoullInterval(n,k,z) result(j)
      real(r8), intent(in)   :: n,k,z
      real(r8)               :: j,dn,dp

      dn = n + z * z
      dp = (1/n) * (k + (z*z)/2)
      j = z * sqrt((dp/dn)*(1-dp))
   !Agresti, A.; Coull, B. A. The American Statistician. 1998, pp 119-126.
   end function getAgrestiCoullInterval

end module maxBasins_m
