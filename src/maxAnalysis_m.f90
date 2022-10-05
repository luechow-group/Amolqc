! Copyright (C) 2015, 2017-2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module maxAnalysis_m

   use kinds_m, only: r8
   use global_m
   use references_m
   use refADT_m
   use parsing_m
   use atom_m, only: atoms_initPositionThresholds
   use sorting_m, only: sort
   use mpiInterface_m, only: myMPIGatherDouble, myMPIGatherInteger, myMPIBarrier, myMPIBcastInteger

   implicit none
   private
   public :: maxana_init, maxana_destroy, maxana_add, maxana_writeResults,maxana_writeParams, &
             maxana_isInitialized, maxana_getDiffMax, maxana_getFirstF

   integer             :: mVerbose=0
   integer             :: mKMax = 30       ! structure list length
   integer             :: mMMax = 0        ! sub list length

   real(r8)              :: mFtol = 1.d-3    ! if the F(max) values of two maxima are <mFtol they are assumed equal
   real(r8)              :: mTolSim = 0      ! max distance for similar structures (bohr)
   real(r8)              :: mTolSame = 0     ! max distance for same structures (bohr)
   real(r8)              :: mTolDist2 = 0    ! max distance for "pos" structures (bohr)
   logical             :: mIsInitialized = .false.
   character(len=3)    :: mMaxMode = ""    ! max mode

   class(referenceContainer), pointer  :: mRC_p => null()         ! list of stored references

contains

   subroutine maxana_init(lines,nl)
   !------------------------------!

      character(len=120), intent(in) :: lines(:)
      integer, intent(in)            :: nl
      integer                        :: iflag,maxFull,iflag1,iflag2,iflag3,idx,nre,nrefs,i,j
      integer                        :: initSuccess
      integer, pointer               :: ire(:,:) => null()
      real(r8)                         :: tolDist2,tolSim,tolSame
      real(r8)                         :: nucThresh,coreThresh,bondThresh
      logical                        :: doSortFreq
      character(len=3)               :: maxMode
      character(len=40)              :: reffile,exclfile

      mVerbose = logmode
      call getinta(lines,nl,'verbose=',mVerbose,iflag)

      mKMax = 30
      call getinta(lines,nl,'kmax=',mKMax,iflag)          ! max # (main list) entries ref file

      maxFull = 5
      call getinta(lines,nl,'mmax=',maxFull,iflag)        ! max # of sublists in ref fie
      mMMax = maxFull
      mFtol = 1.d-3
      call getdbla(lines,nl,'tol_fctn=',mFtol,iflag)      ! same function value if abs diff is smaller
      tolSim = 1.d-1
      call getdbla(lines,nl,'tol_sim=',tolSim,iflag)      ! similar structure if max distance is smaller
      tolSim = tolSim / bohr2angs
      mTolSim = tolSim
      tolSame = 1.d-2
      call getdbla(lines,nl,'tol_same=',tolSame,iflag)    ! same structure if max distance is smaller
      tolSame = tolSame / bohr2angs
      mTolSame = tolSame
      tolDist2 = tolSim*bohr2angs
      call getdbla(lines,nl,'tol_simmax=',tolDist2,iflag) ! type=pos: same position if distance is smaller
      tolDist2 = tolDist2 / bohr2angs
      mTolDist2 = tolDist2

      nucThresh = 0.01d0
      bondThresh = 0.4d0
      coreThresh = 2.d0

      call getstra(lines,nl,'ref_file=',reffile,iflag1)
      call getstra(lines,nl,'excl_file=',exclfile,iflag2)

      doSortFreq = finda(lines,nl,'sort_freq')


      initSuccess = 0
      maxMode = "str"
      call getstra(lines,nl,'max_mode=',maxMode,iflag)
      mMaxMode = maxMode
      if (MASTER) then
         if  (maxMode=="val") then
            call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,elemLength=1,valueThreshold=mFtol)
         else if (maxMode=="vst") then
            call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,elemLength=maxFull,valueThreshold=mFtol, &
                                     sameThreshold=tolSame)
         else if (maxMode=="str" ) then
            if (iflag2 == 0) then
               call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,elemLength=maxFull,exclFile=exclfile, &
                               simThreshold=tolSim,sameThreshold=tolSame)
            else
               call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,elemLength=maxFull, &
                               simThreshold=tolSim,sameThreshold=tolSame)
            end if
         else if (maxMode=="ctr") then
            if (iflag1 /= 0) call abortp("$find_maxima: max_mode=ctr requires ref_file argument")
            if (iflag2 == 0) then
               call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,elemLength=maxFull,refFile=reffile, &
                               exclFile=exclfile,simThreshold=tolSim,sameThreshold=tolSame)
            else
               call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,elemLength=maxFull,refFile=reffile, &
                               simThreshold=tolSim,sameThreshold=tolSame)
            end if
         else if (maxMode=="stp") then
            call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,elemLength=maxFull,simThreshold=tolSim, &
                                     sameThreshold=tolSame)
         else if (maxMode=="pos") then
            call create_pos()     ! internal subroutine
         else
            call abortp("$init_max_search: illegal value for max_mode=[val|vst|str|ctr|stp|pos]")
         end if

         if (associated(mRC_p)) initSuccess = 1
      end if

      call myMPIBcastInteger(initSuccess,1)
      mIsInitialized = (initSuccess == 1)

      ! default values above "maxMode". Values here have preference over values in excl_file
      call getdbla(lines,nl,'nuc_thresh=',nucThresh,iflag1)
      call getdbla(lines,nl,'bond_thresh=',bondThresh,iflag2)
      call getdbla(lines,nl,'core_thresh=',coreThresh,iflag3)
      if (iflag1==0.or.iflag2==0.or.iflag3==0) then
         call atoms_initPositionThresholds(nucThresh,coreThresh,bondThresh)
      end if

      if (MASTER .and. mVerbose >= 2) call maxana_writeParams(iul)

   contains

      subroutine create_pos()
         call getstra(lines,nl,'ref_file=',reffile,iflag)
         call assert(iflag == 0,'$init_max_search: max_mode=pos requires ref_file')
         nre = 0
         call getinta(lines,nl,'ignore_ref_elecs=',nre,iflag)
         if (iflag==0 .and. nre>0) then
            do idx=1,nl
               if (index(lines(idx),'ignore_ref_elecs')>0) exit
            end do
            call assert(nl > idx,'$init_max_search: not enough lines for ignore_ref_elecs')
            idx = idx + 1
            read(lines(idx),*) nrefs
            call assert(nrefs>0,"$init_max_search: line after ignore_ref_elecs must contain # of references")
            call assert(nl >= idx+nrefs,'$init_max_search: not enough lines for ignore_ref_elecs')
            allocate(ire(nre,nrefs))
            do i=1,nrefs
               read(lines(idx+i),*) (ire(j,i),j=1,nre)
               call sort(ire(:,i))
            end do
         else
            allocate(ire(1,1))
            ire(1,1) = 0
         end if
         call getstra(lines,nl,'excl_file=',exclfile,iflag2)

         if (iflag2 == 0) then ! excl_file given

            call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,refFile=reffile,exclFile=exclfile, &
                     simThreshold=tolSim,simmaxThreshold=tolDist2)

         else if (ire(1,1) /= 0) then ! ignore reference electron list given

            call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,refFile=reffile,ire=ire, &
                                simThreshold=tolSim,simmaxThreshold=tolDist2)

         else  ! neither ire nore exclFile given

            call create_refContainer(mRC_p,maxMode,mKMax,doSortFreq,refFile=reffile, &
                                simThreshold=tolSim,simmaxThreshold=tolDist2)

         end if
      end subroutine create_pos

   end subroutine maxana_init


   subroutine maxana_destroy()
      if (associated(mRC_p)) call mRC_p%destroy()
   end subroutine maxana_destroy


   function maxana_isInitialized() result(res)
      logical res
      res = mIsInitialized
   end function maxana_isInitialized


   integer function maxana_getDiffMax()
      if (associated(mRC_p)) then
         maxana_getDiffMax = mRC_p%getSize()
      else
         maxana_getDiffMax = 0
      end if
   end function maxana_getDiffMax


   real(r8) function maxana_getFirstF()
      if (associated(mRC_p)) then
         maxana_getFirstF = mRC_p%getFirstF()
      else
         maxana_getFirstF = 0
      end if
   end function maxana_getFirstF


   subroutine maxana_writeParams(iu)
      integer, intent(in) :: iu
      character(len=240)   :: s2

      write(iu,'(A/)') '    maximum analysis parameters:'

      write(iu,'(1X,A21,3X,A3)') " maximum list mode =",mMaxMode
      write(iu,'(2(1X,A21,G12.2,3X))') " same max func tol =",mFtol
      write(iu,'(1X,A21,I6,9X,A21,I6)') 'nmax =',mKMax,' mmax=',mMMax
      write(iu,'(2(1X,A21,F12.4,2X))') " tol_sim (A) =",mTolSim*bohr2angs," tol_same (A) =",mTolSame*bohr2angs
      write(iu,'(2(1X,A21,F12.4,2X))') " tol_fctn =",mFtol," tol_simmax (A) =",mTolDist2*bohr2angs
      write(iu,*)
      if (associated(mRC_p)) then
        call mRC_p%getInfoMessage(s2)
        write(iu,'(2x,a/)') trim(s2)
      else
        write(iu,'(1x,a/)') " WARNING: reference container not allocated"
      endif
   end subroutine maxana_writeParams


   subroutine maxana_add(x,y,z,f,iflag)
     ! insert maximum x,y,z with function value f into the
     ! list of lists reference container datastructure
     ! note: only MASTER holds data structure
     real(r8),intent(in)   :: x(:),y(:),z(:)
     real(r8), intent(in)  :: f
     integer, intent(in) :: iflag    ! =0 no error in max search, >0 error in max search (see nextpsimax)
     real(r8)              :: vec(3*ne+1)
     real(r8)              :: vec_rcv((3*ne+1)*nproc)
     integer             :: int_rcv(nproc)
     integer             :: ierr,i,n
     type(reference)     :: r

     if (MASTER) call assert(associated(mRC_p),"maxana_add: reference container not allocated")

     vec(1:ne) = x(1:ne)
     vec(ne+1:2*ne) = y(1:ne)
     vec(2*ne+1:3*ne) = z(1:ne)
     vec(3*ne+1) = f
     vec_rcv = 0d0
     int_rcv = 0
     n = 3*ne+1

     call myMPIGatherDouble(vec,n,vec_rcv,ierr)
     call myMPIGatherInteger(iflag,1,int_rcv,ierr)

     if (MASTER) then
        call r%new(ne)
        do i=0,nproc-1
           r%x = vec_rcv(i*n+1:i*n+ne)
           r%y = vec_rcv(i*n+ne+1:i*n+2*ne)
           r%z = vec_rcv(i*n+2*ne+1:i*n+3*ne)
           r%f = vec_rcv(i*n+n)
           if (int_rcv(i+1) == 0) then
              call mRC_p%insert(r)
           endif
        end do
        call r%destroy()
     end if

     call myMPIBarrier(ierr)

   end subroutine maxana_add



   subroutine maxana_writeResults()
      if (MASTER) then
         call assert(associated(mRC_p),"maxana_writeResults: reference container not allocated")

         write(iul,'(2a)') ' maximum list mode: ',mMaxMode
         call mRC_p%writeShortList(iul,'Max')
         if (mRC_p%getElemLength() > 1) then
            call mRC_p%writeFullList(baseName,"Max")
         else
            call mRC_p%writeSimpleList(baseName,"Max")
         endif
      endif
   end subroutine maxana_writeResults

end module maxAnalysis_m
