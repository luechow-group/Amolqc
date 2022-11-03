! Copyright (C) 2013-2016 Kaveh Haghighi Mood
! Copyright (C) 2014-2015, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optParamsVNL2SOL_m

   ! implements nl2sol variance minimization of a fixed sample
   ! without weighting

   use kinds_m, only: r8
   use global_m
   use subloop_m, only: subloop
   use rwSample_m
   use elocAndPsiTermsLM_m
   use wfParameters_m
   use multiDetParam_m
   use moParam_m
   use elocData_m
   use jastrowParamData_m
   use waveFunction_m, only: writeWF
   use ecp_m, only: EcpType
   use utils_m, only: intToStr
   use random_m, only: mygran
   use nl2sol_i, only: dfault, itsmry, nl2itr
   use mpiInterface_m, only: myMPIBcastInteger, myMPIBcastDouble, myMPIGatherDoubleV

   implicit none

   private
   public :: varmin2_optimizeSample

contains


   subroutine varmin2_optimizeSample(lines,nl,WFP,sample,converged)
   !--------------------------------------------------------------!

   ! optimize the variance for a fixed sample using NL2SOL

   integer, intent(in)               :: nl
   character(len=120), intent(in)    :: lines(nl)
   type(WFParamDef), pointer           :: WFP
   type(RWSample), intent(inout)     :: sample  ! (fixed) sample for optimization
   logical, intent(out)              :: converged
   integer                           :: nParam,i,ierr,sSize,rSize
   real(r8), allocatable               :: x(:)       ! parameter vector
   real(r8)                            :: eOld, varOld, varRefOld,emean
   real(r8)                            :: var,varRef, optERef, eRef
   integer, parameter  :: d = 27, j = 33, nfcall = 6      ! NL2SOL parameters
   integer, parameter  ::  nfgcal = 7, r = 50, toobig = 2 ! NL2SOL parameters
   logical             :: strted                          ! NL2SOL
   integer,allocatable ::   iv(:)
   integer             :: nf,d1, r1, j1,ncallf=0           ! NL2SOL
   real(r8),allocatable  :: v(:)
   integer             :: sampleSizes(nproc)
   integer             :: eqIter, eqStep, max_iter, iv16, npJ, npCI, npMO,optIter
   character(len=80)   :: subName,fname
   logical             :: doWriteWF, fixed
   logical             :: nliterPresent,iv16Present,eRefPresent
   logical             :: noiseInit=.false.! initialize mo params with a random noise
   real(r8)              :: noiseCoeff
   type(ElocAndPsiTermsLM)           :: EPsiTLM
   type(WFParamDerivTerms)           :: wfpDT
   type(EcpType)                     :: ecp
   logical                           :: exitSubloop = .false.
   converged = .true.   !! not yet implemented
   nliterPresent=.false.; iv16Present=.false.; eRefPresent=.true.

   call internal_readInput()
   if (logmode>=2) then
       write(iul,'(/A/A)') '   - -  varmin (nl2sol) optimization  - -'
   end if

  eqStep = 1
  do

   if(doWriteWF) call getPlusOptIter(optIter)
   call ElocAndPsiTermsLM_create(EPsiTLM,eRef,WFP)
    if (logmode>=2) then

      if (.not.fixed .and. eRefPresent) then
            write(iul,'(1X,A,F14.4)')' E_ref = ', optERef
         else if (.not.fixed) then
            write(iul,'(1X,A)')' E_ref = sample Emean'
         else
            write(iul,'(1X,A,F14.4)')' E_ref = ',optERef
         end if
   end if

   nParam = ElocAndPsiTermsLM_nParams(EPsiTLM)
   npJ  =  ElocAndPsiTermsLM_nJParams (EPsiTLM)
   npCI =  ElocAndPsiTermsLM_nCIParams(EPsiTLM)
   npMO  = ElocAndPsiTermsLM_nMOParams(EPsiTLM)
   if(.not. allocated(wfpDT%ELi)) allocate(wfpDT%fi(nParam),wfpDT%fij(nParam,nParam),wfpDT%ELi(nParam))
   WFP => ElocAndPsiTermsLM_getWFP(EPsiTLM)
   sSize=getSampleSize(sample)
   rSize=getSampleSizeAllNodes(sample)
   call getSampleSizes(sample, sampleSizes)
   if (rSize<=nParam) call abortp('varmin2(NL2SOL):You need at least&
      & nParams+1 SampleSize*CPUs to use NL2SOL. Icrease sample size or CPUs')
   call assert(nParam>0,'varmin2_optimizeSample: no parameters')
   allocate(x(nParam),v(93+rSize*nParam+3*rSize+nParam*(3*nParam+33)/2),iv(60+nParam))   !v is needed for nl2sol, jcb is Jacobian Matrix, iv is control variable for nl2sol

   ! iteration of individual nl2sol call with newly equilibrated samples


      call dfault ( iv, v ) ! defaults for nl2sol
      iv(14)=0
      iv(15)=0
      iv(19)=0
      iv(20)=0
      iv(21)=-1
      iv(22)=0
      iv(23)=0
      iv(24)=0
      iv(17)=400   ! Max function calls
      iv(16)=0
      if(iv16Present)  iv(16)=iv16

      if(nliterPresent) iv(18)=max_iter
      d1 = 94 + 2*rSize + ( nParam * ( 3 * nParam + 31 ) ) / 2
      iv(d) = d1
      r1 = d1 + nParam
      iv(r) = r1
      j1 = r1 + rSize
      iv(j) = j1
      strted = .true.
      sSize=getSampleSize(sample)
      if (logmode >= 2) write(iul,*) ' starting wf parameter optimization with optType=',WFP%optType
      call ElocAndPsiTermsLM_reset(EPsiTLM)
      call internal_calcEPsiTerms()
      eOld = ElocAndPsiTermsLM_EmeanALL(EPsiTLM)

      if(.not.eRefPresent) then
       optERef=eOld
       ! In cases that user do not specify E_ref using sample emean as first E_ref
       ! is not usualy so 0.02 is added to improve initial value.
       if (eqStep==1) optERef=optERef*1.02

     else
       if(.not.fixed .and. eqStep>1) optERef=eOld
     endif
      varOld = ElocAndPsiTermsLM_varALL(EPsiTLM)
      varRefOld = ElocAndPsiTermsLM_varRefALL(EPsiTLM)

       !! add gaussian whit noise for starting param for mo opt
      if(  noiseInit .eqv. .true. .and. npMO > 0 .and. eqStep==1) then
        if(MASTER) then
          x = wfparams_get(WFP)
           do i=npJ+npCI+1,nParam
             x(i)=x(i)+noiseCoeff*mygran()
           enddo
         endif
        call myMPIBcastDouble(x,nParam)
        call wfparams_set(WFP,x)
        call ElocAndPsiTermsLM_reset(EPsiTLM)
        call internal_calcEPsiTerms()
      endif

      if (MASTER .and. logmode>=2) then
         if (logmode >= 2) then
            x = wfparams_get(WFP)
            write(iul,*) 'initial parameters:'
            write(iul,'(10F10.4)') x
            write(iul,'(3(A,F13.5),A,G11.3)') ' initial values: Emean = ',eOld,' var = ',varOld, &
                 ' varRef = ',varRef
         end if
      end if
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NL2SOL begin
     if ( iv(1) /= 0 .and. iv(1) /= 12 ) go to 40! <><><><><><> goto branch

            strted = .false.
            iv(nfcall) = 1
            iv(nfgcal) = 1
     10   continue

         nf = iv(nfcall)

         ncallf=ncallf+1
         call calcR(v(r1),sample,WFP%optType,sSize,rSize,sampleSizes,optERef)
         if ( strted ) then

           if ( nf <= 0 ) then
             iv(toobig) = 1
           end if

           go to 40           ! <><><><><><>goto branch

         end if

         if ( nf <= 0 ) then
           iv(1) = 13
           call itsmry ( v(d1), iv, nParam, v, x )
         end if
     30    continue

         call calcJ(v(j1),nParam,sample,WFP%optType,sSize,rSize,sampleSizes,npJ,npCI,npMO)

         if ( iv(nfgcal) == 0 ) then
           iv(1) = 15
           call itsmry ( v(d1), iv, nParam, v, x )
         end if
         strted = .true.
    40   continue
      if (master) then
         x = wfparams_get(WFP)
         call nl2itr ( v(d1), iv, v(j1), rSize, rSize, nParam, v(r1), v, x )
      end if
          call myMPIBcastInteger(iv(1),1)
          call myMPIBarrier(ierr)
          call myMPIBcastDouble(x,nParam)
          call wfparams_set(WFP,x)
          call ElocAndPsiTermsLM_reset(EPsiTLM)
          call internal_calcEPsiTerms()
     if ( iv(1) == 2 ) then
       go to 30  !<><><><><><><><> goto branch
     end if

     if ( iv(1) < 2 ) then
       go to 10 !<><><><><><><><> goto branch
     end if

         if (MASTER) then   ! NL2SOL messages
          if (logmode >= 2) then
             if ( iv(1) == 3 ) then

              write ( iul, * ) ' '
              write ( iul, * ) 'NL2SOL:    Parameters-convergence.'

              else if ( iv(1) == 4 ) then

              write ( iul, * ) ' '
              write ( iul, * ) 'NL2SOL:    Variance convergence.'

              else if ( iv(1) == 5 ) then

              write ( iul, * ) ' '
              write ( iul, * ) 'NL2SOL:    Parameters and Variance convergence.'

              else if ( iv(1) == 6 ) then

              write ( iul, * ) ' '
              write ( iul, * ) 'NL2SOL:     Absolute function convergence.'

              else if ( iv(1) == 7 ) then

              write ( iul, * ) ' '
              write ( iul, * ) 'NL2SOL:     Singular convergence.'

              else if ( iv(1) == 8 ) then

              write ( iul, * ) ' '
              write ( iul, * ) 'NL2SOL:    False convergence. Check the Jacobian Matrix! '

              else if ( iv(1) == 9 ) then

              write ( iul, * ) ' '
              write ( iul, * ) 'NL2SOL:    Function evaluation limit.'

              else if ( iv(1) == 10 ) then

              write ( iul, * ) ' '
              write ( iul, * ) 'NL2SOL:    Iteration limit.'
           endif
         endif

            x = wfparams_get(WFP)
            call wfparams_set(WFP,x,.true.) ! normalize CI coeffs
            x = wfparams_get(WFP)

            if (logmode >= 2 .and. npCI>0) then
                write ( iul, * ) ' '
                write ( iul, * ) 'Parameters are normalized after optimization: '
                write ( iul, * ) ' '
            endif


            if (logmode >= 2) then
                write ( iul, * ) ' '
                write ( iul, * ) 'Parameters after optimization: '
                write(iul,'(10F10.4)') x
            endif
         endif

         emean = ElocAndPsiTermsLM_EmeanALL(EPsiTLM)
         var = ElocAndPsiTermsLM_varALL(EPsiTLM)
         varRef = ElocAndPsiTermsLM_varRefALL(EPsiTLM)

         !call myMPIBcastDouble(optERef,1)
         if (logmode >= 2) write(iul,'(3(A,F13.5),A,G11.3)') ' Emean(after opt) =',emean,' var = ',var,' varRef = ',varRef

      if (doWritewf) then
         fname = trim(baseName)//'-'//trim(intToStr(optIter))//'.wf'
         call writeWF(fname,.false.,ecp)
      end if

   if (eqStep >= eqIter) exit

      eqStep = eqStep + 1
      ! call "subroutine" subName in .in or macro subName.cmd that should contain
      ! code for equilibrating the sample with the new wave function
      call subloop(subName,sample,exitSubloop)
      if (exitSubloop) exit
   deallocate(x,v,iv)
   end do



   call setCurrentResult(emean,0.d0,var)

   call ElocAndPsiTermsLM_destroy(EPsiTLM)


   contains

      subroutine internal_readInput()
      !-----------------------------!

      integer                           :: iflag
      character(len=3)                  :: s

      eRef=0.0
      fixed=.true.
      call getdbla(lines,nl,'E_ref=',eRef,iflag)
      if (iflag /= 0) then
           eRefPresent = .false.
           fixed = .false.
      endif

      optERef = eRef
      call getstra(lines,nl,'E_ref_',s,iflag)

      if (iflag==0) then
        if (s=='fix') then
          fixed = .true.
        elseif(s=='adp') then
          fixed = .false.
        else
          write(iul,*) " Warning E_ref_",s," is not a defined keyword. continue with fixed E_ref."
         endif
      endif
      if(fixed .and. .not.eRefPresent) &
        call abortp('varmin2: Please specify E_ref if you want to use fixed one.')

      if(fixed) optERef = eRef
      max_iter = 0
      call getinta(lines,nl,'max_iter=',max_iter,iflag)
      if (iflag == 0) nliterPresent=.true.
      call getinta(lines,nl,'NL2SOL_D_mode=',iv16,iflag)
      if (iflag == 0) iv16Present=.true.
      doWriteWF = finda(lines,nl,'write_wf')
      subName = 'equilibrate'
      call getstra(lines,nl,'eq_call=',subName,iflag)
      eqIter = 0
      call getinta(lines,nl,'eq_iter=',eqIter,iflag)

      call getdbla(lines,nl,'mo_noise_coeff=',noiseCoeff,iflag)
      if (iflag == 0) then
           noiseInit = .true.
      endif

      end subroutine internal_readInput


      subroutine internal_calcEPsiTerms()
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for linear method
         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)

         rwp => getFirst(sample)
         do
            call pos(rwp,x,y,z)
            call eConfigArray_set(ec,1,x,y,z)
            call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
            call resetTo_without_Calc(rwp,x,y,z) ! reset rw
            call ElocAndPsiTermsLM_add(EPsiTLM,wfpDT)
         if (.not.isNext(sample)) exit
            rwp => getNext(sample)
         enddo

      end subroutine internal_calcEPsiTerms


   subroutine calcR(r,sample,optType,sSize,rSize,rCnts,optERef)
   !----------------------------------------------------------!
   ! calculates residuals for nl2sol
   implicit none
   type(RWSample), intent(inout)        :: sample
   character(len=9), intent(in)         :: optType
   integer,intent(in)                   :: sSize
   integer,intent(in)                   :: rSize
   integer,intent(in)                   :: rCnts(:)
   real(r8), intent(in)                   :: optERef
   real(r8),intent(out)                   :: r(rSize)
   integer                              :: ierr
   real(r8),allocatable                   :: elSend(:),elRec(:)
     allocate(elSend(sSize), elRec(rSize), stat=ierr)
     call assert(ierr == 0, "allocation failed in calcR")
     call calcElocArr(sample,optType,sSize,elSend)
     elSend(:)=elSend(:)-optERef
     call myMPIGatherDoubleV(elSend,sSize,elRec,rCnts,ierr)
     r=elRec
     deallocate(elSend, elRec)
    return
   end subroutine calcR


   subroutine calcJ(j,nppp,sample,optType,sSize,rSize,rCnts,npJ,npCI,npMO)
   !----------------------------------------------------------------!
   ! calculates Jacobian matrix
   ! caller is responsible for allocation of j array
    type(RWSample), intent(inout)  :: sample
    character(len=9), intent(in)   :: optType
    integer,intent(in)             :: sSize, rSize
    integer,intent(in)             :: rCnts(:)
    integer, intent(in)            :: npJ,npCI,npMO ! # of Jastrow, CI and MO parameters
    integer,intent(in)             :: nppp
    real(r8),intent(inout)           :: j(rSize,nppp)
    real(r8), allocatable            :: blockSend(:,:)
    real(r8), allocatable            :: tblockSend(:,:)
    real(r8), allocatable            :: blockRec(:,:)
    real(r8), allocatable            :: tblockRec(:)
    integer                        :: ierr
    integer                        :: pCnts(size(rCnts))

      allocate(blockSend(sSize,nppp),tblockSend(nppp,sSize), blockRec(nppp,rSize), &
        tblockRec(nppp*rSize), stat=ierr)

      call assert(ierr == 0, "allocation in calcJ failed")

      pCnts(:) = nppp * rCnts(:)

        call calcJblock(sample,optType,nppp,sSize,npJ,npCI,npMO,blockSend)

        tblockSend=transpose(blockSend)
      call  myMPIGatherDoubleV(tblockSend,nppp*sSize,tblockRec,pCnts,ierr)
      if(MASTER) then
        blockRec=reshape(tblockRec, (/nppp,rSize/))
        j=transpose(blockRec)
      endif

      deallocate(blockSend, blockRec, tblockRec,tblockSend)
    end subroutine calcJ


    subroutine calcElocArr(sample,optType,sSize,elA)
    !----------------------------------------------!
    !Calculates local energy for each sample and return it as an array.
       type(RWSample), intent(inout)        :: sample
       character(len=9), intent(in)         :: optType
       integer,intent(in)                   :: sSize
       real(r8),intent(inout)                 :: elA(sSize)
       real(r8)                               :: x(ne),y(ne),z(ne)
       type(RandomWalker), pointer          :: rwp
       integer                              :: i
       type(eConfigArray)                   :: ec

       call eConfigArray_new(ec,ne,1)
        rwp => getFirst(sample)
          i=0
         do
          i=i+1
            call pos(rwp,x,y,z)
            call eConfigArray_set(ec,1,x,y,z)
            call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
            elA(i)=elEloc(1)
          if (.not.isNext(sample)) exit
            rwp => getNext(sample)
         enddo
       return
   end subroutine calcElocArr

   subroutine calcJblock(sample,optType,nppp,sSize,npJ,npCI,npMO,blk)
   !-----------------------------------------------------------!
   type(RWSample), intent(inout)  :: sample
   character(len=9), intent(in)   :: optType
   integer,intent(in)             :: nppp
   integer,intent(in)             :: sSize
   integer, intent(in)            :: npJ,npCI,npMO     ! # of Jastrow and CI parameters
   real(r8),intent(inout)           :: blk(sSize,nppp)
   integer                        :: i
   type(RandomWalker), pointer    :: rwp
   real(r8)                         :: x(ne),y(ne),z(ne)
   type(eConfigArray)  :: ec
   i=0

    call eConfigArray_new(ec,ne,1)

      rwp => getFirst(sample)
      do
         i=i+1
         call pos(rwp,x,y,z)
         !call resetTo(rwp,x,y,z,optType)  ! recalculate and reset walker
         call eConfigArray_set(ec,1,x,y,z)
         call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)

            blk(i,:)=wfpDT%Eli(:)

      if (.not.isNext(sample)) exit
         rwp => getNext(sample)
      enddo
   end subroutine calcJblock

   end subroutine varmin2_optimizeSample

end module optParamsVNL2SOL_m
