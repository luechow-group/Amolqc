! Copyright (C) 1996, 2013-2015, 2018 Arne Luechow
! Copyright (C) 2014-2017 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optParamsVarmin_m

   ! implements Levenberg-Marquardt variance minimization of a fixed sample
   ! without weighting

   use kinds_m, only: r8
   use global_m
   use subloop_m, only: subloop
   use rWSample_m
   use elocAndPsiTermsLM_m
   use wfParameters_m
   use parsing_m
   use waveFunction_m, only: writeWF
   use ecp_m, only: EcpType
   use utils_m, only: intToStr
   use mpiInterface_m, only: myMPIBcastDouble

   implicit none

   private
   public :: varmin1_optimizeSample


contains


   subroutine varmin1_optimizeSample(lines,nl,WFP,sample,converged)
   !--------------------------------------------------------------!

   ! optimize the variance for a fixed sample w.r.t. the reference
   ! energy (optERef) using a Levenberg-Marquardt-type algorithm
   integer, intent(in)                  :: nl
   character(len=120), intent(in)       :: lines(nl)
   type(WFParamDef), pointer              :: WFP
   type(RWSample), intent(inout)        :: sample  ! (fixed) sample for optimization
   logical, intent(out)                 :: converged
   integer                              :: np
   real(r8), allocatable                  :: p(:)       ! parameter vector
   real(r8), allocatable                  :: delta_p(:) ! change of parameter vector
   real(r8), allocatable                  :: g(:),b(:),H(:,:)  ! gradient and Hessian
   real(r8) emean, var, eOld, varOld, lambda,varRef,varRefOld,optERef,optLMlambda
   real(r8)                               :: eRef
   integer lwork,iter,i,j,info,eqIter,eqStep,lmIter,optIter
   logical eRefPresent, fixed
   integer, allocatable                 :: ipiv(:)
   real(r8), allocatable                  :: work(:)
   character(len=80)                    :: subName,fname
   logical                              :: doWriteWF
   type(ElocAndPsiTermsLM)              :: EPsiTLM
   type(EcpType)                         :: ecp
   logical                              :: exitSubloop = .false.
   eRefPresent = .true.
   converged = .true.  !! not yet implemented

   call internal_readInput()

   if (logmode>=2) then
      write(iul,'(/A/A)') '   - -  Levenberg-Marquardt variance minimization  - -',' with parameters:'
      if (fixed .and. .not.eRefPresent) then
         write(iul,'(1X,A,2(A,I4),A,G11.3)') &
           ' E_ref = Emean',' max_iter = ',eqIter,' micro_iter = ',lmIter,' lambda = ',optLMlambda
      else
         write(iul,'(1X,A,F15.5,2(A,I4),A,G11.3)') &
          ' E_ref = ',optERef,' max_iter = ',eqIter,' micro_iter = ',lmIter,' lambda = ',optLMlambda
      endif
   endif

   call ElocAndPsiTermsLM_create(EPsiTLM,eRef,WFP)

   np = ElocAndPsiTermsLM_nParams(EPsiTLM)
   WFP => ElocAndPsiTermsLM_getWFP(EPsiTLM)
   call assert(np>0,'varmin1_optimizeSample: no parameters')
   allocate(p(np),delta_p(np),b(np),g(np),H(np,np),ipiv(np),work(np*np))
   eqStep = 1
   do
    if(doWriteWF) call getPlusOptIter(optIter)
     p = 0; delta_p = 0; b = 0; g = 0; H = 0
     lambda = optLMlambda
     eOld = 0; varOld = 0; varRefOld = 0; emean = 0; var = 0; varRef = 0

     if (logmode >= 2) write(iul,*) ' starting wf parameter optimization with optType=',WFP%optType



     call ElocAndPsiTermsLM_reset(EPsiTLM)
     call internal_calcEPsiTerms()
     eOld = ElocAndPsiTermsLM_EmeanALL(EPsiTLM)
     if (.not.eRefPresent) then
         optERef=eOld
         ! In cases that user do not specify E_ref using sample emean as first E_ref
         ! is not usualy so 0.02 is added to improve initial value.
         if (eqStep==0) optERef=optERef*1.02

     else
      if (.not.fixed .and. eqStep>0) optERef=eOld
     end if
     varOld = ElocAndPsiTermsLM_varALL(EPsiTLM)
     varRefOld = ElocAndPsiTermsLM_varRefALL(EPsiTLM)

     if (MASTER .and. logmode>=2) then
        p = wfparams_get(WFP)
        if (logmode >= 2) then
           write(iul,*) 'initial parameters:'
           write(iul,'(10F10.4)') p
           write(iul,'(3(A,F13.5),A,G11.3)') ' initial values: Emean= ',eOld,' var= ',varOld, &
                ' varRef = ',varRef,' lambda= ',lambda
        end if
     end if

     do iter=1,lmIter

        g = ElocAndPsiTermsLM_VxALL(EPsiTLM)
        H = ElocAndPsiTermsLM_VxyALL(EPsiTLM)

        if (logmode >= 3) then
           write(iul,*) ' gradient of variance:'
           write(iul,'(15G10.3)') g
           write(iul,*) ' approx. Hessian of variance:'
           do i=1,np
              write(iul,'(15G10.3)') H(i,:)
           enddo
        endif

        if (MASTER) then
           do i=1,np
              H(i,i) = H(i,i)*(1+lambda)
           enddo
           p = wfparams_get(WFP)
           lwork = np*np
           b = -g
           call DSYSV('L',np,1,H,np,ipiv,b,np,work,lwork,info)
           if (info/=0 .and. logmode>=2) write(iul,*) 'DSYSV failed: INFO=',info

           if (logmode >= 3) then
              write(iul,*) ' DSYSV result b:'
              write(iul,'(10F10.4)') b
           endif

           delta_p = b
           p = p + delta_p

           if (logmode >= 2) then
              write(iul,'(a,i5)') ' new parameter vector in iteration ',iter
              write(iul,'(10F10.4)') p
           endif
        endif

        call myMPIBcastDouble(p,np)
        call wfparams_set(WFP,p)

        call ElocAndPsiTermsLM_reset(EPsiTLM)
        call internal_calcEPsiTerms()
        emean = ElocAndPsiTermsLM_EmeanALL(EPsiTLM)
        var = ElocAndPsiTermsLM_varALL(EPsiTLM)
        varRef = ElocAndPsiTermsLM_varRefALL(EPsiTLM)

        if (varRef < varRefOld) then
           lambda = lambda / 10.d0
        else
           lambda = lambda * 10.d0
        endif
        if (logmode >= 2) write(iul,'(3(A,F13.5),A,G11.3)') ' Emean =',emean,' var = ',var,' varRef = ',varRef, &
           ' new lambda:', lambda
     end do

     call setCurrentResult(emean,0.d0,var)
     call wfparams_set(WFP,p,.true.) ! normalization for ci coeffs
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

   end do

   deallocate(p,delta_p,b,g,H,ipiv,work)

   call ElocAndPsiTermsLM_destroy(EPsiTLM)


   contains


      subroutine internal_readInput()
      !-----------------------------!
      integer iflag
      logical                       :: found
      character(len=3)              :: s

      fixed = .true.
      eRef=0.0
      call getdbla(lines,nl,'E_ref=',eRef,iflag)
      if (iflag /= 0) then
        eRefPresent=.false.
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
          write(iul,*) "!!! Warning E_ref_",s," is not a defined keyword. Continuing with fixed E_ref."
        endif
      endif
      if (fixed .and. .not.eRefPresent) call abortp('lm: Please specify E_ref.')
      if (fixed) optERef = eRef

      lmIter = 3
      call getinta(lines,nl,'max_iter=',lmIter,iflag)
      optLMlambda = 0.001d0
      call getdbla(lines,nl,'lambda=',optLMlambda,iflag)
      doWriteWF = finda(lines,nl,'write_wf')
      subName = 'equilibrate'
      call getstra(lines,nl,'eq_call=',subName,iflag)
      eqIter = 0
      call getinta(lines,nl,'eq_iter=',eqIter,iflag)

      end subroutine internal_readInput

      subroutine internal_calcEPsiTerms()
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for linear method
         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(WFParamDerivTerms) :: wfpDT
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)

         rwp => getFirst(sample)
         do
            call pos(rwp,x,y,z)
            call eConfigArray_set(ec,1,x,y,z)
            call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
            !!!!call resetTo(rwp,x,y,z,WFP%optType)  ! recalculate and reset walker
            call resetTo_without_Calc(rwp,x,y,z) ! reset rw
            call ElocAndPsiTermsLM_add(EPsiTLM,wfpDT)
         if (.not.isNext(sample)) exit
            rwp => getNext(sample)
         enddo

      end subroutine internal_calcEPsiTerms


   end subroutine varmin1_optimizeSample


end module optParamsVarmin_m
