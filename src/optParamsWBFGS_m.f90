! Copyright (C) 2011, 2014-2015, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optParamsWBFGS_m

   use kinds_m, only: r8
   use global_m
   use rWSample_m
   use elocAndPsiTermsWEBFGS_m
   use wfParameters_m
   use lbfgsb3_m, only: setulb
   use mpiInterface_m, only: myMPIBcastDouble, myMPIBcastString
   implicit none

   private
   public :: eminWBFGS_optimizeSample

contains

   subroutine eminWBFGS_optimizeSample(lines,nl,WFP,sample,converged)
   !----------------------------------------------------------------!
   integer, intent(in)                  :: nl
   character(len=120), intent(in)       :: lines(nl)
   type(WFParamDef), pointer              :: WFP
   type(RWSample), intent(inout)        :: sample  ! (fixed) sample for optimization
   logical, intent(out)                 :: converged
   integer                              :: np
   real(r8), allocatable                  :: p(:),p0(:)       ! parameter vector
   real(r8), allocatable                  :: g(:),g1(:)             ! gradient
   real(r8), allocatable                  :: wk(:),wkEL(:),wELk(:)
   !!!real(r8), allocatable                  :: s(:),y(:),bs(:),HB(:,:)
   integer optMode,maxiter,iflag
   real(r8) e0,eps,eRef,var,stdDev
   real(r8) w,wEL
   type(ElocAndPsiTermsWEBFGS)          :: EPsiTWBFGS


   ! LBFGS-B variables
   integer, parameter :: nmax=1024
   integer, parameter :: mmax = 17
   integer, parameter :: lenwa = 2*mmax*nmax + 4*nmax + 11*mmax*mmax + 8*mmax
!     nmax  is the dimension of the largest problem to be solved.
!     mmax  is the maximum number of limited memory corrections.
!     lenwa is the corresponding real workspace required.
   character*60     task, csave
   logical          lsave(4)
   integer          i,m,iprint
   integer          nbd(nmax), iwa(3*nmax), isave(44)
   real(r8)           ff, factr, pgtol
   real(r8)           xx(nmax), l(nmax), u(nmax), gg(nmax), dsave(29)
   real(r8)           wa(lenwa)

   ! LBFGS-B settings
   iprint = -1    ! We suppress the default output.
   factr  = 0.0d0 ! We suppress both code-supplied stopping tests because the
   pgtol  = 0.0d0 ! user is providing his own stopping criteria.
   m      =  5

   optMode = 1
   call getinta(lines,nl,'mode=',optMode,iflag)
   maxIter = 20
   call getinta(lines,nl,'max_iter=',maxIter,iflag)
   eps = 1.d-4

   eRef = 0.d0
   call ElocAndPsiTermsWEBFGS_create(EPsiTWBFGS,eRef,wfp)

   if (logmode>=2) then
      write(iul,'(/A/A)') '   - - weighted energy minimization using BFGS: initialization - -', &
      ' with:'
      write(iul,*) ' mode = ',optMode
   endif

   np = ElocAndPsiTermsWEBFGS_nParams(EPsiTWBFGS)
   call assert(np <= nmax,'eminWBFGS_optimizeSample: too many parameters (change nmax)')
   WFP => ElocAndPsiTWEBFGS_getWFP(EPsiTWBFGS)
   call assert(np>0,'eminWBFGS_optimizeSample: no parameters')
   allocate(p(np),p0(np),g(np),g1(np))
   allocate(wk(np),wkEL(np),wELk(np))
   !!!allocate(s(np),y(np),bs(np))

   p = 0; p0 = 0; g = 0
   wk = 0; wkEL = 0; wELk = 0

   if (logmode >= 2) write(iul,*) ' weighted EBFGS wf parameter optimization  -- optType=',WFP%optType

   if (MASTER) then
      p = wfparams_get(WFP)
      write(iul,*) ' start params, np=',np
      write(iul,'(10f12.4)') p
      p0 = p
   end if

   call ElocAndPsiTermsWEBFGS_reset(EPsiTWBFGS)
   call internal_calcEPsiTerms()
   call ElocAndPsiTermsWEBFGS_resultALL(EPsiTWBFGS,w,wEL,wk,wkEL,wELk)

   if (MASTER) then
      p = wfparams_get(WFP)
      write(iul,*) 'initial params:'
      write(iul,'(10f12.5)') p
      p0 = p
      e0 = wEL/w
      g = wkEL/w - e0*wk/w
      g1 = wkEL/w - e0*wk/w + wELk/w
      write(iul,*) ' start: E_VMC=',e0,' mean abs grad=',sum(abs(g))/np,' max abs grad=',maxval(abs(g))
      write(iul,*) ' start: E_VMC=',e0,' mean abs grad=',sum(abs(g1))/np,' max abs grad=',maxval(abs(g1))
      write(iul,'(10g12.3)') g 
   end if

   ! redo calculation with unchanged parameters
   call myMPIBcastDouble(p,np)
   call wfparams_set(WFP,p)
   call ElocAndPsiTermsWEBFGS_reset(EPsiTWBFGS)
   call internal_calcEPsiTerms()
   call ElocAndPsiTermsWEBFGS_resultALL(EPsiTWBFGS,w,wEL,wk,wkEL,wELk)

   if (MASTER) then
      e0 = wEL/w
      g = wkEL/w - e0*wk/w
      g1 = wkEL/w - e0*wk/w + wELk/w
      write(iul,*) ' start1: E_VMC=',e0,' mean abs grad=',sum(abs(g))/np,' max abs grad=',maxval(abs(g))
      write(iul,*) ' start1: E_VMC=',e0,' mean abs grad=',sum(abs(g1))/np,' max abs grad=',maxval(abs(g1))
      write(iul,'(10g12.3)') g 
   end if

   if (MASTER) then
      p = p + 0.01d0
      write(iul,*) ' new params, np=',np
      write(iul,'(10f12.4)') p
      write(iul,*) ' mean abs delta p=',sum(abs(p-p0))/np,' max abs delta p=',maxval(abs(p-p0))
   end if

   call myMPIBcastDouble(p,np)
   call wfparams_set(WFP,p)
   call ElocAndPsiTermsWEBFGS_reset(EPsiTWBFGS)
   call internal_calcEPsiTerms()
   call ElocAndPsiTermsWEBFGS_resultALL(EPsiTWBFGS,w,wEL,wk,wkEL,wELk)

   if (MASTER) then
      e0 = wEL/w
      g = wkEL/w - e0*wk/w
      g1 = wkEL/w - e0*wk/w + wELk/w
      write(iul,*) ' step1: E_VMC=',e0,' mean abs grad=',sum(abs(g))/np,' max abs grad=',maxval(abs(g))
      write(iul,*) ' step1: E_VMC=',e0,' mean abs grad=',sum(abs(g1))/np,' max abs grad=',maxval(abs(g1))
      write(iul,'(10g12.3)') g 
   end if

   if (MASTER) then
      p = p + 0.01d0
      write(iul,*) ' new params, np=',np
      write(iul,'(10f12.4)') p
      write(iul,*) ' mean abs delta p=',sum(abs(p-p0))/np,' max abs delta p=',maxval(abs(p-p0))
   end if

   call myMPIBcastDouble(p,np)
   call wfparams_set(WFP,p)
   call ElocAndPsiTermsWEBFGS_reset(EPsiTWBFGS)
   call internal_calcEPsiTerms()
   call ElocAndPsiTermsWEBFGS_resultALL(EPsiTWBFGS,w,wEL,wk,wkEL,wELk)

   if (MASTER) then
      e0 = wEL/w
      g = wkEL/w - e0*wk/w
      g1 = wkEL/w - e0*wk/w + wELk/w
      write(iul,*) ' step2: E_VMC=',e0,' mean abs grad=',sum(abs(g))/np,' max abs grad=',maxval(abs(g))
      write(iul,*) ' step2: E_VMC=',e0,' mean abs grad=',sum(abs(g1))/np,' max abs grad=',maxval(abs(g1))
      write(iul,'(10g12.3)') g 
   end if

   ! note the reverse calling mode of the optimizer
   task = 'START'
   do 
      if (MASTER) then     
         call setulb(np,m,p,l,u,nbd,e0,g,factr,pgtol,wa,iwa,task,iprint, &
                  csave,lsave,isave,dsave)
      end if
      call myMPIBcastString(task,60)

      if (task(1:2) == 'FG') then ! request for function f and gradient g at x
         if (logmode>=2) write(iul,*) ' new f and g ...'
         call myMPIBcastDouble(p,np)
         call wfparams_set(WFP,p)
         call ElocAndPsiTermsWEBFGS_reset(EPsiTWBFGS)
         call internal_calcEPsiTerms()
         call ElocAndPsiTermsWEBFGS_resultALL(EPsiTWBFGS,w,wEL,wk,wkEL,wELk)
         if (MASTER) then
            e0 = wEL/w
            g = wkEL/w - e0*wk/w
            !!!g = wkEL/w - e0*wk/w + wELk/w
            if (logmode>=2) then
               write(iul,*) ' step params:'
               write(iul,'(10f12.5)') p
               write(iul,'(a,100f12.5)') 'terms:',w,wEL,wkEL,wELk
               write(iul,*) ' step: E_VMC=',e0,' mean abs grad=',sum(abs(g))/np,' max abs grad=',maxval(abs(g))
            end if
         end if 
      else if (task(1:5) == 'NEW_X') then   
         if (isave(34) >= maxiter) then
            task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
         else if (dsave(13) < eps) then
            ! dsave(13) contains [norm of] projected gradient
            task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
            ! the current iteration number, isave(30)
         end if
         if (logmode>=2) then 
            write (*,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate', &
               isave(30),'nfg =',isave(34),'E0 =',e0,'|proj g| =',dsave(13)
            write(iul,*) ' parameter mean distance: ',sum(abs(p-p0))/np 
         end if
         if (task(1:4) == 'STOP') then
            if (logmode>=2) then
               write (iul,*) task  
               write (iul,*) 'Final p='
               write (iul,'((1x,1p, 6(1x,d11.4)))') (p(i), i = 1, np)
            end if
         end if
      else
         if (iprint <= -1 .and. task(1:4) /= 'STOP' .and. logmode>=2) write(iul,*) task
         exit
      end if
   end do

   ! update sample recalculate with new parameters!
   call recalculateSample(sample)

   call getSampleEnergyAndVarianceAllNodes(sample,e0,var,stdDev)
   write(iul,*) 'sample mean=',e0,' sample var=',var,' sample sigma=',stdDev


   if (logmode>=2) write(iul,*) 'L-BFGS-B: done'

   converged = .false.

   deallocate(p,p0,g)
   !!deallocate(s,y,bs)
   deallocate(wk,wkEL,wELk)

   call ElocAndPsiTermsWEBFGS_destroy(EPsiTWBFGS)

   contains

      subroutine internal_calcEPsiTerms()
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for linear method
         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(WFParamDerivTerms) :: wfpDT
         type(eConfigArray)  :: ec
         real(r8) phi0, U0

         call eConfigArray_new(ec,ne,1)

         rwp => getFirst(sample)
         do
            call pos(rwp,x,y,z)
            call eConfigArray_set(ec,1,x,y,z)
            call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
            phi0 = phi(rwp); U0 = ju(rwp)
            !!!!call resetTo(rwp,x,y,z,WFP%optType)  ! recalculate and reset walker
            call ElocAndPsiTermsWEBFGS_add(EPsiTWBFGS,phi0,U0,wfpDT)
         if (.not.isNext(sample)) exit
            rwp => getNext(sample)
         enddo

      end subroutine internal_calcEPsiTerms

   end subroutine eminWBFGS_optimizeSample


end module optParamsWBFGS_m
