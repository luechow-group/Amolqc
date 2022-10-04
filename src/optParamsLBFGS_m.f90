! Copyright (C) 2011, 2015, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optParamsLBFGS_m

   use kinds_m, only: r8
   use global_m
   use rWSample_m
   use elocAndPsiTermsEBFGS_m
   use wfParameters_m
   use qmc_m, only: qmc_init, qmc_run
   use lbfgsb3_m, only: setulb
   use mpiInterface_m, only: myMPIBcastInteger, myMPIBcastDouble

   implicit none

   private
   public :: eminLBFGS_optInit, eminLBFGS_destroy, eminLBFGS_optimizeSample

   integer :: mOptMaxIter = 3
   integer :: mOptMode = 0

   type(ElocAndPsiTermsEBFGS), save :: mEPsiTEBFGS


contains


   !----------------------------------------!
   subroutine eminLBFGS_optInit(lines,nl,wfp)
   !----------------------------------------!

   integer, intent(in)                  :: nl
   character(len=120), intent(in)       :: lines(nl)
   type(WFParamDef), pointer              :: wfp
   integer iflag
   logical                       :: found,finda
   real(r8)                        :: eRef

   mOptMode = 1
   call getinta(lines,nl,'mode=',mOptMode,iflag)
   mOptMaxIter = 5
   call getinta(lines,nl,'max_iter=',mOptMaxIter,iflag)

   eRef = 0.d0
   call ElocAndPsiTermsEBFGS_create(mEPsiTEBFGS,eRef,wfp)

   if (logmode>=2) then
      write(iul,'(/A/A)') '   - - energy minimization using L-BFGS: initialization - -', &
      ' setting parameters:'
      write(iul,*) ' mode = ',mOptMode
   endif

   end subroutine eminLBFGS_optInit


   !----------------------------!
   subroutine eminLBFGS_destroy()
   !----------------------------!

   call ElocAndPsiTermsEBFGS_destroy(mEPsiTEBFGS)
   end subroutine eminLBFGS_destroy


   !-----------------------------------------!
   subroutine eminLBFGS_optimizeSample(sample)
   !-----------------------------------------!

   ! optimize the energy for a fixed sample using the L-BFGS alg. (L-BFGSB v2.3)

   type(RWSample), intent(inout)        :: sample  ! (fixed) sample for optimization

   type(WFParamDef), pointer              :: WFP
   integer                              :: np
   integer                              :: res(1)
   real(r8), allocatable                  :: p(:),p1(:),p2(:)       ! parameter vector
   real(r8), allocatable                  :: delta_p(:) ! change of parameter vector
   real(r8), allocatable                  :: g(:)       ! gradient
   real(r8) e0,var
   real(r8), allocatable                  :: fi(:),ELi(:),fiEL(:)

   integer, parameter :: nmax=1024
   integer, parameter :: mmax = 17
   integer, parameter :: lenwa = 2*mmax*nmax + 4*nmax + 11*mmax*mmax + 8*mmax

!     nmax  is the dimension of the largest problem to be solved.
!     mmax  is the maximum number of limited memory corrections.
!     lenwa is the corresponding real workspace required.
   character(len=60)     task, csave
   logical          lsave(4)
   integer          i,n,m,iprint
   integer          nbd(nmax), iwa(3*nmax), isave(44)
   real(r8)           ff, factr, pgtol
   real(r8)           xx(nmax), l(nmax), u(nmax), gg(nmax), dsave(29)
   real(r8)           wa(lenwa)
   np = ElocAndPsiTermsEBFGS_nParams(mEPsiTEBFGS)
   WFP => ElocAndPsiTermsEBFGS_getWFP(mEPsiTEBFGS)
   allocate(p(np),p1(np),p2(np),delta_p(np),g(np))
   allocate(fi(np),ELi(np),fiEL(np))
   e0 = 0; g = 0

   if (logmode >= 3) write(iul,*) 'LBFGS wf parameter optimization  -- optType=',WFP%optType

   p = wfparams_get(WFP)
   call ElocAndPsiTermsEBFGS_reset(mEPsiTEBFGS)
   call calcEPsiTerms(sample,WFP)
   e0 = ElocAndPsiTermsEBFGS_EmeanALL(mEPsiTEBFGS)

   if (logmode>=3) write(iul,*) ' L-BFGS: e0=',e0

   call ElocAndPsiTermsEBFGS_resultALL(mEPsiTEBFGS,fi,ELi,fiEL)
   ! calculate gradient estimator
   if (MASTER) then
      g = 2*( fiEL - e0*fi )
      if (logmode>=3) write(iul,'(15G11.3)') g
   endif

   ! testing: repeat
   call ElocAndPsiTermsEBFGS_reset(mEPsiTEBFGS)
   call calcEPsiTerms(sample,WFP)
   e0 = ElocAndPsiTermsEBFGS_EmeanALL(mEPsiTEBFGS)

   if (logmode>=3) write(iul,*) ' repeat: L-BFGS: e0=',e0

   call ElocAndPsiTermsEBFGS_resultALL(mEPsiTEBFGS,fi,ELi,fiEL)
   if (MASTER) then

      g = 2*( fiEL - e0*fi )

      if (logmode>=3) write(iul,'(15G11.3)') g

      ! L-BFGS-B init
      iprint = -1    ! We suppress the default output
      factr  = 0.0d0 ! own stopping criteria
      pgtol  = 0.0d0
      m      =  5  ! number of limited memory corrections: should be parameter

      ! We now specify nbd which defines the bounds on the variables:
      ! l   specifies the lower bounds, u   specifies the upper bounds.
      do i = 1, np
         nbd(i) = 2
         l(i)   = -100.d0
         u(i)   = +100.d0
      enddo

      task = 'START'
   endif

   ! L-BFGS-B loop
   ! note: only MASTER has e0 and g, and calls bfgs routine
   do
      if (MASTER) then
         call setulb(np,m,p,l,u,nbd,e0,g,factr,pgtol,wa,iwa,task,iprint, &
                     csave,lsave,isave,dsave)

         if (logmode>=3) write(iul,'(2a,g15.5/100g10.3)') 'LBFGS-call:',task,e0,p

         if (task(1:2) == 'FG') then
            res(1)=1
         else if (task(1:5) == 'NEW_X') then
            res(1)=2
         else
            res(1)=0
         endif
      endif

      call myMPIBcastInteger(res,1)
      if (res(1)==1) then ! task(1:2) == 'FG'
         call myMPIBcastDouble(p,np)
         call wfparams_set(WFP,p)

         if (MASTER .and. logmode>=3) write(iul,'(a/100g10.3)') 'new parameters:',p

         call qmc_run(sample)   ! for equilibration necessary, alternatively reweighting
         call ElocAndPsiTermsEBFGS_reset(mEPsiTEBFGS)
         call calcEPsiTerms(sample,WFP)
         e0 = ElocAndPsiTermsEBFGS_EMeanALL(mEPsiTEBFGS)

         if (logmode >= 3) write(iul,*) ' L-BFGS: FG: E=',e0

         call ElocAndPsiTermsEBFGS_resultALL(mEPsiTEBFGS,fi,ELi,fiEL)
         if (MASTER) then
            g = 2*( fiEL - e0*fi )
            if (logmode>=3) write(iul,'(15G11.3)') g
         endif

      else if (res(1)==2) then ! task(1:5)=='NEW_X'
         if (MASTER) then
            if (dsave(13) <= 1.d-10*(1.0d0 + abs(ff))) then
               task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
            else if (isave(34) >= mOptMaxIter) then
               task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
               ! the current iteration number, isave(30),
               ! the total number of f and g evaluations, isave(34),
               ! the value of the objective function f,
               ! the norm of the projected gradient,  dsave(13)
            endif
            write (998,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate', &
                   isave(30),'nfg =',isave(34),'e0 =',e0,'|proj g| =',dsave(13)
            if (task(1:4) == 'STOP') then
               write (iul,*) 'L-BFGS: task=',task
               write (iul,*) 'Final parameters p='
               write (iul,'(10G13.5)') (p(i), i = 1, np)
            endif
         endif

      else
         if (MASTER) then
            if (iprint <= -1 .and. task(1:4) /= 'STOP') write(iul,*) 'L-BFGS err:',task
         endif
         exit
      endif
   enddo

   deallocate(p,p1,p2,g)
   deallocate(fi,ELi,fiEL)

   end subroutine eminLBFGS_optimizeSample


   subroutine calcEPsiTerms(sample,WFP)
   !-----------------------------------
      type(RWSample), intent(inout)     :: sample
      type(WFParamDef), intent(in)      :: WFP
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
         call ElocAndPsiTermsEBFGS_add(mEPsiTEBFGS,wfpDT)
      if (.not.isNext(sample)) exit
         rwp => getNext(sample)
      enddo
   end subroutine calcEPsiTerms


end module optParamsLBFGS_m
