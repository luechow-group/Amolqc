! Copyright (C) 2011, 2015, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optParamsBFGS_m

   use kinds_m, only: r8
   use global_m
   use rwSample_m
   use elocAndPsiTermsENR_m
   use wfParameters_m
   use qmc_m, only : qmc_run
   use mpiInterface_m, only: myMPIBcastDouble, myMPIBcastString
   implicit none

   private
   public :: eminBFGS_optInit, eminBFGS_destroy, eminBFGS_optimizeSample

   integer :: mOptMaxIter = 3
   integer :: mOptMode = 0

   type(ElocAndPsiTermsENR), save :: mEPsiTENR

contains

   !---------------------------------------!
   subroutine eminBFGS_optInit(lines,nl,wfp)
   !---------------------------------------!

   integer, intent(in)                  :: nl
   character(len=120), intent(in)       :: lines(nl)
   type(WFParamDef), pointer              :: wfp
   integer iflag
   real(r8)                        :: eRef

   mOptMode = 1
   call getinta(lines,nl,'mode=',mOptMode,iflag)
   mOptMaxIter = 5
   call getinta(lines,nl,'max_iter=',mOptMaxIter,iflag)

   eRef = 0.d0
   call ElocAndPsiTermsENR_create(mEPsiTENR,eRef,wfp)

   if (logmode>=2) then
      write(iul,'(/A/A)') '   - - energy minimization using BFGS: initialization - -', &
      ' setting parameters:'
      write(iul,*) ' mode = ',mOptMode
   endif

   end subroutine eminBFGS_optInit


   !---------------------------!
   subroutine eminBFGS_destroy()
   !---------------------------!
   call ElocAndPsiTermsENR_destroy(mEPsiTENR)
   end subroutine eminBFGS_destroy


   !----------------------------------------!
   subroutine eminBFGS_optimizeSample(sample)
   !----------------------------------------!

   type(RWSample), intent(inout)        :: sample  ! (fixed) sample for optimization

   type(WFParamDef), pointer              :: WFP
   integer                              :: np
   real(r8), allocatable                  :: p(:),p1(:),p2(:)       ! parameter vector
   real(r8), allocatable                  :: delta_p(:) ! change of parameter vector
   real(r8), allocatable                  :: g(:),g1(:),bb(:),H(:,:)  ! gradient and Hessian
   real(r8) e0
   real(r8), allocatable                  :: fi(:),ELi(:),fiEL(:)
   real(r8), allocatable                  :: fifj(:,:),fifjEL(:,:),fiELj(:,:),fij(:,:),fijEL(:,:)
   real(r8), allocatable                  :: A(:,:),B(:,:),D(:,:)
   real(r8), allocatable                  :: s(:),y(:),bs(:),HB(:,:)
   integer lwork,i,j,info,iter,liter,maxiter
   integer, allocatable                 :: ipiv(:)
   real(r8), allocatable                  :: work(:)
   real(r8) x0,x1,x2,e1,e2,alpha,tmp1,tmp2,xmin,emin

   np = ElocAndPsiTermsENR_nParams(mEPsiTENR)
   WFP => ElocAndPsiTermsENR_getWFP(mEPsiTENR)
   call assert(np>0,'eminNR_optimizeSample: no parameters')
   allocate(p(np),p1(np),p2(np),delta_p(np),bb(np),g(np),g1(np),H(np,np),ipiv(np),work(np*np))
   allocate(fi(np),ELi(np),fiEL(np))
   allocate(fifj(np,np),fifjEL(np,np),fiELj(np,np),fij(np,np),fijEL(np,np))
   allocate(A(np,np),B(np,np),D(np,np))
   allocate(s(np),y(np),bs(np),HB(np,np))

   p = 0; delta_p = 0

   if (logmode >= 3) write(iul,*) ' BFGS wf parameter optimization  -- optType=',WFP%optType

   ! approximate Hessian in BFGS algorithm is initially a unit matrix
   HB = 0.d0
   do i=1,np
      HB(i,i) = 1.d0
   enddo

   call ElocAndPsiTermsENR_reset(mEPsiTENR)
   call calcEPsiTerms(sample,WFP)
   e0 = ElocAndPsiTermsENR_EmeanALL(mEPsiTENR)
   if (logmode>=3) write(iul,*) ' BFGS: e0=',e0
   call ElocAndPsiTermsENR_resultALL(mEPsiTENR,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)

   if (MASTER) then
      g = 2*( fiEL - e0*fi )

      if (logmode>=3) write(iul,'(15G11.3)') g
   endif

   ! repeat as test
   call ElocAndPsiTermsENR_reset(mEPsiTENR)
   call calcEPsiTerms(sample,WFP)
   e0 = ElocAndPsiTermsENR_EmeanALL(mEPsiTENR)
   if (logmode>=3) write(iul,*) ' BFGS: e0=',e0
   call ElocAndPsiTermsENR_resultALL(mEPsiTENR,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)

   if (MASTER) then
      g = 2*( fiEL - e0*fi )

      if (logmode>=3) write(iul,'(15G11.3)') g
   endif

   do iter=1,mOptMaxiter

      if (MASTER) then

         H = calcExactHessian()

         if (logmode >= 3) then
            write(iul,*) ' BFGS: energy = ',e0
            write(iul,*) ' BFGS: gradient:'
            write(iul,'(15G11.3)') g
            write(iul,*) ' BFGS: exact H:'
            do i=1,np
               write(iul,'(15G11.3)') H(i,:)
            enddo
         endif

         lwork = np*np
         bb = -g
         ! NR step
         call DSYSV('L',np,1,H,np,ipiv,bb,np,work,lwork,info)

         if (info/=0 .and. logmode>=2) write(iul,*) 'exact H -- DSYSV failed: INFO=',info
         if (logmode >= 3) then
            write(iul,*) ' DSYSV result b:'
            write(iul,'(10F10.4)') bb
         endif

         bb = -g
         ! NR step with approx Hessian
         H = HB   ! keep HB for updates

         if (logmode >= 3) then
            write(iul,*) ' BFGS: approx H:'
            do i=1,np
               write(iul,'(15G11.3)') H(i,:)
            enddo
         endif

         call DSYSV('L',np,1,H,np,ipiv,bb,np,work,lwork,info)

         if (info/=0 .and. logmode>=2) write(iul,*) 'approx H -- DSYSV failed: INFO=',info
         if (logmode >= 3) then
            write(iul,*) ' DSYSV result b:'
            write(iul,'(10F10.4)') bb
         endif

         alpha = 1.d0 / (2.d0 * maxval(abs(bb)))
         if (logmode >= 3) write(iul,*) 'alpha=',alpha

      endif

      ! line search along direction bb
      maxiter = 2
      do liter = 1,maxiter

         if (MASTER) then
            p = wfparams_get(WFP)
            x0 = 0
            x1 = alpha
            p1 = p + x1*bb

            if (logmode >= 3) write(iul,'(15G11.3)') p1
         endif
         call myMPIBcastDouble(p1,np)
         call wfparams_set(WFP,p1)
         !!! equilibrate sample with new parameters !!!
         call qmc_run(sample)

         call calcEPsiTerms(sample,WFP)
         ! should use grad information for minimization
         e1 = ElocAndPsiTermsENR_EMeanALL(mEPsiTENR)

         if (MASTER) then
            if (logmode >= 3) write(iul,*) liter,'x1=',x1,'e1=',e1

            if (e1 < e0) then
               x2 = 2*alpha
            else
               alpha = alpha / 2
               x2 = alpha
            endif
            p2 = p + x2*bb
         endif
         call myMPIBcastDouble(p2,np)
         call wfparams_set(WFP,p2)
         !!! equilibrate sample with new parameters !!!
         call qmc_run(sample)

         call calcEPsiTerms(sample,WFP)
         e2 = ElocAndPsiTermsENR_EMeanALL(mEPsiTENR)

         if (MASTER) then
            if (logmode >= 3) write(iul,*) liter,'x2=',x2,'e2=',e2

            call parabolaMin(x0,e0,x1,e1,x2,e2,xmin,emin,alpha)

            if (logmode >= 3) write(iul,*) liter,'xmin=',xmin,'emin=',emin,'alpha=',alpha

            p = p + xmin*bb
         endif
         call myMPIBcastDouble(p,np)
         call wfparams_set(WFP,p)
         !!! equilibrate sample with new parameters !!!
         call qmc_run(sample)
         call calcEPsiTerms(sample,WFP)
         e0 = ElocAndPsiTermsENR_EMeanALL(mEPsiTENR)

         if (logmode >= 3) write(iul,*) liter,'xmin=',x1,'e0=',e1
      enddo

      call ElocAndPsiTermsENR_resultALL(mEPsiTENR,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)

      if (MASTER) then
         s = xmin*bb
         g1 = 2*( fiEL - e0*fi )
         y = g1 - g
         bs = matmul(HB,s)
         tmp1 = dot_product(y,s)
         tmp2 = dot_product(s,bs)
         do i=1,np
            HB(i,:) = HB(i,:) + y(i)*y(:)/tmp1 - bs(i)*bs(:)/tmp2
         enddo
         g = g1

         if (logmode >= 3) then
            write(iul,*) ' after line search p='
            write(iul,'(10F10.4)') p
            write(iul,*) ' new gradient g='
            write(iul,'(10F10.4)') g
         endif
      endif

   enddo

   deallocate(p,delta_p,bb,g,g1,H,ipiv,work)
   deallocate(fi,ELi,fiEL)
   deallocate(fifj,fifjEL,fiELj,fij,fijEL)
   deallocate(A,B,D)

   contains

      function calcExactHessian()
         real(r8) :: calcExactHessian(np,np)
         A = 2*( fijEL - fij*e0 - fifjEL + fifj*e0 )
         do j=1,np
            B(:,j) = -2*( fi(:)*g(j) + fi(j)*g(:) )
            D(:,j) = -fi(:)*ELi(j) - fi(j)*ELi(:)
         enddo
         B = B + 4*( fijEL - fij*e0 )
         D = D + fiELj + transpose(fiELj)
         calcExactHessian = A + B + D
      end function


   end subroutine eminBFGS_optimizeSample


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
         call ElocAndPsiTermsENR_add(mEPsiTENR,wfpDT)
      if (.not.isNext(sample)) exit
         rwp => getNext(sample)
      enddo
   end subroutine calcEPsiTerms


   subroutine parabolaMin(x0,y0,x1,y1,x2,y2,xmin,ymin,alpha)
      real(r8), intent(in)    :: x0,y0,x1,y1,x2,y2
      real(r8), intent(inout) :: xmin,ymin,alpha
      integer, parameter :: nmax=3
      integer :: info,n
      integer :: ipiv(nmax)
      real(r8)  :: a(nmax,nmax), b(nmax)
      real(r8) d,d1,d2
      n = nmax
      a(1,:) = (/ 1d0 , x0 , x0**2   /)
      a(2,:) = (/ 1d0 , x1 , x1**2   /)
      a(3,:) = (/ 1d0 , x2 , x2**2   /)
      b(:) = (/ y0, y1, y2 /)
      call dgesv(n,1,a,nmax,ipiv,b,n,info)
      if (info /= 0) then
         write(iul,*) ' parabolaMin: error '
      endif
      xmin = - 0.5d0*b(2)/b(3)
      ymin = (4*b(1)*b(3) - b(2)**2) / (4*b(3))

      d1 = abs(x0 - x1); d2 = abs(x0 - x2)
      if (d1 < d2) then
         d = d2
      else
         d = d1
      endif
      if (abs(x0 - xmin) > max(d1,d2)) then
         xmin = x0 + sign(2*d,xmin-x0)
         ymin = b(1) + b(2)*xmin + b(3)*xmin**2
      endif
      alpha = max(alpha/2d0, abs(x0 - xmin))
   end subroutine


end module optParamsBFGS_m
