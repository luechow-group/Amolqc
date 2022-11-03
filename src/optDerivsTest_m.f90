! Copyright (C) 2010, 2015, 2018 Arne Luechow
! Copyright (C) 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optDerivsTest_m

use kinds_m, only: r8
use global_m
use rwSample_m
use wfParameters_m
use jastrowParamData_m
use multiDetParam_m
use moParam_m
use eloc_m
use utils_m, only: numDerivative
use random_m, only: myran,init_ran

implicit none

private
public :: optimizeTest

contains

   !--------------------------------------!
   subroutine optimizeTest(lines,nl,sample)
   !--------------------------------------!
   integer, intent(in)            :: nl
   character(len=120), intent(in) :: lines(nl)
   type(RWSample), intent(inout)  :: sample
   integer                        :: points=7  ! points for differentiation rule
   integer                        :: optMode,iflag,verbose,oldLogmode
   character(len=9)               :: optType      ! 'jastrow'|'ci'
   type(WFParamDef)               :: WFP
   real(r8)                         :: h=1.d-5 ! denominator for numerical differentiation
   real(r8)                         :: tol     !tolorence for comparing num and analytical drivatives
   logical                        :: printRep

   call assert(ne>0, ' optimizeTest: wave function must be initialized')
   call assert(getSampleSize(sample) > 0,' optimizeTest: no sample available')

   write(iul,'(/a/)') '!!!!!  optimizeTest requires $ecp(...,no_random_rotation)  !!!!!'
   printRep = .false.
   tol = 1.d-8
   call getdbla(lines,nl,'tol=',tol,iflag)
   if (iflag==0) printRep = .true.
   h=1.d-5
   call getdbla(lines,nl,'h=',h,iflag)
   points = 7
   call getinta(lines,nl,'rule=',points,iflag)

   optType = 'jastrow'
   call getstra(lines,nl,'params=',optType,iflag)
   optMode = 1
   call getinta(lines,nl,'optmode=',optMode,iflag)

   oldLogmode = logmode
   call getinta(lines,nl,'verbose=',verbose,iflag)
   if (iflag==0) logmode = verbose
   if (optType=='mo' .or. optType=='mo+ci' .or. optType=='jas+mo' .or. optType=='jas+mo+ci') then
      call wfparams_init(WFP,optType,optMode,lines,nl)
   else
      call wfparams_init(WFP,optType,optMode)
   endif

   call optimizeTestGen(WFP,sample,h,points,tol,printRep)

   call wfparams_destroy(WFP)

   logmode = oldLogmode

   end subroutine optimizeTest


   subroutine optimizeTestGen(WFP,sample,h,points,tol,printRep)
   !---------------------------------------------!

   type(WFParamDef), intent(in)  :: WFP
   type(RWSample), intent(inout) :: sample
   real(r8), intent(in)            :: h            ! denominator for numerical differentiation
   integer, intent(in)           :: points       ! points for differentiation rule
   real(r8), intent(in)            :: tol          ! tolorence for comparing num and analytical drivatives
   logical, intent(in)           :: printRep

   integer                       :: np,pnt,i
   real(r8)                        :: x(ne),y(ne),z(ne),u0, EL0, phi0
   real(r8)                        :: ELknum, Uk, phik
   real(r8), allocatable           :: p(:),p0(:),f(:),g(:),uk0(:),ukgrad0(:,:),uklapl0(:)
   real(r8), allocatable           :: uknum(:),ukgradnum(:,:),uklaplnum(:)
   real(r8), allocatable           :: ELk(:),fi0(:)
   character(len=9)              :: optType
   type(RandomWalker), pointer   :: rwp
   type(eConfigArray)            :: ec
   type(WFParamDerivTerms)       :: wfpDT
   logical                       :: Psierr, ELkerr

   Psierr = .false.
   ELkerr = .false.
   np = WFP%nParams
   optType = WFP%optType
   allocate(p(np),p0(np),f(-points/2:points/2),g(-points/2:points/2))
   allocate(uk0(np),ukgrad0(3*ne,np),uklapl0(np))
   allocate(uknum(np),ukgradnum(3*ne,np),uklaplnum(np))
   allocate(ELk(np),fi0(np))

   rwp => getFirst(sample)
   call pos(rwp,x,y,z)
   call eConfigArray_new(ec,ne,1)
   call eConfigArray_set(ec,1,x,y,z)
   call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)

   u0 = elU(1)
   phi0 = elPhi(1)
   EL0 = elEloc(1)
   fi0 = wfpDT%fi
   ELk = wfpDT%ELi

!    ulapl0 = elUlapl(1)
!    uk0 = uk
!    ukgrad0 = ukgrad
!    uklapl0 = uklapl

   p = wfparams_get(WFP)
   p0 = p
   write(iul,*) ' * * *  parameter derivative test  * * * '
   write(iul,*)
   write(iul,'(a,i3,a,g20.10)') 'rule=',points,' h=',h
   write(iul,*) 'testing Psi_k/Psi_0:'
   write(iul,*) ' (p, numerical, computed, abs. error)'
   do i=1,np
      f(0) = phi0
      g(0) = u0
      if (logmode>2) write(iul,'(i3,2g20.12)') 0,f(0),g(0)
      do pnt=1,points/2
         p(i) = p0(i) + pnt*h
         call wfparams_set(WFP,p)
         call eloc(0,ec,optType)
         f(pnt) = elPhi(1)
         g(pnt) = elU(1)
         if (logmode>2) write(iul,'(i3,2g20.12)') pnt,f(pnt),g(pnt)
         p(i) = p0(i) - pnt*h
         call wfparams_set(WFP,p)
         call eloc(0,ec,optType)
         f(-pnt) = elPhi(1)
         g(-pnt) = elU(1)
         if (logmode>2) write(iul,'(i3,2g20.12)') pnt,f(-pnt),g(-pnt)
      enddo
      phik = numDerivative(f,h)
      Uk   = numDerivative(g,h)

      write(iul,*) i,phik/phi0 + Uk, fi0(i), phik/phi0 + Uk - fi0(i)
      if (phik/phi0 + Uk - fi0(i)>tol) psiErr = .true.
!!      call assertEqualRelative(uk0(i),uknum(i),1.d-5, &
!!                               '(optimizeTestJastrow): Computed vs Numerical U_k failed')
!!      call assertEqualRelative(uklapl0(i),uklaplnum(i),1.d-5,  &
!!                               '(optimizeTestJastrow): Computed vs Numerical Lapl U_k failed')
      p = p0
   enddo

   ! test of EL_k
   write(iul,*) ' * * *  parameter derivative test  * * * '
   write(iul,*)
   write(iul,*) 'testing EL_k:'
   write(iul,*) ' (p, numerical, computed, abs. error)'
   p = p0
   call wfparams_set(WFP,p)
   do i = 1, np
      f(0) = EL0
      do pnt = 1, points/2
         p(i) = p0(i) + pnt*h
         call wfparams_set(WFP,p)
         call eloc(0,ec,optType)
         f(pnt) = elEloc(1)
         p(i) = p0(i) - pnt*h
         call wfparams_set(WFP,p)
         call eloc(0,ec,optType)
         f(-pnt) = elEloc(1)
      enddo
      ELknum = numDerivative(f,h)
      write(iul,*) i,ELknum,ELk(i),(ELknum-ELk(i))
      if (abs(ELknum-ELk(i))> tol) elkErr = .true.
      !!!call assertEqualRelative(ELk(i),ELknum,1.d-5, &
      !!!                         '(optimizeTestJastrow): Computed vs Numerical EL_k failed')
      p = p0
   enddo
   if (printRep) then
      write(iul,*) ""
      if (.not.elkErr .and. .not.psiErr) then
         write(iul,*) '           all tests passed!'
      else 
         if (psiErr) then
            write(iul,*) '!!!!!!!!!!!!Psi_k/Psi_0 test failed!!!!!!!!!!!!'
         endif
         write(iul,*) ""
         if(elkErr) then
            write(iul,*) '!!!!!!!!!!!!!!!EL_k tests failed!!!!!!!!!!!!!!!'
         endif
      endif
      write(iul,*) ""
   endif
   call eConfigArray_destroy(ec)
   end subroutine optimizeTestGen


end module optDerivsTest_m
