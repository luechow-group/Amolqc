! Copyright (C) 2014-2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module elocAndPsiTermsWEBFGS_m
   ! implements ElocAndPsiTerms base class using delegation
   ! for BFGS  energy minimization
   use kinds_m, only: r8
   use global_m
   use wfParameters_m
   use elocAndPsiTermsBase_m
   use mpiInterface_m, only: myMPIReduceSumInteger, myMPIReduceSumDouble
   implicit none

   type ElocAndPsiTermsWEBFGS
      private
      type(ElocAndPsiTermsBase) :: EPTB
      real(r8)          :: w = 0              ! weights (Psi/Psi0)**2
      real(r8)          :: wEL = 0            ! weighted Eloc
      real(r8), pointer :: wk(:) => null()    ! deriv of w wrt p_k
      real(r8), pointer :: wkEL(:) => null()  ! deriv of w wrt p_k x Eloc
      real(r8), pointer :: wELk(:) => null()  ! w x deriv of Eloc
   end type ElocAndPsiTermsWEBFGS

private
public :: ElocAndPsiTermsWEBFGS, &
          ElocAndPsiTermsWEBFGS_create, ElocAndPsiTermsWEBFGS_destroy, ElocAndPsiTermsWEBFGS_reset, ElocAndPsiTermsWEBFGS_nData, &
          ElocAndPsiTermsWEBFGS_add, ElocAndPsiTermsWEBFGS_EMean, ElocAndPsiTermsWEBFGS_result, ElocAndPsiTermsWEBFGS_nParams, &
          ElocAndPsiTermsWEBFGS_EMeanALL, ElocAndPsiTermsWEBFGS_varALL, ElocAndPsiTermsWEBFGS_resultALL, ElocAndPsiTWEBFGS_getWFP

contains

   !----------------------------------------------------!
   subroutine ElocAndPsiTermsWEBFGS_create(this,eRef,wfp)
   !----------------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   real(r8), intent(in)                        :: eRef
   type(WFParamDef), pointer                   :: wfp

   call abortp("ElocAndPsiTermsWEBFGS: needs check of maths in derivs (_add): Jastrow only?")

   call ElocAndPsiTermsBase_create(this%EPTB,eRef,wfp)
   allocate(this%wk(wfp%nParams),this%wkEL(wfp%nParams),this%wELk(wfp%nParams))
   this%wk = 0; this%wkEL = 0; this%wELk = 0
   end subroutine ElocAndPsiTermsWEBFGS_create

   !--------------------------------------------!
   subroutine ElocAndPsiTermsWEBFGS_destroy(this)
   !--------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   deallocate(this%wk,this%wkEL,this%wELk)
   end subroutine ElocAndPsiTermsWEBFGS_destroy

   !------------------------------------------!
   subroutine ElocAndPsiTermsWEBFGS_reset(this)
   !------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   call ElocAndPsiTermsBase_reset(this%EPTB)
   this%w = 0; this%wEL = 0
   this%wk = 0; this%wkEL = 0; this%wELk = 0
   end subroutine ElocAndPsiTermsWEBFGS_reset

   !------------------------------------------------!
   integer function ElocAndPsiTermsWEBFGS_nData(this)
   !------------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   ElocAndPsiTermsWEBFGS_nData = ElocAndPsiTermsBase_nData(this%EPTB)
   end function ElocAndPsiTermsWEBFGS_nData

   !---------------------------------------------------!
   integer function ElocAndPsiTermsWEBFGS_nDataALL(this)
   !---------------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   ElocAndPsiTermsWEBFGS_nDataALL = ElocAndPsiTermsBase_nDataALL(this%EPTB)
   end function ElocAndPsiTermsWEBFGS_nDataALL

   !-------------------------------------------------!
   integer function ElocAndPsiTermsWEBFGS_nParams(this)
   !-------------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   ElocAndPsiTermsWEBFGS_nParams = ElocAndPsiTermsBase_nParams(this%EPTB)
   end function ElocAndPsiTermsWEBFGS_nParams

   !-------------------------------------------------------!
   function ElocAndPsiTWEBFGS_getWFP(this) result(wfp_p)
   !-------------------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   type(WFParamDef), pointer :: wfp_p
   wfp_p => ElocAndPsiTermsBase_getWFP(this%EPTB)
   end function ElocAndPsiTWEBFGS_getWFP


   subroutine ElocAndPsiTermsWEBFGS_add(this,phi0,U0,wfpDT)
   !----------------------------------------------------
      type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
      real(r8), intent(in)                      :: phi0    ! phi(p_0) original parameters
      real(r8), intent(in)                      :: U0      ! U(p_0) original parameters
      type(WFParamDerivTerms), intent(in)     :: wfpDT 
      real(r8) w
      integer np

      call ElocAndPsiTermsBase_add(this%EPTB,wfpDT)
      np = this%EPTB%wfp%nParams

      w = exp(2d0*(wfpDT%U - U0))
      this%w = this%w + w
      this%wEL = this%wEL + w * wfpDT%eloc

      if (np /= size(wfpDT%ELi)) call abortp("ElocAndPsiTermsWEBFGS_add: inconsistent sizes")

      !!!this%wk(:) = this%wk(:) + 2d0 * w * wfpDT%Ui(:)
      !!!this%wkEL(:) = this%wkEL(:) + 2d0 * w * wfpDT%Ui(:) * wfpDT%eloc
      this%wELk(:) = this%wELk(:) + w * wfpDT%Eli(:)
   end subroutine ElocAndPsiTermsWEBFGS_add


   !----------------------------------------------!
   real(r8) function ElocAndPsiTermsWEBFGS_EMean(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   ElocAndPsiTermsWEBFGS_EMean = ElocAndPsiTermsBase_Emean(this%EPTB)
   end function ElocAndPsiTermsWEBFGS_EMean

   !--------------------------------------------------!
   real(r8) function ElocAndPsiTermsWEBFGS_EMeanALL(this)
   !--------------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   ElocAndPsiTermsWEBFGS_EMeanALL = ElocAndPsiTermsBase_EmeanALL(this%EPTB)
   end function ElocAndPsiTermsWEBFGS_EMeanALL

   !---------------------------------------------!
   real(r8) function ElocAndPsiTermsWEBFGS_var(this)
   !---------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   ElocAndPsiTermsWEBFGS_var = ElocAndPsiTermsBase_var(this%EPTB)
   end function ElocAndPsiTermsWEBFGS_var

   !-----------------------------------------------!
   real(r8) function ElocAndPsiTermsWEBFGS_varALL(this)
   !-----------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   ElocAndPsiTermsWEBFGS_varALL = ElocAndPsiTermsBase_varALL(this%EPTB)
   end function ElocAndPsiTermsWEBFGS_varALL

   !-----------------------------------------------------------------!
   subroutine ElocAndPsiTermsWEBFGS_resultALL(this,w,wEL,wk,wkEL,wELk)
   !-----------------------------------------------------------------!
   !
   ! note: copying of array unnecessary but required for reshape;
   ! better solution: how to transfer arrays (contiguous memory) with MPI?
   type(ElocAndPsiTermsWEBFGS), intent(inout) :: this
   real(r8), intent(out) :: w,wEL,wk(:),wkEL(:),wELk(:)  ! for notation see type def.
   integer nData(1),nTotalData(1),np,n
   real(r8) vsum(this%EPTB%wfp%nParams)
   np = this%EPTB%wfp%nParams
   nData(1) = this%EPTB%nData
   call myMPIReduceSumInteger(nData,nTotalData,1)
   if (mytid==0) then
      n = nTotalData(1)
   else
      n = 1
   endif

   ! reduce scalars
   w = 0
   call myMPIReduceSumDouble(this%w,w,1)
   w = w / n
   wEL = 0
   call myMPIReduceSumDouble(this%wEL,wEL,1)
   wEL = wEL / n
   ! reduce and calc vectors
   vsum = 0
   call myMPIReduceSumDouble(this%wk,vsum,np)
   wk = vsum / n
   call myMPIReduceSumDouble(this%wkEL,vsum,np)
   wkEL = vsum / n
   call myMPIReduceSumDouble(this%wELk,vsum,np)
   wELk = vsum / n

   end subroutine ElocAndPsiTermsWEBFGS_resultALL

   !-------------------------------------------------------------!
   subroutine ElocAndPsiTermsWEBFGS_result(this,w,wEL,wk,wkEL,wELk)
   !-------------------------------------------------------------!
   type(ElocAndPsiTermsWEBFGS), intent(in) :: this
   real(r8), intent(out) :: w,wEL,wk(:),wkEL(:),wELk(:)  ! for notation see type def.
   integer n
   n = this%EPTB%nData

   w   = this%w / n
   wEL = this%wEL / n
   wk  = this%wk / n
   wkEL= this%wkEL / n
   wELk= this%wELk / n

   end subroutine ElocAndPsiTermsWEBFGS_result

end module elocAndPsiTermsWEBFGS_m
