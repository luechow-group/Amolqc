! Copyright (C) 2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module elocAndPsiTermsEBFGS_m
   ! implements ElocAndPsiTerms base class using delegation
   ! for BFGS  energy minimization
   use kinds_m, only: r8
   use global_m
   use wfParameters_m
   use elocAndPsiTermsBase_m
   use mpiInterface_m, only: myMPIReduceSumInteger, myMPIReduceSumDouble
   implicit none

   type ElocAndPsiTermsEBFGS
      private
      type(ElocAndPsiTermsBase) :: EPTB
      real(r8), pointer :: fi(:) => null()     ! notation: Psi -> f, i.e. fi = Psi_i/Psi_0
      real(r8), pointer :: fiEL(:) => null()   ! fiEL = Psi_i/Psi_0*E_L
      real(r8), pointer :: ELi(:) => null()
   end type ElocAndPsiTermsEBFGS

private
public :: ElocAndPsiTermsEBFGS, &
          ElocAndPsiTermsEBFGS_create, ElocAndPsiTermsEBFGS_destroy, ElocAndPsiTermsEBFGS_reset, ElocAndPsiTermsEBFGS_nData, &
          ElocAndPsiTermsEBFGS_add, ElocAndPsiTermsEBFGS_EMean, ElocAndPsiTermsEBFGS_result, ElocAndPsiTermsEBFGS_nParams, &
          ElocAndPsiTermsEBFGS_EMeanALL, ElocAndPsiTermsEBFGS_varALL, ElocAndPsiTermsEBFGS_resultALL, ElocAndPsiTermsEBFGS_getWFP

contains

   !---------------------------------------------------!
   subroutine ElocAndPsiTermsEBFGS_create(this,eRef,wfp)
   !---------------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   real(r8), intent(in)                        :: eRef
   type(WFParamDef), pointer                   :: wfp

   call ElocAndPsiTermsBase_create(this%EPTB,eRef,wfp)
   allocate(this%fi(wfp%nParams),this%ELi(wfp%nParams),this%fiEL(wfp%nParams))
   this%fi = 0; this%ELi = 0; this%fiEL = 0
   end subroutine ElocAndPsiTermsEBFGS_create

   !-------------------------------------------!
   subroutine ElocAndPsiTermsEBFGS_destroy(this)
   !-------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   deallocate(this%fi,this%ELi,this%fiEL)
   end subroutine ElocAndPsiTermsEBFGS_destroy

   !-----------------------------------------!
   subroutine ElocAndPsiTermsEBFGS_reset(this)
   !-----------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   call ElocAndPsiTermsBase_reset(this%EPTB)
   this%fi = 0; this%Eli = 0; this%fiEL = 0
   end subroutine ElocAndPsiTermsEBFGS_reset

   !-----------------------------------------------!
   integer function ElocAndPsiTermsEBFGS_nData(this)
   !-----------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   ElocAndPsiTermsEBFGS_nData = ElocAndPsiTermsBase_nData(this%EPTB)
   end function ElocAndPsiTermsEBFGS_nData

   !--------------------------------------------------!
   integer function ElocAndPsiTermsEBFGS_nDataALL(this)
   !--------------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   ElocAndPsiTermsEBFGS_nDataALL = ElocAndPsiTermsBase_nDataALL(this%EPTB)
   end function ElocAndPsiTermsEBFGS_nDataALL

   !-------------------------------------------------!
   integer function ElocAndPsiTermsEBFGS_nParams(this)
   !-------------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   ElocAndPsiTermsEBFGS_nParams = ElocAndPsiTermsBase_nParams(this%EPTB)
   end function ElocAndPsiTermsEBFGS_nParams

   !------------------------------------------------------!
   function ElocAndPsiTermsEBFGS_getWFP(this) result(wfp_p)
   !------------------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   type(WFParamDef), pointer :: wfp_p
   wfp_p => ElocAndPsiTermsBase_getWFP(this%EPTB)
   end function ElocAndPsiTermsEBFGS_getWFP


   subroutine ElocAndPsiTermsEBFGS_add(this,wfpDT)
   !---------------------------------------------!
      type(ElocAndPsiTermsEBFGS), intent(inout) :: this
      type(WFParamDerivTerms), intent(in)     :: wfpDT   
      integer np

      call ElocAndPsiTermsBase_add(this%EPTB,wfpDT)
      np = this%EPTB%wfp%nParams

      if (np /= size(wfpDT%fi)) call abortp("ElocAndPsiTermsENR_add: inconsistent sizes")

      this%fi(:) = this%fi(:) + wfpDT%fi(:)
      this%ELi(:) = this%ELi(:) + wfpDT%ELi(:) 
      this%fiEL(:) = this%fiEL(:) + wfpDT%fi(:) * wfpDT%eloc 
   end subroutine ElocAndPsiTermsEBFGS_add


   !----------------------------------------------!
   real(r8) function ElocAndPsiTermsEBFGS_EMean(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   ElocAndPsiTermsEBFGS_EMean = ElocAndPsiTermsBase_Emean(this%EPTB)
   end function ElocAndPsiTermsEBFGS_EMean

   !-------------------------------------------------!
   real(r8) function ElocAndPsiTermsEBFGS_EMeanALL(this)
   !-------------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   ElocAndPsiTermsEBFGS_EMeanALL = ElocAndPsiTermsBase_EmeanALL(this%EPTB)
   end function ElocAndPsiTermsEBFGS_EMeanALL

   !--------------------------------------------!
   real(r8) function ElocAndPsiTermsEBFGS_var(this)
   !--------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   ElocAndPsiTermsEBFGS_var = ElocAndPsiTermsBase_var(this%EPTB)
   end function ElocAndPsiTermsEBFGS_var

   !-----------------------------------------------!
   real(r8) function ElocAndPsiTermsEBFGS_varALL(this)
   !-----------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   ElocAndPsiTermsEBFGS_varALL = ElocAndPsiTermsBase_varALL(this%EPTB)
   end function ElocAndPsiTermsEBFGS_varALL

   !---------------------------------------------------------!
   subroutine ElocAndPsiTermsEBFGS_resultALL(this,fi,ELi,fiEL)
   !---------------------------------------------------------!
   !
   ! note: copying of array unnecessary but required for reshape;
   ! better solution: how to transfer arrays (contiguous memory) with MPI?
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   real(r8), intent(out) :: fi(:),ELi(:),fiEL(:)  ! for notation see type def.
   integer nData(1),nTotalData(1),np,np2,n
   real(r8) vsum(this%EPTB%wfp%nParams)
   np = this%EPTB%wfp%nParams
   np2 = np**2
   nData(1) = this%EPTB%nData
   call myMPIReduceSumInteger(nData,nTotalData,1)
   if (mytid==0) then
      n = nTotalData(1)
   else
      n = 1
   endif

   ! reduce and calc vectors
   vsum = 0
   call myMPIReduceSumDouble(this%fi,vsum,np)
   fi = vsum / n
   call myMPIReduceSumDouble(this%ELi,vsum,np)
   ELi = vsum / n
   call myMPIReduceSumDouble(this%fiEL,vsum,np)
   fiEL = vsum / n

   end subroutine ElocAndPsiTermsEBFGS_resultALL

   !------------------------------------------------------!
   subroutine ElocAndPsiTermsEBFGS_result(this,fi,ELi,fiEL)
   !------------------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(in) :: this
   real(r8), intent(out) :: fi(:),ELi(:),fiEL(:)  ! for notation see type def.
   integer n,np
   n = this%EPTB%nData
   np = size(fi)

   fi = this%fi / n
   ELi = this%ELi / n
   fiEL = this%fiEL / n

   end subroutine ElocAndPsiTermsEBFGS_result

   !------------------------------------------------------!
   subroutine ElocAndPsiTermsEBFGS_gradALL(this,e0,fi,fiEL)
   !------------------------------------------------------!
   !
   ! energy and terms for gradient only
   type(ElocAndPsiTermsEBFGS), intent(inout) :: this
   real(r8), intent(out) :: e0
   real(r8), intent(out) :: fi(:),fiEL(:)  ! for notation see type def.
   integer nData(1),nTotalData(1),np,n
   real(r8) vsum(this%EPTB%wfp%nParams)
   np = this%EPTB%wfp%nParams
   nData(1) = this%EPTB%nData

   e0 = ElocAndPsiTermsBase_EmeanALL(this%EPTB)

   call myMPIReduceSumInteger(nData,nTotalData,1)
   if (mytid==0) then
      n = nTotalData(1)
   else
      n = 1
   endif

   ! reduce and calc vectors
   vsum = 0
   call myMPIReduceSumDouble(this%fi,vsum,np)
   fi = vsum / n
   call myMPIReduceSumDouble(this%fiEL,vsum,np)
   fiEL = vsum / n

   end subroutine ElocAndPsiTermsEBFGS_gradALL

   !---------------------------------------------------!
   subroutine ElocAndPsiTermsEBFGS_grad(this,e0,fi,fiEL)
   !---------------------------------------------------!
   type(ElocAndPsiTermsEBFGS), intent(in) :: this
   real(r8), intent(out) :: e0
   real(r8), intent(out) :: fi(:),fiEL(:)  ! for notation see type def.
   integer n,np
   n = this%EPTB%nData
   np = size(fi)

   e0 = ElocAndPsiTermsBase_EMean(this%EPTB)
   fi = this%fi / n
   fiEL = this%fiEL / n

   end subroutine ElocAndPsiTermsEBFGS_grad

end module elocAndPsiTermsEBFGS_m
