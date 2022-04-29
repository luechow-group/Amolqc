! Copyright (C) 2013, 2015 Kaveh Haghighi Mood
! Copyright (C) 2015, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module elocAndPsiTermsLM_m
   ! implements ElocAndPsiTerms abstract class
   ! for Levenberg Marquardt-type optimization, using delegation from Base class
   use kinds_m, only: r8
   use error_m, only: assert
   use wfParameters_m
   use elocAndPsiTermsBase_m
   use mpiInterface_m, only: myMPIReduceSumDouble, myMPIReduceSumInteger

   implicit none


   type, public :: ElocAndPsiTermsLM
      type(ElocAndPsiTermsBase) :: EPTB
      real(r8), pointer :: vxSum(:) => null()
      real(r8), pointer :: vxySum(:,:) => null()
   end type ElocAndPsiTermsLM

private
public :: ElocAndPsiTermsLM_create, ElocAndPsiTermsLM_destroy, ElocAndPsiTermsLM_reset, &
          ElocAndPsiTermsLM_nData, ElocAndPsiTermsLM_nDataALL, ElocAndPsiTermsLM_nParams, &
          ElocAndPsiTermsLM_add, ElocAndPsiTermsLM_Emean, ElocAndPsiTermsLM_EmeanALL, &
          ElocAndPsiTermsLM_var, ElocAndPsiTermsLM_varALL,ElocAndPsiTermsLM_varRef, &
          ElocAndPsiTermsLM_varRefALL, ElocAndPsiTermsLM_getWFP, &
          ElocAndPsiTermsLM_Vx, ElocAndPsiTermsLM_Vxy, ElocAndPsiTermsLM_VxALL, ElocAndPsiTermsLM_VxyALL, &
          ElocAndPsiTermsLM_nCIParams,ElocAndPsiTermsLM_njParams,ElocAndPsiTermsLM_nMOParams

contains

   !------------------------------------------------!
   subroutine ElocAndPsiTermsLM_create(this,eRef,wfp)
   !------------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   real(r8), intent(in)                     :: eRef     ! reference energy
   type(WFParamDef), pointer                :: wfp

   call ElocAndPsiTermsBase_create(this%EPTB,eRef,wfp)
   allocate(this%vxSum(this%EPTB%wfp%nParams))
   this%vxSum = 0
   allocate(this%vxySum(this%EPTB%wfp%nParams,this%EPTB%wfp%nParams))
   this%vxySum = 0
   end subroutine ElocAndPsiTermsLM_create

   !----------------------------------------!
   subroutine ElocAndPsiTermsLM_destroy(this)
   !----------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   deallocate(this%vxSum)
   deallocate(this%vxySum)
   end subroutine ElocAndPsiTermsLM_destroy

   !--------------------------------------!
   subroutine ElocAndPsiTermsLM_reset(this)
   !--------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   call ElocAndPsiTermsBase_reset(this%EPTB)
   this%vxSum = 0
   this%vxySum = 0
   end subroutine ElocAndPsiTermsLM_reset

   !---------------------------------------------!
   integer function ElocAndPsiTermsLM_nData(this)
   !---------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_nData = ElocAndPsiTermsBase_nData(this%EPTB)
   end function ElocAndPsiTermsLM_nData

   !-----------------------------------------------!
   integer function ElocAndPsiTermsLM_nDataALL(this)
   !-----------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_nDataALL = ElocAndPsiTermsBase_nDataALL(this%EPTB)
   end function ElocAndPsiTermsLM_nDataALL

   !----------------------------------------------!
   integer function ElocAndPsiTermsLM_nParams(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_nParams = ElocAndPsiTermsBase_nParams(this%EPTB)
   end function ElocAndPsiTermsLM_nParams
   !----------------------------------------------!
   integer function ElocAndPsiTermsLM_nCIParams(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_nCIParams = this%EPTB%wfp%nCIParams
   end function ElocAndPsiTermsLM_nCIParams
   !----------------------------------------------!
   integer function ElocAndPsiTermsLM_nJParams(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_nJParams = this%EPTB%wfp%nJParams
   end function ElocAndPsiTermsLM_nJParams
   !----------------------------------------------!
   integer function ElocAndPsiTermsLM_nMOParams(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_nMOParams = this%EPTB%wfp%nMOParams
   end function ElocAndPsiTermsLM_nMOParams

   !---------------------------------------------------!
   function ElocAndPsiTermsLM_getWFP(this) result(wfp_p)
   !---------------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   type(WFParamDef), pointer :: wfp_p
   wfp_p => ElocAndPsiTermsBase_getWFP(this%EPTB)
   end function ElocAndPsiTermsLM_getWFP


   subroutine ElocAndPsiTermsLM_add(this,wfpDT)
   !------------------------------------------!
      type(ElocAndPsiTermsLM), intent(inout) :: this
      type(WFParamDerivTerms), intent(in)     :: wfpDT   
      integer np,j

      call ElocAndPsiTermsBase_add(this%EPTB,wfpDT)
      np = this%EPTB%wfp%nParams

      call assert(np == size(wfpDT%ELi), "ElocAndPsiTermsLM_add: inconsistent sizes")

      this%vxSum(:) = this%vxSum(:) + (wfpDT%eloc - this%EPTB%eRef) * wfpDT%ELi(:) 

      do j = 1, np
         this%vxySum(:,j) = this%vxySum(:,j) + wfpDT%ELi(:) * wfpDT%ELi(j)
      enddo
   end subroutine ElocAndPsiTermsLM_add



   !-------------------------------------------!
   real(r8) function ElocAndPsiTermsLM_EMean(this)
   !-------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_EMean = ElocAndPsiTermsBase_Emean(this%EPTB)
   end function ElocAndPsiTermsLM_EMean

   !----------------------------------------------!
   real(r8) function ElocAndPsiTermsLM_EMeanALL(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_EMeanALL = ElocAndPsiTermsBase_EmeanALL(this%EPTB)
   end function ElocAndPsiTermsLM_EMeanALL

   !-----------------------------------------!
   real(r8) function ElocAndPsiTermsLM_var(this)
   !-----------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_var = ElocAndPsiTermsBase_var(this%EPTB)
   end function ElocAndPsiTermsLM_var

   !--------------------------------------------!
   real(r8) function ElocAndPsiTermsLM_varRef(this)
   !--------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_varRef = ElocAndPsiTermsBase_varRef(this%EPTB)
   end function ElocAndPsiTermsLM_varRef

   !--------------------------------------------!
   real(r8) function ElocAndPsiTermsLM_varALL(this)
   !--------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_varALL = ElocAndPsiTermsBase_varALL(this%EPTB)
   end function ElocAndPsiTermsLM_varALL

   !-----------------------------------------------!
   real(r8) function ElocAndPsiTermsLM_varRefALL(this)
   !-----------------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   ElocAndPsiTermsLM_varRefALL = ElocAndPsiTermsBase_varRefALL(this%EPTB)
   end function ElocAndPsiTermsLM_varRefALL

   !------------------------------------!
   function ElocAndPsiTermsLM_VxALL(this)
   !------------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   real(r8)  :: ElocAndPsiTermsLM_VxALL(this%EPTB%wfp%nParams)
   integer nData(1),nTotalData(1)
   real(r8) VxTotalSum(this%EPTB%wfp%nParams)
   VxTotalSum = 0
   call myMPIReduceSumDouble(this%VxSum,VxTotalSum,this%EPTB%wfp%nParams)
   nData(1) = this%EPTB%nData
   call myMPIReduceSumInteger(nData,nTotalData,1)
   ElocAndPsiTermsLM_VxALL = 2*VxTotalSum/nTotalData(1)
   end function ElocAndPsiTermsLM_VxALL

   !---------------------------------!
   function ElocAndPsiTermsLM_Vx(this)
   !---------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   real(r8) :: ElocAndPsiTermsLM_Vx(this%EPTB%wfp%nParams)
   ElocAndPsiTermsLM_Vx = 2*this%VxSum / this%EPTB%nData
   end function ElocAndPsiTermsLM_Vx


   !-------------------------------------!
   function ElocAndPsiTermsLM_VxyALL(this)
   !-------------------------------------!
   !
   ! note: copying of array unnecessary but required for reshape;
   ! better solution: how to transfer arrays (contiguous memory) with MPI?
   type(ElocAndPsiTermsLM), intent(inout) :: this
   real(r8) :: ElocAndPsiTermsLM_VxyALL(this%EPTB%wfp%nParams,this%EPTB%wfp%nParams)
   integer nData(1),nTotalData(1),np,np2
   real(r8) VxySum(this%EPTB%wfp%nParams,this%EPTB%wfp%nParams),VxyTotalSum(this%EPTB%wfp%nParams,this%EPTB%wfp%nParams)
   real(r8) total(this%EPTB%wfp%nParams**2),vector(this%EPTB%wfp%nParams**2)
   np = this%EPTB%wfp%nParams
   np2 = np**2
   vector = reshape(this%VxySum, (/ np2 /) )
   total = 0
   call myMPIReduceSumDouble(vector,total,np2)
   VxyTotalSum = reshape(total, (/ np,np /) )
   nData(1) = this%EPTB%nData
   call myMPIReduceSumInteger(nData,nTotalData,1)
   ElocAndPsiTermsLM_VxyALL = 2*VxyTotalSum/nTotalData(1)
   end function ElocAndPsiTermsLM_VxyALL

   !----------------------------------!
   function ElocAndPsiTermsLM_Vxy(this)
   !----------------------------------!
   type(ElocAndPsiTermsLM), intent(inout) :: this
   real(r8)                                 :: ElocAndPsiTermsLM_Vxy(size(this%VxySum,1),size(this%VxySum,2))
   ElocAndPsiTermsLM_Vxy = 2*this%VxySum / this%EPTB%nData
   end function ElocAndPsiTermsLM_Vxy

end module elocAndPsiTermsLM_m
