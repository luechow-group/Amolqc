! Copyright (C) 2013, 2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module elocAndPsiTermsBase_m
   ! implements ElocAndPsiTerms base class
   use kinds_m, only: r8
   use error_m
   use wfParameters_m
   use mpiInterface_m, only: myMPIAllReduceSumInteger, myMPIReduceSumDouble,&
           myMPIReduceSumInteger, myMPIBcastDouble

   implicit none

   type, public :: ElocAndPsiTermsBase
      real(r8)  :: elocSum  = 0
      real(r8)  :: eloc2Sum = 0
      real(r8)  :: eRef = 0
      integer :: nData = 0
      type(WFParamDef), pointer :: wfp => null()
   end type ElocAndPsiTermsBase

private
public :: ElocAndPsiTermsBase_create, ElocAndPsiTermsBase_reset, ElocAndPsiTermsBase_eRef, &
          ElocAndPsiTermsBase_nData, ELocAndPsiTermsBase_nDataALL, ElocAndPsiTermsBase_nParams, &
          ElocAndPsiTermsBase_getWFP, ElocAndPsiTermsBase_add, ElocAndPsiTermsBase_Emean, &
          ElocAndPsiTermsBase_EmeanALL, ElocAndPsiTermsBase_var, ElocAndPsiTermsBase_varALL, &
          ElocAndPsiTermsBase_varRef, ElocAndPsiTermsBase_varRefALL

contains

   !--------------------------------------------------!
   subroutine ElocAndPsiTermsBase_create(this,eRef,wfp)
   !--------------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   real(r8), intent(in)                       :: eRef     ! reference energy
   type(WFParamDef), pointer                  :: wfp

   call assert(wfp%nParams>0,'ElocAndPsiTermsBase_create: positive number of parameters required')

   this%eRef = eRef
   this%elocSum = 0
   this%eloc2Sum = 0
   this%nData = 0
   this%wfp => wfp
   end subroutine ElocAndPsiTermsBase_create

   !----------------------------------------!
   subroutine ElocAndPsiTermsBase_reset(this)
   !----------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   this%elocSum = 0
   this%eloc2Sum = 0
   this%nData = 0
   end subroutine ElocAndPsiTermsBase_reset

   !--------------------------------------------!
   real(r8) function ElocAndPsiTermsBase_eRef(this)
   !--------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   ElocAndPsiTermsBase_eRef = this%eRef
   end function ElocAndPsiTermsBase_eRef

   !----------------------------------------------!
   integer function ElocAndPsiTermsBase_nData(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   ElocAndPsiTermsBase_nData = this%nData
   end function ElocAndPsiTermsBase_nData

   !-------------------------------------------------!
   integer function ElocAndPsiTermsBase_nDataALL(this)
   !-------------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   integer nTotalData(1),nData(1)
   nData(1) = this%nData
   nTotalData(1) = 0
   call myMPIAllReduceSumInteger(nData,nTotalData,1)
   ElocAndPsiTermsBase_nDataALL = nTotalData(1)
   end function ElocAndPsiTermsBase_nDataALL

   !------------------------------------------------!
   integer function ElocAndPsiTermsBase_nParams(this)
   !------------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   ElocAndPsiTermsBase_nParams = this%wfp%nParams
   end function ElocAndPsiTermsBase_nParams

   !-----------------------------------------------------!
   function ElocAndPsiTermsBase_getWFP(this) result(wfp_p)
   !-----------------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   type(WFParamDef), pointer :: wfp_p
   wfp_p => this%wfp
   end function ElocAndPsiTermsBase_getWFP


   subroutine ElocAndPsiTermsBase_add(this,wfpDT)
   !---------------------------------------------
      type(ElocAndPsiTermsBase), intent(inout) :: this
      type(WFParamDerivTerms), intent(in)      :: wfpDT
      this%elocSum = this%elocSum + wfpDT%eloc
      this%eloc2Sum = this%eloc2Sum + wfpDT%eloc**2
      this%nData = this%nData + 1
   end subroutine ElocAndPsiTermsBase_add


   !---------------------------------------------!
   real(r8) function ElocAndPsiTermsBase_EMean(this)
   !---------------------------------------------!
   type(ElocAndPsiTermsBase), intent(in) :: this
   ElocAndPsiTermsBase_EMean = this%elocSum / this%nData
   end function ElocAndPsiTermsBase_EMean

   !------------------------------------------------!
   real(r8) function ElocAndPsiTermsBase_EMeanALL(this)
   !------------------------------------------------!
   type(ElocAndPsiTermsBase), intent(in) :: this
   integer nData(1),nTotalData(1)
   real(r8) ESum(1), ETotalSum(1)
   ESum(1) = this%elocSum
   ETotalSum(1) = 0
   call myMPIReduceSumDouble(ESum,ETotalSum,1)
   nData(1) = this%nData
   nTotalData(1) = 1
   call myMPIReduceSumInteger(nData,nTotalData,1)
   ETotalSum(1) = ETotalSum(1)/nTotalData(1)
   call myMPIBcastDouble(ETotalSum,1)
   ElocAndPsiTermsBase_EMeanALL = ETotalSum(1)
   end function ElocAndPsiTermsBase_EMeanALL

   !-------------------------------------------!
   real(r8) function ElocAndPsiTermsBase_var(this)
   !-------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   ElocAndPsiTermsBase_var = this%eloc2Sum/this%nData - (this%elocSum/this%nData)**2
   end function ElocAndPsiTermsBase_var

   !----------------------------------------------!
   real(r8) function ElocAndPsiTermsBase_varRef(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   ElocAndPsiTermsBase_varRef = this%eloc2Sum/this%nData    &
      - 2*this%elocSum/this%nData*this%eRef + this%eRef**2
   end function ElocAndPsiTermsBase_varRef

   !----------------------------------------------!
   real(r8) function ElocAndPsiTermsBase_varALL(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   integer nData(1),nTotalData(1)
   real(r8) ESum(1), ETotalSum(1)
   real(r8) E2Sum(1), E2TotalSum(1)
   ESum(1) = this%elocSum
   ETotalSum(1) = 0
   call myMPIReduceSumDouble(ESum,ETotalSum,1)
   E2Sum(1) = this%eloc2Sum
   E2TotalSum(1) = 0
   call myMPIReduceSumDouble(E2Sum,E2TotalSum,1)
   nData(1) = this%nData
   nTotalData(1) = 1
   call myMPIReduceSumInteger(nData,nTotalData,1)
   E2TotalSum(1) = E2TotalSum(1)/nTotalData(1) - (ETotalSum(1)/nTotalData(1))**2
   call myMPIBcastDouble(E2TotalSum,1)
   ElocAndPsiTermsBase_varALL = E2TotalSum(1)
   end function ElocAndPsiTermsBase_varALL

   !-------------------------------------------------!
   real(r8) function ElocAndPsiTermsBase_varRefALL(this)
   !-------------------------------------------------!
   type(ElocAndPsiTermsBase), intent(inout) :: this
   integer nData(1),nTotalData(1)
   real(r8) ESum(1), ETotalSum(1)
   real(r8) E2Sum(1), E2TotalSum(1)
   ESum(1) = this%elocSum
   ETotalSum(1) = 0
   call myMPIReduceSumDouble(ESum,ETotalSum,1)
   E2Sum(1) = this%eloc2Sum
   E2TotalSum(1) = 0
   call myMPIReduceSumDouble(E2Sum,E2TotalSum,1)
   nData(1) = this%nData
   nTotalData(1) = 1
   call myMPIReduceSumInteger(nData,nTotalData,1)
   E2TotalSum(1) = E2TotalSum(1)/nTotalData(1) &
      - 2*ETotalSum(1)/nTotalData(1)*this%eRef + this%eRef**2
   call myMPIBcastDouble(E2TotalSum,1)
   ElocAndPsiTermsBase_varRefALL = E2TotalSum(1)
   end function ElocAndPsiTermsBase_varRefALL

end module elocAndPsiTermsBase_m
