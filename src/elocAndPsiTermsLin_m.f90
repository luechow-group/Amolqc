! Copyright (C) 2013, 2015 Arne Luechow
! Copyright (C) 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module elocAndPsiTermsLin_m
   ! implements ElocAndPsiTerms base class using delegation
   ! for linear energy minimization
   use kinds_m, only: r8
   use global_m
   use wfParameters_m
   use elocAndPsiTermsBase_m
   use mpiInterface_m, only: myMPIReduceSumInteger, myMPIReduceSumDouble
   implicit none

   type ElocAndPsiTermsLin
      private
      type(ElocAndPsiTermsBase) :: EPTB
      real(r8), pointer :: fi(:) => null()     ! notation: Psi -> f, i.e. fi = Psi_i/Psi_0
      real(r8), pointer :: fiEL(:) => null()   ! fiEL = Psi_i/Psi_0*E_L
      real(r8), pointer :: ELi(:) => null()
      real(r8), pointer :: fifj(:,:) => null()    ! fifj = Psi_i/Psi_0 * Psi_j/Psi_0
      real(r8), pointer :: fifjEL(:,:) => null()
      real(r8), pointer :: fiELj(:,:) => null()
   end type ElocAndPsiTermsLin

private
public :: ElocAndPsiTermsLin, ElocAndPsiTermsLin_eRef, &
          ElocAndPsiTermsLin_create, ElocAndPsiTermsLin_destroy, ElocAndPsiTermsLin_reset, ElocAndPsiTermsLin_nData, &
          ElocAndPsiTermsLin_add, ElocAndPsiTermsLin_EMean, ElocAndPsiTermsLin_EMeanALL, &
          ElocAndPsiTermsLin_result, ElocAndPsiTermsLin_resultALL, ElocAndPsiTermsLin_nParams, &
          ElocAndPsiTermsLin_getWFP, ElocAndPsiTermsLin_var, ElocAndPsiTermsLin_varALL, &
          ElocAndPsiTermsLin_nPCI, ElocAndPsiTermsLin_nPMO,  ElocAndPsiTermsLin_nPJ

contains

   !-------------------------------------------------!
   subroutine ElocAndPsiTermsLin_create(this,eRef,wfp)
   !-------------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   real(r8), intent(in)                      :: eRef
   type(WFParamDef), pointer                 :: wfp

   call ElocAndPsiTermsBase_create(this%EPTB,eRef,wfp)
   allocate(this%fi(wfp%nParams),this%ELi(wfp%nParams),this%fiEL(wfp%nParams),this%fifj(wfp%nParams,wfp%nParams))
   allocate(this%fifjEL(wfp%nParams,wfp%nParams),this%fiELj(wfp%nParams,wfp%nParams))
   this%fi = 0; this%ELi = 0; this%fiEL = 0; this%fifj = 0
   this%fifjEL = 0; this%fiELj = 0
   end subroutine ElocAndPsiTermsLin_create

   !----------------------------------------!
   subroutine ElocAndPsiTermsLin_destroy(this)
   !----------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   deallocate(this%fi,this%ELi,this%fiEL,this%fifj)
   deallocate(this%fifjEL,this%fiELj)
   end subroutine ElocAndPsiTermsLin_destroy

   !---------------------------------------!
   subroutine ElocAndPsiTermsLin_reset(this)
   !---------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   call ElocAndPsiTermsBase_reset(this%EPTB)
   this%fi = 0; this%Eli = 0; this%fiEL = 0; this%fifj = 0
   this%fifjEL = 0; this%fiELj = 0
   end subroutine ElocAndPsiTermsLin_reset

   !-------------------------------------------!
   real(r8) function ElocAndPsiTermsLin_eRef(this)
   !-------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_eRef = ElocAndPsiTermsBase_eRef(this%EPTB)
   end function ElocAndPsiTermsLin_eRef

   !---------------------------------------------!
   integer function ElocAndPsiTermsLin_nData(this)
   !---------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_nData = ElocAndPsiTermsBase_nData(this%EPTB)
   end function ElocAndPsiTermsLin_nData

   !-----------------------------------------------!
   integer function ElocAndPsiTermsLin_nDataALL(this)
   !-----------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_nDataALL = ElocAndPsiTermsBase_nDataALL(this%EPTB)
   end function ElocAndPsiTermsLin_nDataALL

   !----------------------------------------------!
   integer function ElocAndPsiTermsLin_nParams(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_nParams = ElocAndPsiTermsBase_nParams(this%EPTB)
   end function ElocAndPsiTermsLin_nParams

   !----------------------------------------------!
   integer function ElocAndPsiTermsLin_nPCI(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_nPCI = this%EPTB%wfp%nCIParams
   end function ElocAndPsiTermsLin_nPCI

   !----------------------------------------------!
   integer function ElocAndPsiTermsLin_nPJ(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_nPJ = this%EPTB%wfp%nJParams
   end function ElocAndPsiTermsLin_nPJ
   !----------------------------------------------!
   integer function ElocAndPsiTermsLin_nPMO(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_nPMO = this%EPTB%wfp%nMOParams
   end function ElocAndPsiTermsLin_nPMO

   !----------------------------------------------------!
   function ElocAndPsiTermsLin_getWFP(this) result(wfp_p)
   !----------------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   type(WFParamDef), pointer :: wfp_p
   wfp_p => ElocAndPsiTermsBase_getWFP(this%EPTB)
   end function ElocAndPsiTermsLin_getWFP


   subroutine ElocAndPsiTermsLin_add(this,wfpDT)
   !-------------------------------------------!
      type(ElocAndPsiTermsLin), intent(inout) :: this
      type(WFParamDerivTerms), intent(in)     :: wfpDT
      integer np,j

      call ElocAndPsiTermsBase_add(this%EPTB,wfpDT)
      np = this%EPTB%wfp%nParams

      if (np /= size(wfpDT%fi)) call abortp("ElocAndPsiTermsLin_add: inconsistent sizes")

      this%fi(:) = this%fi(:) + wfpDT%fi(:)
      this%ELi(:) = this%ELi(:) + wfpDT%ELi(:)
      this%fiEL(:) = this%fiEL(:) + wfpDT%fi(:) * wfpDT%eloc

      do j=1,np
         this%fifj(:,j) = this%fifj(:,j) + wfpDT%fi(:) * wfpDT%fi(j)
         this%fifjEL(:,j) = this%fifjEL(:,j) + wfpDT%fi(:) * wfpDT%fi(j) * wfpDT%eloc
         this%fiELj(:,j) = this%fiELj(:,j) + wfpDT%fi(:) * wfpDT%ELi(j)
      enddo
   end subroutine ElocAndPsiTermsLin_add


   !-------------------------------------------!
   real(r8) function ElocAndPsiTermsLin_EMean(this)
   !-------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_EMean = ElocAndPsiTermsBase_Emean(this%EPTB)
   end function ElocAndPsiTermsLin_EMean

   !----------------------------------------------!
   real(r8) function ElocAndPsiTermsLin_EMeanALL(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_EMeanALL = ElocAndPsiTermsBase_EmeanALL(this%EPTB)
   end function ElocAndPsiTermsLin_EMeanALL

   !-----------------------------------------!
   real(r8) function ElocAndPsiTermsLin_var(this)
   !-----------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_var = ElocAndPsiTermsBase_var(this%EPTB)
   end function ElocAndPsiTermsLin_var

   !--------------------------------------------!
   real(r8) function ElocAndPsiTermsLin_varALL(this)
   !--------------------------------------------!
   type(ElocAndPsiTermsLin), intent(inout) :: this
   ElocAndPsiTermsLin_varALL = ElocAndPsiTermsBase_varALL(this%EPTB)
   end function ElocAndPsiTermsLin_varALL

   !-----------------------------------------------------------------------------------!
   subroutine ElocAndPsiTermsLin_resultALL(this,fi,ELi,fiEL,fifj,fifjEL,fiELj)
   !-----------------------------------------------------------------------------------!
   !
   ! note: copying of array unnecessary but required for reshape;
   ! better solution: how to transfer arrays (contiguous memory) with MPI?
   type(ElocAndPsiTermsLin), intent(inout) :: this
   real(r8), intent(out) :: fi(:),ELi(:),fiEL(:)  ! for notation see type def.
   real(r8), intent(out) :: fifj(:,:),fifjEL(:,:),fiELj(:,:)
   integer nData(1),nTotalData(1),np,np2,n,i,j
   real(r8) total(this%EPTB%wfp%nParams**2),vector(this%EPTB%wfp%nParams**2),vsum(this%EPTB%wfp%nParams)
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

   ! now the matrices with reshape to vector
   total = 0
   vector = reshape(this%fifj, (/ np2 /) )
   call myMPIReduceSumDouble(vector,total,np2)
   fifj = reshape(total, (/ np,np /) )
   fifj = fifj / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fifj(i,j) = fifj(j,i)
      enddo
   enddo

   vector = reshape(this%fifjEL, (/ np2 /) )
   call myMPIReduceSumDouble(vector,total,np2)
   fifjEL = reshape(total, (/ np,np /) )
   fifjEL = fifjEL / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fifjEL(i,j) = fifjEL(j,i)
      enddo
   enddo


   vector = reshape(this%fiELj, (/ np2 /) )
   call myMPIReduceSumDouble(vector,total,np2)
   fiELj = reshape(total, (/ np,np /) )
   fiELj = fiELj / n
   ! this is unsymmetric

   end subroutine ElocAndPsiTermsLin_resultALL

   !--------------------------------------------------------------------------------!
   subroutine ElocAndPsiTermsLin_result(this,fi,ELi,fiEL,fifj,fifjEL,fiELj)
   !--------------------------------------------------------------------------------!
   type(ElocAndPsiTermsLin), intent(in) :: this
   real(r8), intent(out) :: fi(:),ELi(:),fiEL(:)  ! for notation see type def.
   real(r8), intent(out) :: fifj(:,:),fifjEL(:,:),fiELj(:,:)
   integer n,i,j,np
   n = this%EPTB%nData
   np = size(fi)
   fi = this%fi / n
   ELi = this%ELi / n
   fiEL = this%fiEL / n
   fifj = this%fifj / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fifj(i,j) = fifj(j,i)
      enddo
   enddo
   fifjEL = this%fifjEL / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fifjEL(i,j) = fifjEL(j,i)
      enddo
   enddo
   fiELj = this%fiELj / n
   end subroutine ElocAndPsiTermsLin_result

end module elocAndPsiTermsLin_m
