! Copyright (C) 2015 Arne Luechow
! Copyright (C) 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module elocAndPsiTermsENR_m
   ! implements ElocAndPsiTerms base class using delegation
   ! for Newton/Raphson  energy minimization
   ! notation follows Toulouse/Umrigar e.g. JCP 2007, and AL's notes on QMC algorithms
   use kinds_m, only: r8
   use global_m
   use wfParameters_m
   use elocAndPsiTermsBase_m
   use mpiInterface_m, only: myMPIReduceSumInteger, myMPIReduceSumDouble
   implicit none

   type ElocAndPsiTermsENR
      private
      type(ElocAndPsiTermsBase) :: EPTB
      real(r8), pointer :: fi(:) => null()       ! notation: Psi -> f, i.e. fi = Psi_i/Psi_0
      real(r8), pointer :: fiEL(:) => null()     ! fiEL = Psi_i/Psi_0*E_L
      real(r8), pointer :: ELi(:) => null()      ! ELi = E_L,i = \partial E_L / \partial p_i
      real(r8), pointer :: fij(:,:) => null()    ! fij = Psi_ij / Psi_0
      real(r8), pointer :: fifj(:,:) => null()   ! fifj = Psi_i/Psi_0 * Psi_j/Psi_0
      real(r8), pointer :: fijEL(:,:) => null()
      real(r8), pointer :: fifjEL(:,:) => null()
      real(r8), pointer :: fiELj(:,:) => null()
   end type ElocAndPsiTermsENR

private
public :: ElocAndPsiTermsENR, &
          ElocAndPsiTermsENR_create, ElocAndPsiTermsENR_destroy, ElocAndPsiTermsENR_reset, ElocAndPsiTermsENR_nData, &
          ElocAndPsiTermsENR_add, ElocAndPsiTermsENR_EMean, ElocAndPsiTermsENR_result, ElocAndPsiTermsENR_nParams, &
          ElocAndPsiTermsENR_EMeanALL, ElocAndPsiTermsENR_varALL, ElocAndPsiTermsENR_resultALL, ElocAndPsiTermsENR_getWFP, &
          ElocAndPsiTermsENR_nPCI

contains

   !-------------------------------------------------!
   subroutine ElocAndPsiTermsENR_create(this,eRef,wfp)
   !-------------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   real(r8), intent(in)                      :: eRef
   type(WFParamDef), pointer                  :: wfp

   call ElocAndPsiTermsBase_create(this%EPTB,eRef,wfp)
   allocate(this%fi(wfp%nParams),this%ELi(wfp%nParams),this%fiEL(wfp%nParams))
   allocate(this%fij(wfp%nParams,wfp%nParams),this%fifj(wfp%nParams,wfp%nParams))
   allocate(this%fijEL(wfp%nParams,wfp%nParams),this%fifjEL(wfp%nParams,wfp%nParams),this%fiELj(wfp%nParams,wfp%nParams))
   this%fi = 0; this%ELi = 0; this%fiEL = 0; this%fij = 0; this%fifj = 0
   this%fijEL = 0; this%fifjEL = 0; this%fiELj = 0
   end subroutine ElocAndPsiTermsENR_create

   !----------------------------------------!
   subroutine ElocAndPsiTermsENR_destroy(this)
   !----------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   deallocate(this%fi,this%ELi,this%fiEL,this%fij,this%fifj)
   deallocate(this%fijEL,this%fifjEL,this%fiELj)
   end subroutine ElocAndPsiTermsENR_destroy

   !--------------------------------------!
   subroutine ElocAndPsiTermsENR_reset(this)
   !--------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   call ElocAndPsiTermsBase_reset(this%EPTB)
   this%fi = 0; this%Eli = 0; this%fiEL = 0; this%fij = 0; this%fifj = 0
   this%fijEL = 0; this%fifjEL = 0; this%fiELj = 0
   end subroutine ElocAndPsiTermsENR_reset

   !---------------------------------------------!
   integer function ElocAndPsiTermsENR_nData(this)
   !---------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   ElocAndPsiTermsENR_nData = ElocAndPsiTermsBase_nData(this%EPTB)
   end function ElocAndPsiTermsENR_nData

   !-----------------------------------------------!
   integer function ElocAndPsiTermsENR_nDataALL(this)
   !-----------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   ElocAndPsiTermsENR_nDataALL = ElocAndPsiTermsBase_nDataALL(this%EPTB)
   end function ElocAndPsiTermsENR_nDataALL

   !----------------------------------------------!
   integer function ElocAndPsiTermsENR_nParams(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   ElocAndPsiTermsENR_nParams = ElocAndPsiTermsBase_nParams(this%EPTB)
   end function ElocAndPsiTermsENR_nParams

   !----------------------------------------------!
   integer function ElocAndPsiTermsENR_nPCI(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   ElocAndPsiTermsENR_nPCI = this%EPTB%wfp%nCIParams
   end function ElocAndPsiTermsENR_nPCI


   !----------------------------------------------------!
   function ElocAndPsiTermsENR_getWFP(this) result(wfp_p)
   !----------------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   type(WFParamDef), pointer :: wfp_p
   wfp_p => ElocAndPsiTermsBase_getWFP(this%EPTB)
   end function ElocAndPsiTermsENR_getWFP


   subroutine ElocAndPsiTermsENR_add(this,wfpDT)
   !-------------------------------------------!
      type(ElocAndPsiTermsENR), intent(inout) :: this
      type(WFParamDerivTerms), intent(in)     :: wfpDT
      integer np,j

      call ElocAndPsiTermsBase_add(this%EPTB,wfpDT)
      np = this%EPTB%wfp%nParams

      if (np /= size(wfpDT%fi)) call abortp("ElocAndPsiTermsENR_add: inconsistent sizes")

      this%fi(:) = this%fi(:) + wfpDT%fi(:)
      this%ELi(:) = this%ELi(:) + wfpDT%ELi(:)
      this%fiEL(:) = this%fiEL(:) + wfpDT%fi(:) * wfpDT%eloc

      this%fij(:,:) = this%fij(:,:) + wfpDT%fij(:,:)
      this%fijEL(:,:) = this%fijEL(:,:) + wfpDT%fij(:,:) * wfpDT%eloc

      do j =1, np
         this%fifj(:,j) = this%fifj(:,j) + wfpDT%fi(:) * wfpDT%fi(j)
         this%fifjEL(:,j) = this%fifjEL(:,j) + wfpDT%fi(:) * wfpDT%fi(j) * wfpDT%eloc
         this%fiELj(:,j) = this%fiELj(:,j) + wfpDT%fi(:) * wfpDT%ELi(j)
      enddo
   end subroutine ElocAndPsiTermsENR_add


   !-------------------------------------------!
   real(r8) function ElocAndPsiTermsENR_EMean(this)
   !-------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   ElocAndPsiTermsENR_EMean = ElocAndPsiTermsBase_Emean(this%EPTB)
   end function ElocAndPsiTermsENR_EMean

   !----------------------------------------------!
   real(r8) function ElocAndPsiTermsENR_EMeanALL(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   ElocAndPsiTermsENR_EMeanALL = ElocAndPsiTermsBase_EmeanALL(this%EPTB)
   end function ElocAndPsiTermsENR_EMeanALL

   !-----------------------------------------!
   real(r8) function ElocAndPsiTermsENR_var(this)
   !-----------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   ElocAndPsiTermsENR_var = ElocAndPsiTermsBase_var(this%EPTB)
   end function ElocAndPsiTermsENR_var

   !--------------------------------------------!
   real(r8) function ElocAndPsiTermsENR_varALL(this)
   !--------------------------------------------!
   type(ElocAndPsiTermsENR), intent(inout) :: this
   ElocAndPsiTermsENR_varALL = ElocAndPsiTermsBase_varALL(this%EPTB)
   end function ElocAndPsiTermsENR_varALL

   !-----------------------------------------------------------------------------------!
   subroutine ElocAndPsiTermsENR_resultALL(this,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)
   !-----------------------------------------------------------------------------------!
   !
   ! note: copying of array unnecessary but required for reshape;
   ! better solution: how to transfer arrays (contiguous memory) with MPI?
   type(ElocAndPsiTermsENR), intent(inout) :: this
   real(r8), intent(out) :: fi(:),ELi(:),fiEL(:)  ! for notation see type def.
   real(r8), intent(out) :: fij(:,:),fifj(:,:),fijEL(:,:),fifjEL(:,:),fiELj(:,:)
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
   vector = reshape(this%fij, (/ np2 /) )
   call myMPIReduceSumDouble(vector,total,np2)
   fij = reshape(total, (/ np,np /) )
   fij = fij / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fij(i,j) = fij(j,i)
      enddo
   enddo

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

   vector = reshape(this%fijEL, (/ np2 /) )
   call myMPIReduceSumDouble(vector,total,np2)
   fijEL = reshape(total, (/ np,np /) )
   fijEL = fijEL / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fijEL(i,j) = fijEL(j,i)
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
   ! fiELj is unsymmetric!!

   end subroutine ElocAndPsiTermsENR_resultALL

   !--------------------------------------------------------------------------------!
   subroutine ElocAndPsiTermsENR_result(this,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)
   !--------------------------------------------------------------------------------!
   type(ElocAndPsiTermsENR), intent(in) :: this
   real(r8), intent(out) :: fi(:),ELi(:),fiEL(:)  ! for notation see type def.
   real(r8), intent(out) :: fij(:,:),fifj(:,:),fijEL(:,:),fifjEL(:,:),fiELj(:,:)
   integer n,np,i,j
   n = this%EPTB%nData
   np = size(fi)

   fi = this%fi / n
   ELi = this%ELi / n
   fiEL = this%fiEL / n

   fij = this%fij / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fij(i,j) = fij(j,i)
      enddo
   enddo

   fifj = this%fifj / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fifj(i,j) = fifj(j,i)
      enddo
   enddo

   fijEL = this%fijEL / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fijEL(i,j) = fijEL(j,i)
      enddo
   enddo

   fifjEL = this%fifjEL / n
   ! complete matrix
   do i=1,np-1
      do j=i+1,np
         fifjEL(i,j) = fifjEL(j,i)
      enddo
   enddo

   fiELj = this%fiELj / n   ! unsymmetric!!

   end subroutine ElocAndPsiTermsENR_result

   !----------------------------------------------------!
   subroutine ElocAndPsiTermsENR_gradALL(this,e0,fi,fiEL)
   !----------------------------------------------------!
   !
   ! energy and terms for gradient only
   type(ElocAndPsiTermsENR), intent(inout) :: this
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

   end subroutine ElocAndPsiTermsENR_gradALL

   !-------------------------------------------------!
   subroutine ElocAndPsiTermsENR_grad(this,e0,fi,fiEL)
   !-------------------------------------------------!
   type(ElocAndPsiTermsENR), intent(in) :: this
   real(r8), intent(out) :: e0
   real(r8), intent(out) :: fi(:),fiEL(:)  ! for notation see type def.
   integer n,np
   n = this%EPTB%nData
   np = size(fi)

   e0 = ElocAndPsiTermsBase_EMean(this%EPTB)
   fi = this%fi / n
   fiEL = this%fiEL / n

   end subroutine ElocAndPsiTermsENR_grad

end module elocAndPsiTermsENR_m
