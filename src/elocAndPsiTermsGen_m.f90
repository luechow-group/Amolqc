! Copyright (C) 2016 Arne Luechow
! Copyright (C) 2017 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module elocAndPsiTermsGen_m
   ! implements ElocAndPsiTerms base class using delegation
   ! for Newton/Raphson  energy minimization
   ! notation follows Toulouse/Umrigar e.g. JCP 2007, and AL's notes on QMC algorithms
   use kinds_m, only: r8
   use global_m
   use wfParameters_m
   use elocAndPsiTermsBase_m
   use mpiInterface_m, only: myMPIReduceSumInteger, myMPIReduceSumDouble
   implicit none

   type ElocAndPsiTermsGen
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
      logical         :: fiValid = .true.      ! all "adds" added fi terms
      logical         :: ELiValid = .true.     ! all "adds" added ELi terms
      logical         :: fijValid = .true.     ! all "adds" added fij terms
   end type ElocAndPsiTermsGen

private
public :: ElocAndPsiTermsGen, &
          ElocAndPsiTermsGen_create, ElocAndPsiTermsGen_destroy, ElocAndPsiTermsGen_reset, ElocAndPsiTermsGen_nData, &
          ElocAndPsiTermsGen_add, ElocAndPsiTermsGen_EMean, ElocAndPsiTermsGen_result, ElocAndPsiTermsGen_nParams, &
          ElocAndPsiTermsGen_EMeanALL, ElocAndPsiTermsGen_varALL, ElocAndPsiTermsGen_resultALL, ElocAndPsiTermsGen_getWFP, &
          ElocAndPsiTermsGen_nPCI,ElocAndPsiTermsGen_nPJ,ElocAndPsiTermsGen_nPMO

contains

   !-------------------------------------------------!
   subroutine ElocAndPsiTermsGen_create(this,eRef,wfp)
   !-------------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   real(r8), intent(in)                      :: eRef
   type(WFParamDef), pointer               :: wfp

   call ElocAndPsiTermsBase_create(this%EPTB,eRef,wfp)
   allocate(this%fi(wfp%nParams),this%ELi(wfp%nParams),this%fiEL(wfp%nParams))
   allocate(this%fij(wfp%nParams,wfp%nParams),this%fifj(wfp%nParams,wfp%nParams))
   allocate(this%fijEL(wfp%nParams,wfp%nParams),this%fifjEL(wfp%nParams,wfp%nParams),this%fiELj(wfp%nParams,wfp%nParams))
   this%fi = 0; this%ELi = 0; this%fiEL = 0; this%fij = 0; this%fifj = 0
   this%fijEL = 0; this%fifjEL = 0; this%fiELj = 0
   end subroutine ElocAndPsiTermsGen_create

   !----------------------------------------!
   subroutine ElocAndPsiTermsGen_destroy(this)
   !----------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   deallocate(this%fi,this%ELi,this%fiEL,this%fij,this%fifj)
   deallocate(this%fijEL,this%fifjEL,this%fiELj)
   end subroutine ElocAndPsiTermsGen_destroy

   !--------------------------------------!
   subroutine ElocAndPsiTermsGen_reset(this)
   !--------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   call ElocAndPsiTermsBase_reset(this%EPTB)
   this%fi = 0; this%Eli = 0; this%fiEL = 0; this%fij = 0; this%fifj = 0
   this%fijEL = 0; this%fifjEL = 0; this%fiELj = 0
   this%fiValid = .true.; this%ELiValid = .true.; this%fijValid = .true.
   end subroutine ElocAndPsiTermsGen_reset

   !---------------------------------------------!
   integer function ElocAndPsiTermsGen_nData(this)
   !---------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_nData = ElocAndPsiTermsBase_nData(this%EPTB)
   end function ElocAndPsiTermsGen_nData

   !-----------------------------------------------!
   integer function ElocAndPsiTermsGen_nDataALL(this)
   !-----------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_nDataALL = ElocAndPsiTermsBase_nDataALL(this%EPTB)
   end function ElocAndPsiTermsGen_nDataALL

   !----------------------------------------------!
   integer function ElocAndPsiTermsGen_nParams(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_nParams = ElocAndPsiTermsBase_nParams(this%EPTB)
   end function ElocAndPsiTermsGen_nParams
   !----------------------------------------------!
   integer function ElocAndPsiTermsGen_nPJ(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_nPJ = this%EPTB%wfp%nJParams
   end function ElocAndPsiTermsGen_nPJ
   !----------------------------------------------!
   integer function ElocAndPsiTermsGen_nPMO(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_nPMO = this%EPTB%wfp%nMOParams
   end function ElocAndPsiTermsGen_nPMO

   !----------------------------------------------!
   integer function ElocAndPsiTermsGen_nPCI(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_nPCI = this%EPTB%wfp%nCIParams
   end function ElocAndPsiTermsGen_nPCI


   !----------------------------------------------------!
   function ElocAndPsiTermsGen_getWFP(this) result(wfp_p)
   !----------------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   type(WFParamDef), pointer :: wfp_p
   wfp_p => ElocAndPsiTermsBase_getWFP(this%EPTB)
   end function ElocAndPsiTermsGen_getWFP


   subroutine ElocAndPsiTermsGen_add(this,wfpDT)
   !-------------------------------------------!
      type(ElocAndPsiTermsGen), intent(inout) :: this
      type(WFParamDerivTerms), intent(in)     :: wfpDT
      integer np,j

      call ElocAndPsiTermsBase_add(this%EPTB,wfpDT)
      np = this%EPTB%wfp%nParams

      if (this%fiValid) then
         if (wfpDT%fiCalc) then
            if (np /= size(wfpDT%fi)) call abortp("ElocAndPsiTermsGen_add: inconsistent sizes")
            this%fi(:) = this%fi(:) + wfpDT%fi(:)
            this%fiEL(:) = this%fiEL(:) + wfpDT%fi(:) * wfpDT%eloc
            do j =1, np
               this%fifj(:,j) = this%fifj(:,j) + wfpDT%fi(:) * wfpDT%fi(j)
               this%fifjEL(:,j) = this%fifjEL(:,j) + wfpDT%fi(:) * wfpDT%fi(j) * wfpDT%eloc
            enddo
         else
            this%fiValid = .false.
         end if
      end if

      if (this%fijValid) then
         if (wfpDT%fijCalc) then
            this%fij(:,:) = this%fij(:,:) + wfpDT%fij(:,:)
            this%fijEL(:,:) = this%fijEL(:,:) + wfpDT%fij(:,:) * wfpDT%eloc
         else
            this%fijValid = .false.
         end if
      end if

      if (this%ELiValid) then
         if (wfpDT%ELiCalc) then
            this%ELi(:) = this%ELi(:) + wfpDT%ELi(:)
            do j =1, np
               this%fiELj(:,j) = this%fiELj(:,j) + wfpDT%fi(:) * wfpDT%ELi(j)
            enddo
         end if
      end if

   end subroutine ElocAndPsiTermsGen_add


   !-------------------------------------------!
   real(r8) function ElocAndPsiTermsGen_EMean(this)
   !-------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_EMean = ElocAndPsiTermsBase_Emean(this%EPTB)
   end function ElocAndPsiTermsGen_EMean

   !----------------------------------------------!
   real(r8) function ElocAndPsiTermsGen_EMeanALL(this)
   !----------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_EMeanALL = ElocAndPsiTermsBase_EmeanALL(this%EPTB)
   end function ElocAndPsiTermsGen_EMeanALL

   !-----------------------------------------!
   real(r8) function ElocAndPsiTermsGen_var(this)
   !-----------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_var = ElocAndPsiTermsBase_var(this%EPTB)
   end function ElocAndPsiTermsGen_var

   !--------------------------------------------!
   real(r8) function ElocAndPsiTermsGen_varALL(this)
   !--------------------------------------------!
   type(ElocAndPsiTermsGen), intent(inout) :: this
   ElocAndPsiTermsGen_varALL = ElocAndPsiTermsBase_varALL(this%EPTB)
   end function ElocAndPsiTermsGen_varALL

   !-----------------------------------------------------------------------------------!
   subroutine ElocAndPsiTermsGen_resultALL(this,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)
   !-----------------------------------------------------------------------------------!
   !
   ! note: copying of array unnecessary but required for reshape;
   ! better solution: how to transfer arrays (contiguous memory) with MPI?
   type(ElocAndPsiTermsGen), intent(inout) :: this
   real(r8), optional, intent(out) :: fi(:)    ! for notation see type def.
   real(r8), optional, intent(out) :: ELi(:)
   real(r8), optional, intent(out) :: fiEL(:)
   real(r8), optional, intent(out) :: fij(:,:)
   real(r8), optional, intent(out) :: fifj(:,:)
   real(r8), optional, intent(out) :: fijEL(:,:)
   real(r8), optional, intent(out) :: fifjEL(:,:)
   real(r8), optional, intent(out) :: fiELj(:,:)
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
   if (present(fi)) then
      vsum = 0.d0
      call myMPIReduceSumDouble(this%fi,vsum,np)
      fi = vsum / n
   end if

   if (present(ELi)) then
      vsum = 0.d0
      call myMPIReduceSumDouble(this%ELi,vsum,np)
      ELi = vsum / n
   end if

   if (present(fiEL)) then
      vsum = 0.d0
      call myMPIReduceSumDouble(this%fiEL,vsum,np)
      fiEL = vsum / n
   end if

   ! now the matrices with reshape to vector
   if (present(fij)) then
      total = 0.d0
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
   end if

   if (present(fifj)) then
      total = 0.d0
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
   end if

   if (present(fijEL)) then
      total = 0.d0
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
   end if

   if (present(fifjEL)) then
      total = 0.d0
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
   end if

   if (present(fiELj)) then
      total = 0.d0
      vector = reshape(this%fiELj, (/ np2 /) )
      call myMPIReduceSumDouble(vector,total,np2)
      fiELj = reshape(total, (/ np,np /) )
      fiELj = fiELj / n
      ! fiELj is unsymmetric!!
   end if

   end subroutine ElocAndPsiTermsGen_resultALL

   !--------------------------------------------------------------------------------!
   subroutine ElocAndPsiTermsGen_result(this,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)
   !--------------------------------------------------------------------------------!
   type(ElocAndPsiTermsGen), intent(in) :: this
   real(r8), optional, intent(out) :: fi(:)    ! for notation see type def.
   real(r8), optional, intent(out) :: ELi(:)
   real(r8), optional, intent(out) :: fiEL(:)
   real(r8), optional, intent(out) :: fij(:,:)
   real(r8), optional, intent(out) :: fifj(:,:)
   real(r8), optional, intent(out) :: fijEL(:,:)
   real(r8), optional, intent(out) :: fifjEL(:,:)
   real(r8), optional, intent(out) :: fiELj(:,:)
   integer n,np,i,j
   n = this%EPTB%nData
   np = size(fi)

   if (present(fi))    fi = this%fi / n
   if (present(ELi))   ELi = this%ELi / n
   if (present(fiEL))  fiEL = this%fiEL / n

   if (present(fij)) then
      fij = this%fij / n
      ! complete matrix
      do i=1,np-1
         do j=i+1,np
            fij(i,j) = fij(j,i)
         enddo
      enddo
   end if

   if (present(fifj)) then
      fifj = this%fifj / n
      ! complete matrix
      do i=1,np-1
         do j=i+1,np
            fifj(i,j) = fifj(j,i)
         enddo
      enddo
   end if

   if (present(fijEL)) then
      fijEL = this%fijEL / n
      ! complete matrix
      do i=1,np-1
         do j=i+1,np
            fijEL(i,j) = fijEL(j,i)
         enddo
      enddo
   end if

   if (present(fifjEL)) then
      fifjEL = this%fifjEL / n
      ! complete matrix
      do i=1,np-1
         do j=i+1,np
            fifjEL(i,j) = fifjEL(j,i)
         enddo
      enddo
   end if

   if (present(fiELj)) then
      fiELj = this%fiELj / n   ! unsymmetric!!
   end if

   end subroutine ElocAndPsiTermsGen_result

   !----------------------------------------------------!
   subroutine ElocAndPsiTermsGen_gradALL(this,e0,fi,fiEL)
   !----------------------------------------------------!
   !
   ! energy and terms for gradient only
   type(ElocAndPsiTermsGen), intent(inout) :: this
   real(r8), intent(out) :: e0
   real(r8), intent(out) :: fi(:),fiEL(:)  ! for notation see type def.
   integer nData(1),nTotalData(1),np,np2,n,i,j
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

   end subroutine ElocAndPsiTermsGen_gradALL

   !-------------------------------------------------!
   subroutine ElocAndPsiTermsGen_grad(this,e0,fi,fiEL)
   !-------------------------------------------------!
   type(ElocAndPsiTermsGen), intent(in) :: this
   real(r8), intent(out) :: e0
   real(r8), intent(out) :: fi(:),fiEL(:)  ! for notation see type def.
   integer n,np,i,j
   n = this%EPTB%nData
   np = size(fi)

   e0 = ElocAndPsiTermsBase_EMean(this%EPTB)
   fi = this%fi / n
   fiEL = this%fiEL / n

   end subroutine ElocAndPsiTermsGen_grad

end module elocAndPsiTermsGen_m
