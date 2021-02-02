! Copyright (C) 2015 Arne Luechow
! Copyright (C) 2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module rdataUpdate_m

   ! Rdata type storing R=(x1..xn,y1..yn,z1..zn) as well
   ! as data connected to R such as distances rai, rij
   ! phi(R), U(R), possibly the parameter derivatives fk, Uk
   ! this module is intended only for ECP calculation
   ! only function values, no electron derivatives
   ! one electron updates (but not moves)

   use kinds_m, only: r8
   use error_m
   use global_m
   use wfData_m, only: calcNucElecDists, calcElecElecDists, atoms, ncenter
   implicit none

   type, public :: RdataUpdate
      ! position data
      real(r8), pointer :: x(:) => null()      ! all electrons cartesian coordinates
      real(r8), pointer :: y(:) => null()      ! all electrons cartesian coordinates
      real(r8), pointer :: z(:) => null()      ! all electrons cartesian coordinates
      real(r8), pointer :: rai(:,:) => null()  ! elec-nuc distance matrix
      real(r8), pointer :: rij(:,:) => null()  ! elec-elec distance matrix
      real(r8)  :: xold=0, yold=0, zold=0        ! save old electron coordinates
      real(r8), pointer :: raiOld(:), rijOld(:)  ! save old electron distances
      real(r8), pointer :: raibarOld(:,:) => null()  ! save old \bar{r_{ai}}^v
      ! result data
      real(r8) :: phi=0.d0, U=0.d0             ! results: phi, U
      real(r8) :: phi0=0.d0, U0=0.d0           ! results: initial electron position: phi, U
      real(r8), pointer :: phik(:) => null()   ! parameter derivative \partial_k\Phi
      real(r8), pointer :: Uk(:) => null()     ! parameter derivative \partial_k U
      real(r8), pointer :: Uk0(:) => null()    ! initial parameter derivative \partial_k U
      real(r8), pointer :: mok(:) => null()    ! parameter derivative w.r.t. orbital rotation k
      ! intermediate data for jastrow updates
      real(r8), pointer :: Fij(:,:) => null()  ! jastrow two electron update terms (ee+een)
                                             ! storage: Fij(j,i) with i<j. Note order!
      real(r8), pointer :: Gi(:) => null()     ! jastrow one electron terms (en+aos)
      real(r8), pointer :: Fijold(:) => null() ! save original values
      real(r8), pointer :: Fijk(:,:,:) => null()  ! Fij for each parameter deriv k
      real(r8), pointer :: Fijkold(:,:) => null() ! save original values
      real(r8), pointer :: Gki(:,:) => null()     ! Gi for each parameter deriv k
      real(r8), pointer :: Gkiold(:) => null()    ! save originial values
      real(r8)  :: Giold=0, Xiold=0
      ! misc
      integer :: ieOld = 0               ! previous electron index for distances
      integer :: ieJasOld = 0            ! previous electron index for Jastrow arrays (Fi,Gi,Hi)
      integer :: enDim1 = 0              ! size of dim 1 in raibarOld
      integer :: enDim2 = 0              ! size of dim 2 in raibarOld
      integer :: npJ1 = 0, npJ2 = 0      ! number of one-electron and two-electron Jastrow parameter derivatives to calculate
      integer :: npCI = 0, npMO = 0      ! CI and MO parameter derivatives to calculate
      logical :: allValid = .false.
      logical :: phiValid = .false.
      logical :: jasValid = .false.
      logical :: ecpValid = .false.
      logical :: paramDerivs = .false.
      logical :: jasParamDerivs = .false.
      logical :: ciParamDerivs = .false.
      logical :: moParamDerivs = .false.
      logical :: allParamValid = .false.
      logical :: phiParamValid = .false.
      logical :: jasParamValid = .false.
      logical :: ecpParamValid = .false.
   contains
      procedure :: init => RdataUpdate_init
      procedure :: initWParamDerivs => RdataUpdate_initWParamDerivs
      procedure :: delete => RdataUpdate_delete
      procedure :: reset => RdataUpdate_reset
      procedure :: initENSize => RdataUpdate_initENSave
      procedure :: doParamDerivs => RdataUpdate_doParamDerivs
      procedure :: doJasParamDerivs => RdataUpdate_doJasParamDerivs
      procedure :: doCIParamDerivs => RdataUpdate_doCIParamDerivs
      procedure :: doMOParamDerivs => RdataUpdate_doMOParamDerivs
      procedure :: markJastrowValid => RdataUpdate_setValidJastrow
      procedure :: restoreElectron => RdataUpdate_restoreElectron
      procedure :: getOldElecIdx => RdataUpdate_getOldElecIdx
      procedure :: updateElectron => RdataUpdate_updateElectron
      procedure :: getUpdateResult => RdataUpdate_getUpdateResult
      procedure :: getNonLocalDerivs => RdataUpdate_getNonLocalDerivs

   end type RdataUpdate

contains

   subroutine RdataUpdate_init(this,x,y,z,rai,rij)
      class(RdataUpdate), intent(inout) :: this
      real(r8), intent(in)            :: x(:), y(:), z(:)
      real(r8), intent(in), optional  :: rai(:,:), rij(:,:)

      integer n,i,j
      call assert(.not.associated(this%x),"RdataUpdate_new: Rdata already allocated")
      call assert(size(x)==ne,"RdataUpdate_new: size(x) must much electron number")
      n = size(x)
      allocate(this%x(n),this%y(n),this%z(n))
      this%x = x; this%y = y; this%z = z

      allocate(this%rai(ncenter,n),this%rij(n,n))
      if (present(rai)) then
         this%rai = rai
      else
         call calcNucElecDists(this%x,this%y,this%z,this%rai)
      end if

      if (present(rij)) then
         this%rij = rij
      else
         call calcElecElecDists(this%x,this%y,this%z,this%rij)
      end if

      allocate(this%raiOld(ncenter),this%rijOld(n))
      this%raiOld = 0; this%rijOld = 0

      allocate(this%Fij(n,n), this%Fijold(n), this%Gi(n))
      this%Fij = 0; this%Fijold = 0; this%Gi = 0

      this%ieOld = 0; this%ieJasOld = 0
   end subroutine RdataUpdate_init

   subroutine RdataUpdate_initWParamDerivs(this,x,y,z,npJ1,npJ2,npCI,npMO,rai,rij)
      class(RdataUpdate), intent(inout) :: this
      real(r8), intent(in)            :: x(:), y(:), z(:)
      real(r8), intent(in), optional  :: rai(:,:), rij(:,:)
      integer, intent(in)           :: npJ1, npJ2, npCI, npMO

      integer n,i,j
      call assert(.not.associated(this%x),"RdataUpdate_new: Rdata already allocated")
      call assert(size(x)==ne,"RdataUpdate_new: size(x) must much electron number")
      n = size(x)
      allocate(this%x(n),this%y(n),this%z(n))
      this%x = x; this%y = y; this%z = z

      allocate(this%rai(ncenter,n),this%rij(n,n))
      if (present(rai)) then
         this%rai = rai
      else
         call calcNucElecDists(this%x,this%y,this%z,this%rai)
      end if

      if (present(rij)) then
         this%rij = rij
      else
         call calcElecElecDists(this%x,this%y,this%z,this%rij)
      end if

      allocate(this%raiOld(ncenter),this%rijOld(n))
      this%raiOld = 0; this%rijOld = 0

      allocate(this%Fij(n,n), this%Fijold(n), this%Gi(n))
      this%Fij = 0; this%Fijold = 0; this%Gi = 0

      this%paramDerivs = .true.
      this%npJ1 = npJ1
      this%npJ2 = npJ2
      this%jasParamDerivs = .true.
      if (npJ1+npJ2==0) this%jasParamDerivs = .false.
      this%npCI = npCI
      this%ciParamDerivs = .true.
      if (npCI==0) this%ciParamDerivs = .false.
      this%moParamDerivs = .true.
      this%npMO = npMO
      if (npMO==0) this%moParamDerivs = .false.

      allocate(this%Uk(npJ1+npJ2),this%Uk0(npJ1+npJ2),this%phik(npCI),this%mok(npMO))
      this%Uk = 0; this%Uk0 = 0; this%phik = 0; this%mok = 0

      if (this%jasParamDerivs) then
         allocate(this%Fijk(n,n,npJ2), this%Fijkold(n,npJ2), this%Gki(npJ1,n), this%Gkiold(npJ1))
      endif

      this%ieOld = 0; this%ieJasOld = 0
   end subroutine RdataUpdate_initWParamDerivs

   subroutine RdataUpdate_delete(this)
      class(RdataUpdate), intent(inout) :: this
      if (associated(this%x)) deallocate(this%x,this%y,this%z)
      if (associated(this%rai)) deallocate(this%rai,this%rij)
      if (associated(this%raiOld)) deallocate(this%raiOld, this%rijOld)
      if (associated(this%Fij)) deallocate(this%Fij,this%Fijold,this%Gi)
      if (associated(this%raibarOld)) deallocate(this%raibarOld)
      if (associated(this%Uk)) deallocate(this%Uk,this%Uk0,this%phik,this%mok)
      if (associated(this%Fijk)) deallocate(this%Fijk,this%Fijkold,this%Gki,this%Gkiold)
   end subroutine RdataUpdate_delete

   subroutine RdataUpdate_reset(this,x,y,z,rai,rij)
      class(RdataUpdate), intent(inout) :: this
      real(r8), intent(in)            :: x(:), y(:), z(:)
      real(r8), intent(in), optional  :: rai(:,:), rij(:,:)

      integer n,i,j
      call assert(associated(this%x),"RdataUpdate_new: Rdata not allocated")
      call assert(size(x)==ne,"RdataUpdate_new: size(x) must match electron number")
      n = size(x)
      this%x = x; this%y = y; this%z = z

      if (present(rai)) then
         this%rai = rai
      else
         call calcNucElecDists(this%x,this%y,this%z,this%rai)
      end if

      if (present(rij)) then
         this%rij = rij
      else
         call calcElecElecDists(this%x,this%y,this%z,this%rij)
      end if

      this%Fij = 0; this%Gi = 0

      if (associated(this%raibarOld)) this%raibarOld = 0

      this%ieOld = 0; this%ieJasOld = 0

      this%allValid = .false.
      this%phiValid = .false.
      this%jasValid = .false.
      this%ecpValid = .false.
      this%allParamValid = .false.
      this%phiParamValid = .false.
      this%jasParamValid = .false.
      this%ecpParamValid = .false.
   end subroutine RdataUpdate_reset

   subroutine RdataUpdate_initENSave(this,dim1,dim2)
      class(RdataUpdate), intent(inout) :: this
      integer, intent(in)               :: dim1
      integer, intent(in)               :: dim2

      if (associated(this%raibarOld)) then
         if (this%enDim1 /= dim1 .or. this%enDim2 /= dim2) then
            deallocate(this%raibarOld)
         endif
      endif
      if (.not.associated(this%raibarOld)) then
         allocate(this%raibarOld(dim1,dim2))
         this%enDim1 = dim1
         this%enDim2 = dim2
         this%raibarOld = 0
      endif
   end subroutine RdataUpdate_initENSave

   logical pure function RdataUpdate_doParamDerivs(this)
      class(RdataUpdate), intent(in) :: this
      RdataUpdate_doParamDerivs = this%paramDerivs
   end function RdataUpdate_doParamDerivs

   logical pure function RdataUpdate_doJasParamDerivs(this)
      class(RdataUpdate), intent(in) :: this
      RdataUpdate_doJasParamDerivs = this%jasParamDerivs
   end function RdataUpdate_doJasParamDerivs

   logical pure function RdataUpdate_doCIParamDerivs(this)
      class(RdataUpdate), intent(in) :: this
      RdataUpdate_doCIParamDerivs = this%ciParamDerivs
   end function RdataUpdate_doCIParamDerivs

   logical pure function RdataUpdate_doMOParamDerivs(this)
      class(RdataUpdate), intent(in) :: this
      RdataUpdate_doMOParamDerivs = this%moParamDerivs
   end function RdataUpdate_doMOParamDerivs

   integer pure function RdataUpdate_getOldElecIdx(this)
      class(RdataUpdate), intent(in) :: this
      RdataUpdate_getOldElecIdx = this%ieOld
   end function RdataUpdate_getOldElecIdx


   subroutine RdataUpdate_setValidJastrow(this)
      class(RdataUpdate), intent(inout) :: this
      this%jasValid = .true.
   end subroutine RdataUpdate_setValidJastrow


   subroutine RdataUpdate_updateElectron(this,ie,x,y,z)
      ! change coord and data (rai, rij)
      ! if new ie: restore previous values save current ie values
      ! this way the Rdu coord are always the original (new/reset call)
      ! except for the current ie
      class(RdataUpdate), intent(inout) :: this
      integer, intent(in)              :: ie
      real(r8), intent(in)               :: x,y,z
      integer n,a,j

      n = size(this%x)

      if (this%ieOld==0) then
         this%xold = this%x(ie); this%yold = this%y(ie); this%zold = this%z(ie)
         this%raiOld(1:ncenter) = this%rai(1:ncenter,ie)
         this%rijOld(1:ie-1) = this%rij(1:ie-1,ie)
         this%rijOld(ie+1:n) = this%rij(ie,ie+1:n)
         this%ieOld = ie
      else if (ie /= this%ieOld) then
         this%x(this%ieOld) =this%xold; this%y(this%ieOld) =this%yold; this%z(this%ieOld) =this%zold;
         this%xold = this%x(ie); this%yold = this%y(ie); this%zold = this%z(ie)
         this%rai(1:ncenter,this%ieOld) = this%raiOld(1:ncenter)
         this%rij(1:this%ieOld-1,this%ieOld) = this%rijOld(1:this%ieOld-1)
         this%rij(this%ieOld,this%ieOld+1:n) = this%rijOld(this%ieOld+1:n)
         this%raiOld(1:ncenter) = this%rai(1:ncenter,ie)
         this%rijOld(1:ie-1) = this%rij(1:ie-1,ie)
         this%rijOld(ie+1:n) = this%rij(ie,ie+1:n)
         this%ieOld = ie
      end if ! if (ie == ieOld) do no restoring!

      this%x(ie) = x; this%y(ie) = y; this%z(ie) = z

      do a=1,ncenter
         this%rai(a,ie) = sqrt( &
            (this%x(ie) - atoms(a)%cx)**2 + &
            (this%y(ie) - atoms(a)%cy)**2 + &
            (this%z(ie) - atoms(a)%cz)**2)
      enddo
      do j=1,ie-1
         this%rij(j,ie) = sqrt((this%x(ie)-this%x(j))**2 + &
            (this%y(ie)-this%y(j))**2 + (this%z(ie)-this%z(j))**2)
      enddo
      do j=ie+1,ne
         this%rij(ie,j) = sqrt((this%x(ie)-this%x(j))**2 + &
              (this%y(ie)-this%y(j))**2 + (this%z(ie)-this%z(j))**2)
      enddo
      this%allValid = .false.
      this%phiValid = .false.
      this%jasValid = .false.
      this%ecpValid = .false.
      this%allParamValid = .false.
      this%phiParamValid = .false.
      this%jasParamValid = .false.
      this%ecpParamValid = .false.
   end subroutine RdataUpdate_updateElectron


   subroutine RdataUpdate_restoreElectron(this)
      ! restore electron coordinates and distances to original
      class(RdataUpdate), intent(inout) :: this
      integer n

      n = size(this%x)

      if (this%ieOld > 0) then
         this%x(this%ieOld)=this%xold; this%y(this%ieOld)=this%yold; this%z(this%ieOld)=this%zold
         this%rai(1:ncenter,this%ieOld) = this%raiOld(1:ncenter)
         this%rij(1:this%ieOld-1,this%ieOld) = this%rijOld(1:this%ieOld-1)
         this%rij(this%ieOld,this%ieOld+1:n) = this%rijOld(this%ieOld+1:n)
      endif
      this%ieOld = 0
      this%allValid = .false.
      this%phiValid = .false.
      this%jasValid = .false.
      this%ecpValid = .false.
      this%allParamValid = .false.
      this%phiParamValid = .false.
      this%jasParamValid = .false.
      this%ecpParamValid = .false.
   end subroutine RdataUpdate_restoreElectron


   function RdataUpdate_getUpdateResult(this) result(res)
      class(RdataUpdate), intent(inout) :: this
      real(r8)                            :: res
      res = this%phi/this%phi0 * exp(this%U-this%U0)
   end function RdataUpdate_getUpdateResult


   subroutine RdataUpdate_getNonLocalDerivs(this,nld)
      class(RdataUpdate), intent(inout) :: this
      real(r8), intent(inout)             :: nld(:)
      integer                           :: idx,npJ,npCI,npMO,np
      npJ = this%npJ1 + this%npJ2
      npMO = this%npMO
      npCI = this%npCI
      idx = npJ+npMO
      np = npJ + npMO + npCI
      nld(1:npJ) = this%Uk * this%phi/this%phi0 * exp(this%U-this%U0)
      nld(npJ+1:idx) = this%mok/this%phi0 * exp(this%U-this%U0)
      nld(idx + 1:np) = this%phik/this%phi0 * exp(this%U-this%U0)

   end subroutine RdataUpdate_getNonLocalDerivs



end module rdataUpdate_m
