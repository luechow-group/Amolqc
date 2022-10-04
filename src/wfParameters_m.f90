! Copyright (C) 2011, 2015-2016 Arne Luechow
! Copyright (C) 2012 Alexander Sturm
! Copyright (C) 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module wfParameters_m

   ! allows abstract access to wave function parameters
   ! as parameter vector
   ! and the definition of the parameters
   ! must be a singleton, since only one wf 'type' can
   ! be initialized in jastrow,mdet,mos at one time
   use kinds_m, only: r8
   use global_m
   use jastrow_m, only: getJastrowParamCnt,getJastrowParamVector,&
                      putJastrowParamVector
   use jastrowParamData_m
   use multiDetParam_m, only: mdetparam_create,mdetparam_destroy,getCIParamsCnt,getCIParamsVector,putCIParamsVector
   use moParam_m, only: moparam_create,moparam_destroy,getMOParamsCnt,getMOParamsVector,putMOParamsVector, &
                      withSymmetriseList, getNParamSymEntries, getParamSymmetriseList
   use utils_m, only: tokenize!, numDerivative, num2Derivative
   implicit none
   private

   type, public :: WFParamDef
      !!private
      character(len=9) :: optType = 'jastrow'   ! jastrow|ci|jas+ci|mo
      integer          :: optMode = 1           ! (optional) 'mode' for opt types
      integer          :: nParams = 0           ! total # params
      integer          :: nJParams = 0          ! # Jastrow params
      integer          :: nCIParams = 0         ! # CI params
      integer          :: nMOParams = 0         ! # MO params
      integer          :: nJ1 = 0               ! # linear Jastrow one electron terms
      integer          :: nJ2 = 0               ! # linear Jastrow two electron terms
      integer          :: nJnl = 0              ! # nonlinear Jastrow terms
   end type WFParamDef


   type, public :: WFParamDerivTerms
      real(r8), allocatable :: fi(:)        ! notation: Psi -> f, i.e. fi = Psi_i/Psi_0
      real(r8), allocatable :: ELi(:)       ! ELi = E_L,i = \partial E_L / \partial p_i
      real(r8), allocatable :: fij(:,:)     ! fij = Psi_ij / Psi_0
      real(r8)              :: eloc         ! E_local, necessary for some derivatives
      real(r8)              :: phi          ! necessary for weighting
      real(r8)              :: U            ! necessary for weighting
      logical             :: fiCalc = .true.
      logical             :: ELiCalc = .true.
      logical             :: fijCalc = .true.
      logical             :: noCalc = .false.
   end type WFParamDerivTerms

   public :: wfparams_init, wfparams_destroy, wfparams_set, wfparams_get, &
             wfparams_change, wfparams_getNJastrowParams, wfparams_getNDetParams, &
             wfparams_getNParams, wfparams_symmetriseMOs, cicoeffs_get


contains


   subroutine wfparams_init(this,optType,optMode,lines,nl)
   !-----------------------------------------------------!
      type(WFParamDef), intent(inout) :: this
      character(len=9), intent(in)  :: optType
      integer, intent(in)           :: optMode
      integer, intent(in), optional :: nl
      character(len=120), intent(in), optional :: lines(:)
      ! "local" instance, do not allocate any globals
      ! this allows the existance more than one local instance at the same time
      ! local instances must not be destroyed
      integer npJ,npCI,npMO, npJ1, npJ2, npJnl

      call assert(optType=='jastrow' .or. optType=='ci' .or. optType=='jas+ci' .or. optType=='mo' &
         .or. optType=='mo+ci' .or. optType=='jas+mo' .or. optType=='jas+mo+ci', &
         "wfparams_init: unknown opttype")
      this%optType = optType
      this%optMode = optMode
      npJ = 0; npJ1 = 0; npJ2 = 0; npJnl = 0; npCI = 0; npMO = 0

      if (optType == 'jastrow' .or. optType=='jas+ci' .or. optType=='jas+mo' .or. optType=='jas+mo+ci') then
         call getJastrowParamCnt(this%optMode,npJ1,npJ2,npJnl)
         npJ = npJ1 + npJ2 + npJnl
         if (logmode>=2) write(iul,*) 'initializing jastrow parameters with np=',npJ,' ne= ',ne
         call jastrowparamdata_create(npJ,ne)
      endif
      if (optType == 'ci' .or. optType=='jas+ci' .or. optType=='jas+mo+ci' .or. optType=='mo+ci') then
         call mdetparam_create(this%optMode)
         npCI = getCIParamsCnt()
         if (logmode>=2) write(iul,*) 'initializing ci parameters with np=',npCI
      endif
      if (optType == 'mo' .or. optType=='mo+ci' .or. optType=='jas+mo' .or. optType=='jas+mo+ci') then
         if (present(lines) .and. present(nl)) then
            call moparam_create(this%optMode,lines,nl)
            npMO = getMOParamsCnt()
            if (logmode>=2) write(iul,*) 'initializing orbital rotation parameters with np=',npMO
         else
            call abortp("wfparams_init: orbital optimization with optType=mo requires lines,nl ")
         endif
      endif
      this%nJ1 = npJ1; this%nJ2 = npJ2; this%nJnl = npJnl
      this%nJParams = npJ; this%nCIParams = npCI; this%nMOParams = npMO
      this%nParams = npJ + npCI + npMO
   end subroutine wfparams_init


   subroutine wfparams_destroy(this)
   !--------------------------------
      type(WFParamDef), intent(inout) :: this
      if (this%optType == 'jastrow' .or. this%optType=='jas+ci' .or. this%optType=='jas+mo' &
         .or. this%optType == 'jas+mo+ci') then
         call jastrowparamdata_destroy()
      endif
      if (this%optType == 'ci' .or. this%optType =='jas+ci' .or. this%optType =='mo+ci' &
          .or. this%optType == 'jas+mo+ci' ) then
         call mdetparam_destroy()
      endif
      if (this%optType == 'mo' .or. this%optType == 'jas+mo' .or. this%optType=='mo+ci' &
          .or. this%optType =='jas+mo+ci' ) then
         call moparam_destroy()
      endif
      this%optType = 'jastrow'
      this%optMode = 1
      this%nJParams = 0; this%nCIParams = 0; this%nMOParams = 0
      this%nJ1 = 0; this%nJ2 = 0; this%nJnl = 0
      this%nParams = 0
   end subroutine wfparams_destroy


   subroutine wfparams_getNJastrowParams(this,npJ,npJ1,npJ2,npJnl)
      type(WFParamDef), intent(in)     :: this
      integer, intent(inout)           :: npJ
      integer, intent(inout), optional :: npJ1, npJ2, npJnl
      npJ = this%nJParams
      if (present(npJ1) .and. present(npJ2) .and. present(npJnl)) then
         npJ1 = this%nJ1
         npJ2 = this%nJ2
         npJnl = this%nJnl
      endif
   end subroutine wfparams_getNJastrowParams


   subroutine wfparams_getNDetParams(this,npCI,npMO)
      type(WFParamDef), intent(in) :: this
      integer, intent(inout)       :: npCI, npMO
      npCI = this%nCIParams
      npMO = this%nMOParams
   end subroutine wfparams_getNDetParams


   integer function wfparams_getNParams(this)
      type(WFParamDef), intent(in) :: this
      wfparams_getNParams = this%nParams
   end function wfparams_getNParams



   subroutine wfparams_set(this,p,normCI)
   !-----------------------------!
      type(WFParamDef), intent(in) :: this
      real(r8), intent(in)         :: p(:)
      logical,optional,intent(in) :: normCI   !for CI coefs
      logical                    :: normilize !for CI coefs
      normilize=.false.
      if(present(normCI)) normilize=.true.


      if (this%optType == 'jastrow') then
         call putJastrowParamVector(this%optMode,p)
      else if (this%optType == 'ci') then
         call putCIParamsVector(this%optMode,p,normilize)
      else if (this%optType == 'jas+ci') then
         call putJastrowParamVector(this%optMode,p(:this%nJParams))
         call putCIParamsVector(this%optMode,p(this%nJparams+1:),normilize)
      else if (this%optType == 'mo') then
         call putMOParamsVector(this%optMode,p)
      else if (this%optType == 'jas+mo') then
         call putJastrowParamVector(this%optMode,p(:this%nJParams))
         call putMOParamsVector(this%optMode,p(this%nJparams+1:))
      else if (this%optType == 'mo+ci') then
         call putMOParamsVector(this%optMode,p(:this%nMOParams))
         call putCIParamsVector(this%optMode,p(this%nMOparams+1:),normilize)
      else if (this%optType == 'jas+mo+ci') then
         call putJastrowParamVector(this%optMode,p(:this%nJParams))
         call putMOParamsVector(this%optMode,p(this%nJparams+1:(this%nJParams+this%nMOParams)))
         call putCIParamsVector(this%optMode,p((this%nJParams+this%nMOParams)+1:),normilize)
      else
         call abortp('wfparams_set: optType not implemented')
      endif
   end subroutine wfparams_set


   function wfparams_get(this)
   !-----------------------------!
      type(WFParamDef), intent(in) :: this
      real(r8)                     :: wfparams_get(this%nParams)
      real(r8)                     :: p(this%nParams)
      if (this%optType == 'jastrow') then
         call getJastrowParamVector(this%optMode,p)
      else if (this%optType == 'ci') then
         call getCIParamsVector(this%optMode,p)
      else if (this%optType == 'jas+ci') then
         call getJastrowParamVector(this%optMode,p(:this%nJParams))
         call getCIParamsVector(this%optMode,p(this%nJparams+1:))
      else if (this%optType == 'mo') then
         call getMOParamsVector(this%optMode,p)
      else if (this%optType == 'jas+mo') then
         call getJastrowParamVector(this%optMode,p(:this%nJParams))
         call getMOParamsVector(this%optMode,p(this%nJparams+1:))
      else if (this%optType == 'mo+ci') then
         call getMOParamsVector(this%optMode,p(:this%nMOParams))
         call getCIParamsVector(this%optMode,p(this%nMOparams+1:))
      else if (this%optType == 'jas+mo+ci') then
         call getJastrowParamVector(this%optMode,p(:this%nJParams))
         call getMOParamsVector(this%optMode,p(this%nJparams+1:(this%nJParams+this%nMOParams)))
         call getCIParamsVector(this%optMode,p((this%nJParams+this%nMOParams)+1:))
      else
         call abortp('wfparams_get: optType not implemented')
      endif
      wfparams_get = p
   end function wfparams_get

   !-----------------------------!
   function cicoeffs_get(this)
   !-----------------------------!
   type(WFParamDef), intent(in) :: this
   real(r8)                     :: cicoeffs_get(this%nParams+1)
   real(r8)                     :: p(this%nParams+1)
   if (this%optType /= 'ci') call abortp('cicoeffs_get: this function is only for pure ci optimization')
   if (this%optMode==3) call abortp('cicoeffs_get: this function can not be used for optmod=2 ')
   call getCIParamsVector(3,p)
   cicoeffs_get = p
   end function cicoeffs_get


   subroutine wfparams_change(lines,nl)
   !----------------------------------!
      integer, parameter :: MAXL=200
      integer, intent(in)            :: nl
      character(len=120), intent(in) :: lines(nl)
      character(len=40)              :: token(MAXL)
      real(r8),allocatable             :: p(:)
      integer optMode,iflag,np,io,i,ii,ntokens
      real(r8) pp
      character(len=9)      :: optType      ! 'jastrow'|'ci'
      type(WFParamDef)        :: WFP

      call assert(ne>0, ' wfparams_change: wave function must be initialized')
      call assert(nl==3, ' wfparams_change: three lines required, 2nd containing parameters')

      optType = 'jastrow'
      call getstra(lines,nl,'type=',optType,iflag)
      optMode = 1
      call getinta(lines,nl,'mode=',optMode,iflag)

      if (logmode >= 2) then
         write(iul,'(/3a,i3)') ' * * * Changing wf parameters of type=',optType,' and mode=',optMode
      endif

      call wfparams_init(WFP,optType,optMode)
      np = WFP%nParams
      call assert(np>0, ' wfParams_change: no parameters')

      allocate(p(np))
      p = wfparams_get(WFP)

      ! read 2nd line in the (free) format i1 p1 i2 p2 ... in pn
      call tokenize(lines(2),token,ntokens)

   !   print*,ntokens
   !   print*,lines(2)
   !   do i=1,ntokens,2
   !      print*,i,token(i),token(i+1)
   !   enddo

      do i=1,ntokens,2
         read(token(i),*) ii
         read(token(i+1),*) pp
         if (ii>=0.and.i<=np) then
            if (logmode >= 2) then
               write(iul,'(a,i5,2(a,g16.6))') ' changing ',ii,'-th parameter from ',p(ii),' to:',pp
            endif
            p(ii)=pp
         endif
      enddo

      call wfparams_set(WFP,p)

      call wfparams_destroy(WFP)
   end subroutine wfparams_change



   !---------------------------------------!
   subroutine wfparams_symmetriseMOs(self,p)
   !---------------------------------------!
      ! symmetrise MOs according to moparams:paramSymmetriseList
      type(WFParamDef), intent(in) :: self
      real(r8),intent(inout)       :: p(:)    ! (delta) parameter vector to symmetrise
      integer offset, idx, i, nsymp
      real(r8) psum, pmean
      integer, allocatable :: list(:)

      if (.not.withSymmetriseList()) return

      call assert(size(p)==self%nParams, "wfparams_symmetriseMOs: length mismatch of parameter vector")

      offset = 0

      if (self%optType == 'mo') then
         offset = 0
      else if (self%optType == 'jas+mo') then
         offset = self%nJParams
      else if (self%optType == 'mo+ci') then
         offset = 0
      else if (self%optType == 'jas+mo+ci') then
         offset = self%nJParams
      endif

      nsymp = getNParamSymEntries()
      do idx = 1, nsymp
         call getParamSymmetriseList(idx,list)
         ! replace all parameters whose index is in list with their mean
         psum = 0.d0
         do i = 1, size(list)
            psum = psum + p(offset+list(i))
         end do
         pmean = psum / size(list)
         do i = 1, size(list)
            p(offset+list(i)) = pmean
         end do
      end do

   end subroutine wfparams_symmetriseMOs


end module wfParameters_m
