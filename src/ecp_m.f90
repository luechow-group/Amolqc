! Copyright (C) 1997 Arne Luechow
! Copyright (C) 2003/2004 Christian Diedrich
! Copyright (C) 2005 Annika Bande
! Copyright (C) 2013 Alexander Sturm
! Copyright (C) 2015 Arne Luechow
! Copyright (C) 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module ecp_m
   !

   ! TODO:
   ! check vnlk_create: jasparamdata_m.f90:30/24, wfparamdata_m.f90:101 optmode!!!
   ! why where ? calls necessary ? data in Rdu?


   ! effective core potentials: initialise, calculate localised V_ECP potential
   ! using Lebedev integration with the full trial wave function
   ! optionally parameter derivatives of V_ECP
   ! AL, 2015
   !
   !!use wfData_m
   use kinds_m, only: r8
   use global_m
   use jastrow_m, only: jasCalc, jasCalcInit, jasCalcUpdate, jasCalcInitWithUk, jasCalcUpdateWithUk
   use aoMo_m, only: calc_aomo
   use mos_m, only: mo1calc
   use multiDet_m, only: mdet2calc, resetToOld
   use random_m, only: myran
   use utils_m, only: double_findInSortedList
   use moParam_m, only:moparam_calcderivsOnlyMok,getMoUpdateMode
   use rdataUpdate_m
   use sphericalIntGrids_m

   implicit none
   private

   integer, parameter :: NOCUTOFF=1, NONLOCALCUTOFF=2, FULLCUTOFF=3
   integer, parameter, public :: TM_NONE=0, TM_SIMPLE=1, TM_SIZE_CONSISTENT=2, TM_SIZE_CONSISTENT1=3
   integer, parameter, public :: TM_W_BEFORE=0, TM_W_AFTER=1, TM_W_MEAN=2

   type, public :: PseudoAtomData
      integer :: a = 0           ! atom index in geometry
      integer :: nCoreElecs = 0  ! # of core elecs replaced by ecp
      integer :: ruleIdx = 0     ! start idx of spherical grid points in AllGridPoints
      integer :: intRule = 0     ! integration rule index (in NGridPointArray )
      integer :: nGridPoints = 4 ! grid points are: AllGridPoints(:,ruleIdx+1) ... AllGridPoints(:,ruleIdx+nGridPoints)
      integer :: lmax = 0        ! l_max
      real(r8)  :: cutoff = 0.d0   ! (distance) cutoff for ECP
      integer, allocatable    :: nterms(:)     ! number of (Gaussian) terms for each l
      integer, allocatable    :: nlk(:,:)      ! nlk (powers of r) for each l
      real(r8), allocatable     :: alk(:,:)      ! Gaussian coeffs A_lk, B_lk for each l
      real(r8), allocatable     :: blk(:,:)
   end type PseudoAtomData

   type, public :: TMoveDataType
      integer :: tmove = TM_NONE         ! T move type
      integer :: tmweight = TM_W_BEFORE  ! T move weighting option
      real(r8)  :: tau = 0.d0              ! time step in T move
      real(r8)  :: ecpLocalW               ! local potential for weighting step
      real(r8)  :: ecpNonlocalW            ! non local potential for weighting step
      logical :: isMoved = .false.       ! T move accepted?
   end type TMoveDataType

   type, public :: EcpType
      private
      ! type(MdetType), pointer     :: mdet_p
      ! type(JastrowType), pointer  :: jastrow_p
      type(PseudoAtomData), allocatable :: pAtom(:)        ! pseudo atom information
      integer                           :: cutoffMode = NOCUTOFF
      real(r8)                            :: cutoffThreshold = 0.d0
      logical                           :: randomRotation = .true.
      logical                           :: isInitialisedPA = .false.
      logical                           :: detOnlyLocalisation = .false.   ! localise ECP without Jastrow
      real(r8), allocatable               :: gridCoords(:,:,:)
   contains
      procedure :: setPseudoatoms => ecp_setPseudoatoms
      procedure :: getPseudoatoms => ecp_getPseudoatoms
      procedure :: getNPseudoatoms => ecp_getNPseudoatoms
      procedure :: isInitialised => ecp_isInitialised
      procedure :: nCoreElecs => ecp_nCoreElecs
      procedure :: setGrid => ecp_setGrid
      procedure :: getGrid => ecp_getGrid
      procedure :: setRandomRotation => ecp_setRandomRotation
      procedure :: isRandomRotation => ecp_isRandomRotation
      procedure :: setDetOnlyLocalisation => ecp_setDetOnlyLocalisation
      procedure :: isDetOnlyLocalisation => ecp_isDetOnlyLocalisation
      procedure :: setCutoffMode => ecp_setCutoffMode
      procedure :: getCutoffMode => ecp_getCutoffMode
      procedure :: setCutoffThreshold => ecp_setCutoffThreshold
      procedure :: getCutoffThreshold => ecp_getCutoffThreshold
      procedure :: setCutoff => ecp_setCutoff
      procedure :: getCutoff => ecp_getCutoff
      procedure :: calculate => ecp_calculate
      procedure :: printDebug => ecp_printDebug
   end type EcpType

   type :: uparr
         logical :: moved=.false.
         real(r8)  :: pos(3)=0
   end type uparr

contains


   subroutine ecp_setPseudoatoms(this,pseudoatoms)
      class(EcpType), intent(inout)    :: this
      type(PseudoAtomData), intent(in) :: pseudoatoms(:)           ! pseudo atom data structure

      if (allocated(this%pAtom)) deallocate(this%pAtom)
      !!!this%pAtom = pseudoatoms     ! default copy constructor with deep copy of components!
      allocate(this%pAtom(size(pseudoatoms)))
      this%pAtom = pseudoatoms

      this%isInitialisedPA = .true.
   end subroutine ecp_setPseudoatoms

   ! using automatic deallocation


   subroutine ecp_getPseudoatoms(this, pseudoatoms)
      class(EcpType), intent(in)          :: this
      type(PseudoAtomData), allocatable, intent(inout) :: pseudoatoms(:)           ! pseudo atom data structure

      if (allocated(pseudoatoms)) then
         if (size(pseudoatoms) /= size(this%pAtom)) deallocate(pseudoatoms)
      end if
      allocate(pseudoatoms(size(this%pAtom)))
      pseudoatoms  = this%pAtom
   end subroutine ecp_getPseudoatoms


   function ecp_getNPseudoatoms(this) result(res)
      class(EcpType), intent(in) :: this
      integer                    :: res
      res = 0
      if (allocated(this%pAtom)) res = size(this%pAtom)
   end function ecp_getNPseudoatoms


   function ecp_isInitialised(this) result(res)
      class(EcpType), intent(in) :: this
      logical                    :: res
      res = this%isInitialisedPA   ! .and. associated(jas) .and. associated(mdet)
   end function ecp_isInitialised


   pure function ecp_nCoreElecs(this,atom,pseudoatom) result(res)
      class(EcpType), intent(in)    :: this
      integer                       :: res          ! returns number of all removed core electrons,
                                                    ! or of one pseodo atom only
      integer, optional, intent(in) :: atom         ! number of atom in geometry
      integer, optional, intent(in) :: pseudoatom   ! number of pseudo atom
      integer p

      if (this%isInitialisedPA) then
         if (present(pseudoatom)) then
            res = this%pAtom(pseudoatom)%nCoreElecs
         else if (present(atom)) then
            res = 0
            do p = 1, size(this%pAtom)
               if (this%pAtom(p)%a == atom) then
                  res = this%pAtom(p)%nCoreElecs
                  exit
               end if
            end do
         else
            res = sum(this%pAtom(:)%nCoreElecs)
         end if
      else
         res = 0
      end if
   end function ecp_nCoreElecs


   subroutine ecp_setGrid(this,nGridPoints,atom,pseudoatom)
      class(EcpType), intent(inout) :: this
      integer, intent(in)           :: nGridPoints  ! number of grid points
      integer, optional, intent(in) :: atom         ! idx of atom in geometry
      integer, optional, intent(in) :: pseudoatom   ! idx of pseudo atom
      integer p, r, rule

      rule = 0
      do r=1,size(NGridPointsArray)
         if (nGridPoints == NGridPointsArray(r)) then
            rule = r
            exit
         endif
      enddo

      if (rule == 0) then
         write(iul,'(a)') ' illegal number of grid points given'
         write(iul,'(a)',advance='no') ' possible values: '
         write(iul,'(20i4)') NGridPointsArray
         call abortp("ecp_setGrid: illegal number of grid points")
      endif

      if (present(pseudoatom)) then
         this%pAtom(pseudoatom)%intRule = rule
         this%pAtom(pseudoatom)%nGridPoints = nGridPoints
      else if (present(atom)) then
         do p=1,size(this%pAtom)
            if (this%pAtom(p)%a == atom) then
               this%pAtom(p)%intRule = rule
               this%pAtom(p)%nGridPoints = nGridPoints
               exit
            endif
         enddo
      else
         this%pAtom(:)%intRule = rule
         this%pAtom(:)%nGridPoints = nGridPoints
      endif
   end subroutine ecp_setGrid

   function ecp_getGrid(this,atom,pseudoatom) result(nGridPoints)
      class(EcpType), intent(in)    :: this
      integer, optional, intent(in) :: atom         ! idx of atom in geometry
      integer, optional, intent(in) :: pseudoatom   ! idx of pseudo atom
      integer                       :: nGridPoints  ! number of grid points
      integer p

      if (present(pseudoatom)) then
         nGridPoints = this%pAtom(pseudoatom)%nGridPoints
      else if (present(atom)) then
         do p=1,size(this%pAtom)
            if (this%pAtom(p)%a == atom) then
               nGridPoints = this%pAtom(p)%nGridPoints
               exit
            endif
         enddo
      else
         call abortp("ecp_getGrid: either atom or pseudoatom argument required")
      endif
   end function ecp_getGrid

   function ecp_isRandomRotation(this) result(res)
      class(EcpType), intent(in) :: this
      logical                    :: res
      res = this%randomRotation
   end function ecp_isRandomRotation


   subroutine ecp_setRandomRotation(this,lval)
      class(EcpType), intent(inout)    :: this
      logical, optional, intent(in)    :: lval
      logical l
      l = .true.
      if (present(lval)) l = lval
      this%randomRotation = l
   end subroutine ecp_setRandomRotation


   subroutine ecp_getCutoffMode(this,cutoffMode)
      class(EcpType), intent(inout)    :: this
      character(len=*), intent(inout)  :: cutoffMode

      if (this%cutoffMode == NOCUTOFF) then
         cutoffMode = "no_cutoff"
      else if (this%cutoffMode == NONLOCALCUTOFF) then
         cutoffMode = "nonlocal_cutoff"
      else if (this%cutoffMode == FULLCUTOFF)  then
         cutoffMode = "full_cutoff"
      else
         call abortp("ecp_getCutoffMode: illegal cutoffMode")
      endif
   end subroutine ecp_getCutoffMode


   subroutine ecp_setCutoffMode(this,cutoffMode)
      class(EcpType), intent(inout)    :: this
      character(len=*), intent(in)     :: cutoffMode

      if (cutoffMode=="no_cutoff") then
         this%cutoffMode = NOCUTOFF
      else if (cutoffMode=="nonlocal_cutoff") then
         this%cutoffMode = NONLOCALCUTOFF
      else if (cutoffMode=="full_cutoff") then
         this%cutoffMode = FULLCUTOFF
      else
         call abortp("ecp_init: illegal cutoffMode argument")
      endif
   end subroutine ecp_setCutoffMode


   function ecp_getCutoff(this,atom,pseudoatom) result(cutoffValue)
      class(EcpType), intent(in)    :: this
      integer, optional, intent(in) :: atom         ! idx of atom in geometry
      integer, optional, intent(in) :: pseudoatom   ! idx of pseudo atom
      real(r8)                        :: cutoffValue
      integer p

      if (present(pseudoatom)) then
         cutoffValue = this%pAtom(pseudoatom)%cutoff
      else if (present(atom)) then
         do p=1,size(this%pAtom)
            if (this%pAtom(p)%a == atom) then
               cutoffValue = this%pAtom(p)%cutoff
               exit
            endif
         enddo
      else
         call abortp("ecp_getCutoffValue: either atom or pseudoatom argument required")
      endif
   end function ecp_getCutoff


   subroutine ecp_setCutoff(this, cutoffValue, atom, pseudoatom)
      class(EcpType), intent(inout)  :: this
      real(r8), intent(in)             :: cutoffValue
      integer, optional, intent(in)  :: atom         ! idx of atom in geometry
      integer, optional, intent(in)  :: pseudoatom   ! idx of pseudo atom
      integer p

      if (present(pseudoatom)) then
         this%pAtom(pseudoatom)%cutoff = cutoffValue
      else if (present(atom)) then
         do p=1,size(this%pAtom)
            if (this%pAtom(p)%a == atom) then
               this%pAtom(p)%cutoff = cutoffValue
               exit
            endif
         enddo
      else
         this%pAtom(:)%cutoff = cutoffValue
      endif

   end subroutine ecp_setCutoff


   function ecp_isDetOnlyLocalisation(this) result(res)
      class(EcpType), intent(in) :: this
      logical                    :: res
      res = this%detOnlyLocalisation
   end function ecp_isDetOnlyLocalisation


   subroutine ecp_setDetOnlyLocalisation(this,val)
      class(EcpType), intent(inout) :: this
      logical, optional, intent(in) :: val
      logical l
      l = .true.
      if (present(val)) l = val
      this%detOnlyLocalisation = l
   end subroutine ecp_setDetOnlyLocalisation


   subroutine ecp_setCutoffThreshold(this,cutoffThreshold)
   !------------------------------------------------------
      ! sets the cutoff threshold value (for V_l) and calculates the
      ! corresponding distance cutoffs for ECPs on the basis of
      ! One cutoff is used for all l channels simultaneously
      ! I would be possible to use individual cutoffs for each channel
      ! this version is based on the original version by Christian Diedrich
      class(EcpType), intent(inout)  :: this
      real(r8), intent(in)             :: cutoffThreshold
      integer :: ll,niter,p
      real(r8) :: tmp,curr_y,y1,y2,dydx
      real(r8), allocatable :: x(:)

      !inter: 1/2 of intervall for numerical derivation
      !trustr: trust radius
      !conv: convergence criterium
      real(r8) :: inter,conv,trustr
      parameter(inter=1d-5, conv=1d-10, trustr=0.2d0)
      integer, parameter :: maxiter=500

      this%cutoffThreshold = cutoffThreshold
      if (this%cutoffMode == NOCUTOFF) this%cutoffMode = NONLOCALCUTOFF

      do p = 1, size(this%pAtom)

         allocate(x(0:this%pAtom(p)%lmax))

         do ll = 0, this%pAtom(p)%lmax
            if (this%pAtom(p)%nterms(ll) == 1) then    ! some ECPs have a zero lmax term
               if (this%pAtom(p)%alk(ll,1) == 0.d0) then
                  x(ll) = 0.d0
                  exit
               endif
            endif

            x(ll) = 2.d0     !start guess taken from Fahy Phys.Rev.B 24,3503(1990)
            niter = 0

            call funk(ll,x(ll),curr_y)
            do
               if (niter > maxiter) then
                  write(iul,'(//a)') 'determination of distance cutoff for'
                  write(iul,*) 'non local ECP failed. niter > ',maxiter,ll
                  call abortp('(cutoff_ecp): ECP cutoff determination failed')
               endif
               niter = niter + 1
               call funk(ll,x(ll),curr_y)
               ! numerical derivative
               x(ll) = x(ll) + inter
               call funk(ll,x(ll),y2)
               x(ll) = x(ll) - 2*inter
               call funk(ll,x(ll),y1)
               x(ll) = x(ll) + inter
               dydx = (y2 - y1)/(2*inter)
               ! search at the descending tail of the function --> derivative should be negative
               if (dydx > 0.d0) then
                  x(ll) = x(ll)+trustr
               else
                  tmp = (cutoffThreshold-curr_y)/dydx + x(ll) ! linear estimate for new x-value
                  ! trust radius
                  if (abs(x(ll)-tmp) > trustr) then
                     ! modify x in the same direction like (cutoffthr-curr_y)/dydx
                     ! but only by the trust radius
                     x(ll) = sign(trustr,tmp-x(ll))+x(ll)
                  else
                     x(ll) = tmp
                  endif
               endif
               ! ensure that we are not in a local minimum
               ! which is very unlikely but not impossible
               if (abs(curr_y-cutoffThreshold) < conv) then
                  tmp = x(ll) + inter
                  call funk(ll,tmp,y1)
                  if (y1 < curr_y) then
                     exit
                  else
                     x(ll) = x(ll) + trustr
                  endif
               endif
            enddo
         enddo    ! ll loop

         if (this%cutoffMode == FULLCUTOFF) then
            this%pAtom(p)%cutoff = maxval(x(0:this%pAtom(p)%lmax))
         else if (this%cutoffMode == NONLOCALCUTOFF) then
            this%pAtom(p)%cutoff = maxval(x(0:this%pAtom(p)%lmax-1))
         else
            call abortp("ecp_setCutoffThreshold: illegal cutoffMode value")
         endif

         deallocate(x)

      enddo  ! p loop

   contains

      subroutine funk(ll,x,fx)     ! calculates |V_l(x)|
         integer :: ll, n
         real(r8) :: x, fx
         fx = 0
         do n = 1, this%pAtom(p)%nterms(ll)
            fx = fx + this%pAtom(p)%alk(ll, n) * exp(-this%pAtom(p)%blk(ll, n)*x*x)  &
                             * x**(this%pAtom(p)%nlk(ll, n) - 2)
         enddo
         fx = abs(fx)
      end subroutine funk

   end subroutine ecp_setCutoffThreshold


   function ecp_getCutoffThreshold(this) result(res)
      class(EcpType), intent(in) :: this
      real(r8) res
      res = this%cutoffThreshold
   end function ecp_getCutoffThreshold


   subroutine ecp_printDebug(this, iu)
      class(EcpType), intent(in) :: this
      integer, intent(in)        :: iu
      integer i
      write(iu,*) 'ecp_printDebug:'
      write(iu,*) 'alloc(pa):,', allocated(this%pAtom)
      if (allocated(this%pAtom)) then
         write(iu,*) 'size(pa):', size(this%pAtom), 'a,nCoreElecs:'
         do i = 1, size(this%pAtom)
            write(iu,*) this%pAtom(i)%a, this%pAtom(i)%nCoreElecs
         end do
      end if
      write(iu,*) 'cutoffmode:',this%cutoffMode
   end subroutine ecp_printDebug


!!!!!!!!!!!!!!!!!!!!!!!   ecp calculation   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   recursive subroutine ecp_calculate(this, Rdu, ecpLocal, ecpNonlocal, ecpNonlocalk, TMoveData)
   !--------------------------------------------------------------------------------------------
      ! calculates ECP for given electron coords Rdu.
      ! optionally calculates parameter derivative of nonlocal potential (VMC) when ecpNonlocalk given
      ! optionally calculates T moves (DMC) when tau is given
      class(EcpType), intent(inout)    :: this
      type(RdataUpdate), intent(inout) :: Rdu           ! data structure for electron update calculations
                                                        ! in: x,y,z coords and phi0, out: if T move: new coords
      real(r8), intent(inout)            :: ecpLocal      ! return local ECP energy V_local
      real(r8), intent(inout)            :: ecpNonlocal   ! and localized non-local contribution V_nonlocal (T move: at final coords)
      real(r8), optional, intent(inout)  :: ecpNonlocalk(:) ! parameter derivative of V_nonlocal
      type(TMoveDataType), optional, intent(inout) :: TMoveData  ! if given, do T move, out: ecpNonlocal for Weighting

      real(r8), allocatable              :: Vkk1(:,:)    ! V_k,k' discretised nonlocal potential (see Casula JCP 132, 154113)
      real(r8), allocatable  ::  ints(:), vl(:)
      real(r8)  rrai, tmp, ecpLocalI, ecpNonlocalI
      integer i, ll, lmax, pa, a, nParams, ngp, idx0, tngp

      type(uparr) :: tmscupdate(size(Rdu%x))

      if (present(ecpNonlocalk) .and. present(TMoveData)) then
         call abortp("ecp_calculate: T moves and parameter derivs not implemented")
      end if

      ecpLocal = 0.d0
      ecpNonlocal = 0.d0
      tngp = sum(this%pAtom(:)%nGridPoints)

      if (present(ecpNonlocalk)) ecpNonlocalk = 0.d0

      if (.not.this%detOnlyLocalisation) then
         if (Rdu%doJasParamDerivs()) then
            call jasCalcInitWithUk(Rdu)
         else
            call jasCalcInit(Rdu)
         endif
      end if

      if (present(TMoveData)) then
         TMoveData%isMoved = .false.
         TMoveData%ecpNonlocalW = 0.d0
         tngp = sum(this%pAtom(:)%nGridPoints)
         allocate(Vkk1(tngp,size(Rdu%x)))     ! (all grid points, # elecs)
         Vkk1 = 0.d0

         ! construct array for possible T move positions: all grid points
         if (.not.allocated(this%gridCoords)) then
            allocate(this%gridCoords(3,tngp,size(Rdu%x)))
         else
            if (.not.(size(this%gridCoords,2)==tngp .and. size(this%gridCoords,3)==size(Rdu%x))) &
               call abortp("ecp_calculate: inconsistent sizes in T move")
         end if
      end if

      ! loop over electrons
      do i = 1, size(Rdu%x)

         ! sum contributions for electron i
         ecpLocalI = 0.d0
         ecpNonlocalI = 0.d0
         idx0 = 0

         ! loop over pseudo atoms
         do pa = 1, size(this%pAtom)
            if (pa > 1) idx0 = idx0 + this%pAtom(pa-1)%nGridPoints
            lmax = this%pAtom(pa)%lmax
            a = this%pAtom(pa)%a
            rrai = Rdu%rai(a,i)

            !!!if (present(TMoveData)) print'(a,3i4,2g12.4)','DBG:ecpcalc:loop:',i,pa,a,rrai,this%pAtom(pa)%cutoff

            if (this%cutoffMode==FULLCUTOFF .and. rrai > this%pAtom(pa)%cutoff) cycle

            allocate(vl(0:lmax))
            ! calculate ecp coefficients v_l for pseudo atom pa
            call getVl(this, rrai, pa, lmax, vl)

            ecpLocalI = ecpLocalI + vl(lmax)      ! local part

            if (this%cutoffMode==NONLOCALCUTOFF .and. rrai > this%pAtom(pa)%cutoff) then
               deallocate(vl)
               cycle
            end if

            allocate(ints(0:lmax-1))

            if (present(ecpNonlocalk)) then ! parameter derivatives

               nParams = size(ecpNonlocalk)                     ! work around for gfortran 4.7
               block
                  real(r8) :: intsvNlk(0:lmax-1,nParams)

                  call ecp_paGridIntegration(this, Rdu, lmax, pa, i, ints, intsvNlk=intsvNlk)

                  do ll = 0, lmax - 1
                     tmp = (2*ll+1)/(4*pi) * vl(ll)
                     ecpNonlocalI = ecpNonlocalI +  tmp*ints(ll)
                     ecpNonlocalk(:) = ecpNonlocalk(:) + tmp*intsvNlk(ll,:)
                  end do
               end block

            else if (present(TMoveData)) then ! T moves

               ngp = this%pAtom(pa)%nGridPoints
               block
                  real(r8) :: Vkk1perElec(ngp,0:lmax-1)

                  Vkk1perElec = 0.d0
                  call ecp_paGridIntegration(this, Rdu, lmax, pa, i, ints, Vkk1perElec=Vkk1perElec, idxVkk=idx0)

                  !!!print'(a,g12.4,7i4)','DBG:ecpcalc:',rrai,a,i,pa,idx0,ngp,lmax,size(Vkk1perElec,1)
                  do ll = 0, lmax - 1
                     tmp = (2*ll+1)/(4*pi) * vl(ll)
                     ecpNonlocalI = ecpNonlocalI + tmp * ints(ll)
                     Vkk1(idx0+1:idx0+ngp,i) = Vkk1(idx0+1:idx0+ngp,i) + tmp * Vkk1perElec(:,ll)
                  enddo
               end block

            else ! neither T moves nor parameter derivatives

               call ecp_paGridIntegration(this, Rdu, lmax, pa, i, ints)

               do ll = 0, lmax - 1
                  ecpNonlocalI = ecpNonlocalI + (2*ll+1)/(4*pi) * vl(ll) * ints(ll)
               end do

            endif

            call Rdu%restoreElectron()
            call calc_aomo(i, Rdu%x, Rdu%y, Rdu%z, Rdu%rai)
            deallocate(ints,vl)

         end do ! PA

         if (present(TMoveData)) then
            if (TMoveData%tmove == TM_SIZE_CONSISTENT) call ecp_sizeConsistentTMove(this, Rdu, Vkk1, i, TMoveData)
            if (TMoveData%tmove == TM_SIZE_CONSISTENT1)  then
               call ecp_sizeConsistentTMove(this, Rdu, Vkk1, i, TMoveData,tmscupdate(i))
            endif
         end if

         ecpLocal = ecpLocal + ecpLocalI
         ecpNonlocal = ecpNonlocal + ecpNonlocalI

      end do ! electron

      if (present(TMoveData)) then
         if (TMoveData%tmove == TM_SIZE_CONSISTENT1) then
            do i = 1, size(Rdu%x)
              if (tmscupdate(i)%moved) then
                call Rdu%updateElectron(i, tmscupdate(i)%pos(1), tmscupdate(i)%pos(2), tmscupdate(i)%pos(3))
                call ecp_psit(this, Rdu, i, withParamDerivs=.false.)
                Rdu%phi0 = Rdu%phi; Rdu%U0 = Rdu%U
             endif
            enddo

         endif
         !!!print*,"DBG:ecpcalc:end:",ecpNonlocal,ecpLocal
         if (TMoveData%tmove == TM_SIMPLE) call ecp_TMoveSimple(this, Rdu, Vkk1, TMoveData)
      end if


   end subroutine ecp_calculate


   subroutine ecp_TMoveSimple(this, Rdu, Vkk1, TMoveData)
   !----------------------------------------------------------------------------------
      class(EcpType), intent(inout)    :: this
      type(RdataUpdate), intent(inout) :: Rdu           ! data structure for electron update calculations
      real(r8), intent(in)               :: Vkk1(:,:)     ! V_k,k' discretised nonlocal potential (see Casula JCP 132, 154113)
      type(TMoveDataType), intent(inout)  :: TMoveData
      real(r8) :: TNLpsum(size(Vkk1,1),size(Vkk1,2))        ! transition probabilities (partial sums)
      integer i, j, nElecs, nGridPoints, lb, ub, right, left
      real(r8) xi, vecElecs(-1:size(Vkk1,2)), vecGridPoints(0:size(Vkk1,1))
      real(r8) psum, noMoveProbability, norm

      nGridPoints = size(Vkk1,1)
      nElecs = size(Vkk1,2)

      !!!print*,"DBG:sum(Vkk1)=",sum(Vkk1), nGridPoints

      psum = 0.d0
      do i = 1, nElecs
         do j = 1, nGridPoints
            !!!print*,"DBG:tau*Vkk1:",i,j,TMoveData%tau*Vkk1(j,i)
            if (Vkk1(j,i) < 0) then
               psum = psum - TMoveData%tau * Vkk1(j,i)
            endif
            TNLpsum(j,i) = psum
         end do
      end do

      norm = 1.d0 + TNLpsum(nGridPoints,nElecs)
      TNLpsum = TNLpsum / norm
      noMoveProbability = 1.d0 / norm
      TNLpsum = TNLpsum + noMoveProbability

      !!!print*,'DBG:TNLpsum:per electron:'
      !!!do i = 1, nElecs
      !!!   print'(i5,100f8.4)', i, TNLpsum(:,i)
      !!!end do

      xi = myran()

      vecElecs(-1) = 0.d0
      vecElecs(0) = noMoveProbability
      vecElecs(1:nElecs) = TNLpsum(nGridPoints,1:nElecs)   ! summed probability for last grid point per electron

      !!!print'(a,100f8.4)','DBG:vecElecs:',vecElecs(-1:nElecs)

      lb = -1; ub = nElecs
      left = -1; right = nElecs
      ! find electron to move or no move
      call double_findInSortedList(lb, ub, vecElecs, xi, left, right)


      !!!print*,'DBG:findElec:', xi, left, right


      if (right > 0) then   ! T move

         ! find grid point (or 0 == nomove)
         i = right
         vecGridPoints(1:nGridPoints) = TNLpsum(1:nGridPoints,i)
         vecGridPoints(0) = 0.d0   ! actually: TNLpsum(nGridPoints,i-1) for i > 1

         !!!print'(a,100f8.4)','DBG:vecGrid:',vecGridPoints(0:nElecs)

         lb = 0; ub = nGridPoints
         left = 0; right = nGridPoints
         call double_findInSortedList(lb, ub, vecGridPoints, xi, left, right)
         j = right

         ! do T move
         call Rdu%updateElectron(i, this%gridCoords(1,j,i), this%gridCoords(2,j,i), this%gridCoords(3,j,i))

         !!!print'(a,2i3,3f15.5)', "DBG:Tmove:coords:",i,j,this%gridCoords(1,j,i), this%gridCoords(2,j,i), this%gridCoords(3,j,i)

         call ecp_psit(this, Rdu, i, withParamDerivs=.false.)
         Rdu%phi0 = Rdu%phi; Rdu%U0 = Rdu%U

         TMoveData%isMoved = .true.

      else

         TMoveData%isMoved = .false.

      end if

   end subroutine ecp_TMoveSimple


   subroutine ecp_sizeConsistentTMove(this, Rdu, Vkk1, i, TMoveData,tmscupdate)
      class(EcpType), intent(inout)    :: this
      type(RdataUpdate), intent(inout) :: Rdu           ! data structure for electron update calculations
      real(r8), intent(in)               :: Vkk1(:,:)     ! V_k,k' discretised nonlocal potential (see Casula JCP 132, 154113)
      integer, intent(in)              :: i             ! current electron
      type(TMoveDataType), intent(inout)  :: TMoveData
      type(uparr),intent(out),optional     :: tmscupdate
      integer lb, ub, left, right, j, nGridPoints
      real(r8) psum, xi
      real(r8) vecGridPoints(-1:size(Vkk1,1))

      if(present(tmscupdate)) then
         tmscupdate%pos=0
         tmscupdate%moved=.false.
      endif
      nGridPoints = size(Vkk1,1)

      vecGridPoints(-1) = 0.d0
      psum = 1.d0
      ! "diagonal element": no move for electron i
      vecGridPoints(0) = psum
      do j = 1, nGridPoints
         if (Vkk1(j,i) < 0) then
            psum = psum - TMoveData%tau * Vkk1(j,i)
         endif
         vecGridPoints(j) = psum
      end do

      vecGridPoints = vecGridPoints / vecGridPoints(nGridPoints)

      xi = myran()

      lb = -1; ub = nGridPoints
      left = -1; right = nGridPoints
      call double_findInSortedList(lb, ub, vecGridPoints, xi, left, right)

      if (right > 0) then
         ! T move
         j = right
         if (present(tmscupdate)) then
            tmscupdate%pos(1:3)=this%gridCoords(1:3,j,i)
            tmscupdate%moved=.true.
         else
         ! do T move
            call Rdu%updateElectron(i, this%gridCoords(1,j,i), this%gridCoords(2,j,i), this%gridCoords(3,j,i))
            call ecp_psit(this, Rdu, i, withParamDerivs=.false.)
            Rdu%phi0 = Rdu%phi; Rdu%U0 = Rdu%U
         end if

         TMoveData%isMoved = .true.
      end if
   end subroutine ecp_sizeConsistentTMove


   subroutine ecp_paGridIntegration(this, Rdu, lmax, pa, ie, ints, intsvNlk, Vkk1perElec, idxVkk)
      ! Y_l0 integrals with spherical grid around pseudo atom pa for electron i
      type(EcpType), intent(inout)     :: this
      type(RdataUpdate), intent(inout) :: Rdu             ! data structure for electron update calculations
      integer, intent(in)              :: lmax, pa, ie   ! current pseudo atom pa and current electron ie
      real(r8), intent(inout)            :: ints(0:)        ! integral values (for each l) by Lebedev quadrature
      real(r8), intent(inout), optional  :: intsvNlk(0:,:)  ! integral values for parameter derivatives
      real(r8), intent(inout), optional  :: Vkk1perElec(:,0:)  ! nonlocal contribution for electron ie at pa
      integer, intent(in), optional    :: idxVkk          ! start index in Vkk (first grid point)

      real(r8) :: phi, theta               ! axis orientation (polar coords)
      real(r8) :: rnew(3)                  ! new electron position
      real(r8) :: vGrid(3),vGridRot(3),ww  ! grid points, rotated grid points and weights
      real(r8) :: res
      real(r8) :: rot(3,3)                 ! random rotation matrix
      real(r8) :: costh, sinth, cosph, sinph
      real(r8) :: costhn                   ! cos(theta')
      real(r8) :: rai,raiVec(3)
      integer i,a,idx0,nParams

      a = this%pAtom(pa)%a

      if (this%randomRotation) then
         phi = myran()*2d0*pi
         theta = myran()*pi
      else                              ! fixed orientation of grid coord system
         phi = 0d0
         theta = 0d0
      endif

      ! random rotation matrix
      costh = cos(theta)
      sinth = sin(theta)
      cosph = cos(phi)
      sinph = sin(phi)
      rot(1,1) = cosph
      rot(1,2) = -sinph*costh
      rot(1,3) = sinph*sinth
      rot(2,1) = sinph
      rot(2,2) = cosph*costh
      rot(2,3) = -cosph*sinth
      rot(3,1) = 0d0
      rot(3,2) = sinth
      rot(3,3) = costh

      rai = Rdu%rai(a,ie)
      raiVec(1)  = Rdu%x(ie) - atoms(a)%cx
      raiVec(2)  = Rdu%y(ie) - atoms(a)%cy
      raiVec(3)  = Rdu%z(ie) - atoms(a)%cz

      ints(0:lmax-1) = 0.d0
      if (present(Vkk1perElec)) Vkk1perElec = 0.d0

      if (present(intsvNlk)) then
         nParams = size(intsvNlk,2)
         call internal_paGridIntegrationWithParamDerivs()
      else
         call internal_paGridIntegrationNoParamDerivs()
      endif

   contains

      subroutine internal_paGridIntegrationNoParamDerivs()
         idx0 = RuleIdx(this%pAtom(pa)%intRule)
         do i = 1, this%pAtom(pa)%nGridPoints
            vGrid = AllGridPoints(:,idx0+i)
            vGridRot = matmul(rot,vGrid)
            ww = AllWeights(idx0+i)

            ! Scale rule points by r_ai
            rnew(1) = atoms(a)%cx + rai*vGridRot(1)
            rnew(2) = atoms(a)%cy + rai*vGridRot(2)
            rnew(3) = atoms(a)%cz + rai*vGridRot(3)

            if (present(Vkk1perElec) .and. present(idxVkk)) then
               !!!print'(a,3i4,3f15.5)',"DBG:gi:",ie,idxVkk+i,i,rnew
               this%gridCoords(:,idxVkk + i,ie) = rnew
            end if

            if (lmax > 1) then ! calculate costhn with original positions
               costhn = dot_product(raiVec,vGridRot) / rai
            endif

            call Rdu%updateElectron(ie, rnew(1), rnew(2), rnew(3))
            call ecp_psit(this, Rdu, ie, withParamDerivs=.false.)

            res = Rdu%getUpdateResult()

            ww = ww * 4.d0*pi

            ! s part
            if(lmax > 0) then
               ints(0) = ints(0) + ww * res
               if (present(Vkk1perElec)) then
                  Vkk1perElec(i,0) = Vkk1perElec(i,0) + ww * res
               endif
            endif

            ! p part
            if(lmax > 1) then
               ints(1) = ints(1) + ww * costhn * res
               if (present(Vkk1perElec)) then
                  Vkk1perElec(i,1) = Vkk1perElec(i,1) + ww * costhn * res
               endif
            endif

            ! d part
            if(lmax > 2) then
               ints(2) = ints(2) + ww * (1.5d0 * costhn * costhn - 0.5d0) * res
               if (present(Vkk1perElec)) then
                  Vkk1perElec(i,2) = Vkk1perElec(i,2) + ww * (1.5d0 * costhn * costhn - 0.5d0) * res
               endif
            endif

         enddo

      end subroutine internal_paGridIntegrationNoParamDerivs

      subroutine internal_paGridIntegrationWithParamDerivs()
         real(r8) :: res,resNlk(nParams)

         intsvNlk(:,:) = 0
         idx0 = RuleIdx(this%pAtom(pa)%intRule)
         do i = 1, this%pAtom(pa)%nGridPoints
            vGrid = AllGridPoints(:,idx0+i)   ! i-th spherical grid point
            vGridRot = matmul(rot,vGrid)      ! rotation of grid point
            ww = AllWeights(idx0+i)           ! corresponding weight

            ! Scale rule points by r_ai
            rnew(1) = atoms(a)%cx + rai*vGridRot(1)
            rnew(2) = atoms(a)%cy + rai*vGridRot(2)
            rnew(3) = atoms(a)%cz + rai*vGridRot(3)

            if (lmax > 1) then ! calculate costhn with original positions
               costhn = dot_product(raiVec,vGridRot) / rai
            endif

            call Rdu%updateElectron(ie, rnew(1), rnew(2), rnew(3))
            call ecp_psit(this, Rdu, ie, withParamDerivs=.true.)

            res = Rdu%getUpdateResult()
            call Rdu%getNonlocalDerivs(resNlk)

            ww = ww * 4.d0*pi

            ! s part
            if(lmax > 0) then
               ints(0) = ints(0) + ww * res
               intsvNlk(0,:) = intsvNlk(0,:) + ww * resNlk(:)
            endif

            ! p part
            if(lmax > 1) then
               ints(1) = ints(1) + ww * costhn * res
               intsvNlk(1,:) = intsvNlk(1,:) + ww * costhn * resNlk(:)
            endif

            ! d part
            if(lmax > 2) then
               ints(2)       = ints(2) + ww * (1.5d0 * costhn * costhn - 0.5d0) * res
               intsvNlk(2,:) = intsvNlk(2,:) + ww * (1.5d0 * costhn * costhn - 0.5d0) * resNlk(:)
            endif
         enddo

      end subroutine internal_paGridIntegrationWithParamDerivs

   end subroutine ecp_paGridIntegration

   subroutine ecp_psit(this,Rdu,ie,withParamDerivs)
      type(EcpType), intent(in)        :: this
      type(RdataUpdate), intent(inout) :: Rdu             ! data structure for electron update calculations
      integer, intent(in)              :: ie
      logical, intent(in)              :: withParamDerivs
      real(r8) phi
      integer optmode
      optmode = 1

      if (withParamDerivs) then
         if (.not.this%detOnlyLocalisation) then
            if (Rdu%npJ1 + Rdu%npJ2 > 0) then
               call jasCalcUpdateWithUk(Rdu, ie)
            else
               call jasCalcUpdate(Rdu, ie)
            endif
         endif
         call calc_aomo(ie, Rdu%x, Rdu%y, Rdu%z, Rdu%rai)
         if (Rdu%npCI > 0) then
            if (Rdu%npMO > 0) then
               call mdet2calc(ie,phi,fk=Rdu%phik,ci_param_mode=optmode,resetDet=.false.,moUpdateMode=getMoUpdateMode())
               call moparam_calcderivsOnlyMok(Rdu%mok,ie)
               call resetToOld(ie)
            else
               call mdet2calc(ie,phi,fk=Rdu%phik,ci_param_mode=optmode,resetDet=.true.)

            endif
         else if (Rdu%npMO > 0) then ! MO without CI
            call mdet2calc(ie,phi,resetDet=.false.,moUpdateMode=getMoUpdateMode())
            call moparam_calcderivsOnlyMok(Rdu%mok,ie)
            call resetToOld(ie)
         else
            call mdet2calc(ie, phi)
         endif
      else
         if (.not.this%detOnlyLocalisation) then
            call jasCalcUpdate(Rdu,ie)
         endif
         call calc_aomo(ie, Rdu%x, Rdu%y, Rdu%z, Rdu%rai)
         call mdet2calc(ie, phi)
      endif
      Rdu%phi = phi
   end subroutine ecp_psit


   subroutine getVl(this, rai, pa, lmax, vl)
      type(EcpType), intent(in) :: this
      real(r8), intent(in)        :: rai         ! elec-nuc distance
      integer, intent(in)       :: pa, lmax   ! current pseudo atom with lmax
      real(r8), intent(inout)     :: vl(0:)      ! ecp coefficients V_l

      real(r8) :: r2
      integer :: ll, n

      vl(0:) = 0.d0
      r2 = rai * rai
      do ll = 0, lmax
         do n = 1, this%pAtom(pa)%nterms(ll)
           vl(ll) = vl(ll) + this%pAtom(pa)%alk(ll, n) * exp(-this%pAtom(pa)%blk(ll, n) * r2) &
                             * rai**(this%pAtom(pa)%nlk(ll, n) - 2)
         enddo
      enddo
   end subroutine getVl


end module ecp_m
