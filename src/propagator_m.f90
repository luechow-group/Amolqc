! Copyright (C) 1996, 2012, 2015-2016 Arne Luechow
! Copyright (C) 2013 Alexander Sturm
! Copyright (C) 2013, 2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module propagator_m

  use kinds_m, only: r8
  use error_m
  use random_m, only: myran, mygran
  use statistics
  use randomWalker_m
  use eloc_m, only: eloc
  use wfData_m, only: atoms
  implicit none

  private
  public :: propagator_init, propagator_reset, &
            propagator_writeParams, propagator_setERef0, propagateAndWeight, &
            propagator_endOfBlock, propagator_timeStep, propagator_setTimeStep, &
            propagator_timeStepAcc, propagator_getInitialTau,setNoExp,getAccRatio, &
            propagator_getTMove, propagator_getTMoveRatio,TMoveControl

  integer, parameter :: mMaxBlockSize = 8
  integer     :: mWeight=0            ! weighting strategy (none=0,Rey=1,Umr,Accept)
  integer     :: mMove=1              ! propagator (Rey=1,Umr)
  integer     :: mBlockDiscard=5           ! blocks to discard initially
  integer     :: mBlockGetTauEff=2         ! blocks (before discard) used in DMC to calc tau_a
  real(r8)      :: mTau=0.d0            ! time step
  real(r8)      :: mDriftScal=1.d0      ! option to scal drift term in VMC
  logical     :: mARStep=.true.       ! perform acceptance rejection step
  logical     :: mRejCross=.false.    ! reject nodal crossing?
  logical     :: mNoExp = .false.
  logical     :: mTauIsSet = .false.
  character*3 :: mMoveType='all'      ! 'all' or 'one' electron moves

  real(r8)           :: mSqrtTau=0

  type(RandomWalker), save :: mRWSave(mMaxBlockSize)      ! additional walker to save old positions
  type(RandomWalker), save :: mRWTMoves(mMaxBlockSize)    ! additional walker for T move positions:
                                                          ! -> allocatable !!
  logical, save            :: mIsMoved(mMaxBlockSize)

  logical          :: sAccumulateTauAcc=.false.
  real(r8)           :: sTauAcc=0
  real(r8)           :: sERef0=0              ! set from qmc
  real(r8)           :: sR2acc=0,sR2all=0     ! r**2 statistics for \tau_a (Reynolds alg)
  real(r8)           :: mAccRatio = 0         ! for VMC-Sample Prerun (S.Manten)
  type(simpleStat),save :: sTauaStat
  !logical     :: mNo_neg_nonlocal = .false.    ! excluds negative non local part in weighting step
  type(simpleStat) :: mTMoveCount
  !!type(simpleStat) :: mDBGELoc
  !!type(simpleStat) :: mDBGERef
  !!type(Doublehistogram) :: mDBGWeights

  integer :: mCounter

  type :: TMoveControl
      integer :: mode = 0     ! T move (none=0, simple=1, all_elecs=2)
      integer :: tmoverej = 0 !0: T move for rejected steps, 1: no T move for rejected steps, 2: T move at proposed position
      integer :: tmovew = 0   ! 0: E_L for weight before T move, 1: after (all) T move(s)
      integer :: tmovecross =0 ! 1: reject all crossings, 1: allow double crossing drift-diff then TMove crossing
      logical :: tmoveddeloc = .false.
   end type TMoveControl

  type(TMoveControl) :: mTMove ! -> mTMove in propagator

contains


  subroutine propagator_reset()
  !---------------------------!
    mWeight=0                 ! weighting strategy (none=0,Rey=1,Umr,Accept)
    mMove=0                   ! propagator (Rey=1,Umr)
    mTmove%mode=0             ! T move
    !mNo_neg_nonlocal = .false.! T move
    mBlockDiscard=5           ! blocks to discard initially
    mBlockGetTauEff=2         ! blocks (before discard) used in DMC to calc tau_a
    mTau=0.d0                 ! time step
    mTauIsSet = .false.
    mDriftScal=1.d0           ! option to scal drift term in VMC
    mARStep=.true.            ! perform acceptance rejection step
    mRejCross=.false.         ! reject nodal crossing?
    mMoveType='all'           ! 'all' or 'one' electron moves
    mNoExp = .false.
    mSqrtTau=0

    sAccumulateTauAcc=.false.
    sTauAcc=0
    sERef0=0                  ! set from qmc
    sR2acc=0;sR2all=0         ! r**2 statistics for \tau_a (Reynolds alg)
    mAccRatio=0
    call reset(sTauaStat)
    call reset(mTMoveCount)
    !!if (logmode >= 3) then
    !!  call createHistogram(mDBGWeights,0.5d0,1.5d0,20)
    !!  call reset(mDBGERef)
    !!  call reset(mDBGELoc)
    !!end if
  end subroutine propagator_reset


  subroutine propagator_init(wgt,move,Tmove,dc,bgte,tau,ds,ar,rc,mt,wb)
  !-------------------------------------------------------------------!
    ! initialize propagator. Only most important parameters here
    ! For the rest create individual set functions if required
    integer, intent(in)            :: wgt,move,dc,bgte,wb
    type(TMoveControl),intent(in)  :: Tmove
    real(r8), intent(in)             :: tau,ds
    logical, intent(in)            :: ar,rc
    character(len=3), intent(in)   :: mt
    integer :: alstat,i

    mWeight         = wgt
    mMove           = move

    if (Tmove%mode > 0 .and. move == 3) call abortp("propagator_init: T moves not allowed with two-level propagator")

      mTmove%mode       = Tmove%mode
      mTMove%tmoverej   = TMove%tmoverej
      mTMove%tmovew     = TMove%tmovew
      mTMove%tmovecross = TMove%tmovecross

    call reset(mTMoveCount)

    !!if (logmode >= 3) then
    !!  call reset(mDBGWeights)
    !!  call reset(mDBGERef)
    !!  call reset(mDBGELoc)
    !!end if

    mBlockDiscard   = dc
    mBlockGetTauEff = bgte
    mTau      = tau
    mTauIsSet = .true.
    sTauAcc   = mTau
    mSqrtTau  = sqrt(mTau)
    mDriftScal= ds
    mARStep   = ar
    mRejCross = rc
    mMoveType = 'all'
    if (mt == 'one') mMoveType = mt
    if (bgte < 0)then
       mBlockGetTauEff = 2
    else
       mBlockGetTauEff = bgte
    end if

    call assert(size(mRWSave) >= wb,'(propagator_init): mRWSave too small')
    if (.not. rw_isAllocated(mRWSave(1))) then
       do i=1,mMaxBlockSize
          call rw_new(mRWSave(i))
       enddo
    endif
    call assert(size(mRWTMoves) >= wb,'(propagator_init): mRWTmoves too small')
    if (.not. rw_isAllocated(mRWTMoves(1))) then
       do i=1,mMaxBlockSize
          call rw_new(mRWTMoves(i))
          mIsMoved(i) = .false.
       enddo
    endif

    mCounter = 0

  end subroutine propagator_init


  subroutine propagator_writeParams(iu)
  !-----------------------------------!
    integer, intent(in) :: iu
    character(len=12)   :: s1,s2,s3
    logical optout

    if (mWeight==0) s1='none'
    if (mWeight==1) s1='Reynolds'
    if (mWeight==2) s1='Umrigar'
    if (mWeight==3) s1='Accept'
    if (mWeight==4) s1='tmsm' !non symmetrical with tau
    if (mWeight==5) s1='tmsc' ! symmetrical with tau
    if (mMove==1)   s2='Reynolds'
    if (mMove==2)   s2='Umrigar'
    if (mMove==3)   s2='TwoLevel'
    if (mMove==4)   s2='Gaussian'
    if (mTmove%mode==0)  s3='none'
    if (mTmove%mode==1)  s3='simple'
    if (mTmove%mode==2)  s3='sc'
    if (mTmove%mode==3)  s3='sc1'

    write(iu,'(A/)') '    propagator parameters:'
    write(iu,'(1X,2(A21,A12,3X))') ' weight =',trim(s1), ' move =',trim(s2)
    if (mWeight > 0) then
      write(iu,'(1X,2(A21,A12,3X))') ' moveType =',mMoveType,' T_moves =',trim(s3)
    else
      write(iu,'(1X,2(A21,A12,3X))') ' moveType =',mMoveType
    end if
    write(iu,'(1X,2(A21,F12.5,3X))') ' tau =',mTau,' drift_scal =',mDriftScal
    write(iu,'(1X,2(A21,L12,3X))') ' AR step =',mARStep,' rej_cross =',mRejCross
    optout = .false.
    if (mMove==2) then
       write(iu,'(1X,A21,L12)',advance='no') ' expon. sampling =', .not. mNoExp
       optout = .true.
    end if
    if (mWeight > 0) then
       write(iu,'(1X,A21,I12,2X)',advance='no') ' T_move_reject =', mTMove%tmoverej
       write(iu,'(1X,A21,I12,2X)') ' T_move_wgt =', mTMove%tmovew
       write(iu,'(1X,A21,I12,2X)',advance='no') ' T_move_cross =', mTMove%tmovecross
       optout = .true.
    end if
    if (optout) write(iu,*)

    write(iu,*)
  end subroutine propagator_writeParams


  subroutine setNoExp()
  !--------------------!
   mNoExp = .true.
  end subroutine setNoExp


  subroutine getAccRatio(ratio)
  !---------------------------------!
   real(r8) :: ratio
   ratio = mAccRatio
  end subroutine getAccRatio


  subroutine propagator_setERef0(E)
  !-------------------------------!
    real(r8), intent(in) :: E
    sERef0 = E
  end subroutine propagator_setERef0


  subroutine propagator_setTimeStep(tau)
  !------------------------------------!
    real(r8), intent(in) :: tau
    mTau = tau
    mSqrtTau = sqrt(tau)
    sTauAcc = tau
    mTauIsSet = .true.
  end subroutine propagator_setTimeStep


  real(r8) function propagator_timeStep()
  !-----------------------------------!

    if (.not.mTauIsSet) then
       mTau = propagator_getInitialTau(0.5d0)
    end if
    propagator_timeStep = mTau
  end function propagator_timeStep


  real(r8) function propagator_timeStepAcc()
  !--------------------------------------!
    propagator_timeStepAcc = sTauAcc
  end function propagator_timeStepAcc


  real(r8) function propagator_getInitialTau(ar)
  !------------------------------------------!
     real(r8), intent(in) :: ar ! desired acceptance ratio
     integer i
     real(r8) Zmax,a

     call assert(ar > 0.d0,"propagator_getInitialTau: positive argument required")
     Zmax = 0
     do i=1,ncenter
        Zmax = max(atoms(i)%za,Zmax)
     enddo
     a = 5*Zmax
     propagator_getInitialTau = - log(ar) / a
  end function propagator_getInitialTau


  subroutine propagator_endOfBlock(b)
    ! tell propagator that is's end of block b
    integer, intent(in) :: b
    if (b==mBlockDiscard-mBlockGetTauEff-1) then
       sAccumulateTauAcc = .true.
       call reset(sTauaStat)
    endif
    if (b==mBlockDiscard) sAccumulateTauAcc = .false.
    if (b>=mBlockDiscard-mBlockGetTauEff .and. b<=mBlockDiscard) then
       call addData(sTauaStat,sR2acc/sR2all*mTau)
       sTauAcc = mean(sTauaStat)
    endif

    !!if (logmode >= 3) then
    !!   write(iull,'(a,3f15.5)') 'ELOC:', mean(mDBGELoc), variance(mDBGELoc), stdDevMean(mDBGELoc)
    !!   write(iull,'(a,3f15.5)') 'EREF:', mean(mDBGERef), variance(mDBGERef), stdDevMean(mDBGERef)
    !!   write(iull,'(a,3f15.5)') 'HIST:',mean(mDBGWeights), variance(mDBGWeights), stdDevMean(mDBGWeights)
    !!   call writeHistogram(mDBGWeights,iull)
    !!end if
  end subroutine propagator_endOfBlock


  subroutine propagateAndWeight(rwb,ERef,mode,accepted)
    type(RandomWalker), intent(inout) :: rwb(:)      ! walker block to be propagated
    real(r8), intent(in)                :: ERef        ! only for weighting required
    integer, intent(in)               :: mode        ! 0: default (i.e. use mARStep)
                                                     ! 1,2: with,without accept step
    real(r8), intent(out)               :: accepted    ! ratio of accepted steps
                                                     ! (0.0,1.0 for one all-electron moves
                                                     ! ratio for one electron moves or block moves)

    if (mMoveType=='all') then          ! all electron moves
       call propagateAndWeightAll(rwb,ERef,mode,accepted)
    else if (mMoveType=='one') then     ! one electron moves
       call assert(size(rwb)==1,'propagateAndWeight: one electron move not yet implemented for walker blocks')
       call propagateAndWeightOne(rwb(1),ERef,mode,accepted)
    else
       call abortp("propagateAndWeight: wrong moveType")
    endif
  end subroutine propagateAndWeight


  subroutine propagator_writeStatus(iu)
    integer, intent(in) :: iu
    continue
  end subroutine propagator_writeStatus


  subroutine propagator_readStatus(iu)
    integer, intent(in) :: iu
    continue
  end subroutine propagator_readStatus


  function propagator_getTMove() result(res)
    integer :: res
    res = mTMove%mode
  end function propagator_getTMove

  function propagator_getTMoveRatio() result(res)
    real(r8) :: res
    res = mean(mTMoveCount)
  end function propagator_getTMoveRatio



!----- private routines --------------------------------------------------



  subroutine propagateAndWeightAll(rwb,ERef,mode,accepted)
  !------------------------------------------------------!

    type(RandomWalker), intent(inout) :: rwb(:)      ! walker block to be propagated
    real(r8), intent(in)                :: ERef        ! only for weighting required
    integer, intent(in)               :: mode        ! 0: default (i.e. use mARStep)
                                                     ! 1,2: with,without accept step
    real(r8), intent(out)               :: accepted    ! 1.0 for accepted, 0.0 for rejected
                                                     ! ratio for walker blocks

    integer i, walkerBlockSize, tmove
    real(r8) xi,mp,s1,sn1, tau
    real(r8) r2(size(rwb))
    real(r8) vvr(size(rwb)),vvnr(size(rwb)),accRatio(size(rwb))
    real(r8) uold(size(rwb)), x(ne), y(ne), z(ne)
    real(r8) dw
    logical isAccepted(size(rwb)), isMoved(size(rwb)), noTmove(size(rwb))
    type(eConfigArray) :: ec

    call assert(size(mRWSave)>=size(rwb), '(propagateAndWeightAll): walker array size too small')

    walkerBlockSize = size(rwb)
    do i = 1, walkerBlockSize
       mRWSave(i) = rwb(i)
    enddo

    ! propagation step
       select case (mMove)
       case (1)
          call reyProp(rwb, accRatio, r2, noTmove)
       case (2)
          call umrProp(rwb, accRatio, r2, vvr, vvnr)
       case (3)
          call twoLevelProp(rwb, accRatio, r2)
       case (4)
          call gaussProp(rwb, accRatio, r2)
       case default
          call abortp("propagateAndWeight: illegal mMove")
       end select

    ! For vmc Sample creation
      mAccRatio = accRatio(1)

    ! acceptance rejection step
    if (mARStep .and. mode==0 .or. mode==1) then
       accepted = 0.d0
       do i = 1, walkerBlockSize
          xi = myran()
          if (accRatio(i) > xi) then
             isAccepted(i) = .true.
             accepted = accepted + 1.d0/walkerBlockSize
             sR2acc = sR2acc + r2(i)
             call resetPersistency(rwb(i))
             if (logmode>=3) then
                mCounter = mCounter + 1
                !!!call eloc_getECPContribs(EL,EECPL,EECPNL)
                !!!write(iull,'(i6,3g18.6)') mCounter, EL, EECPL, EECPNL
             end if
             if(mTmove%mode > 0 .and. mTMove%tmoveddeloc) call resetTrw(rwb(i))
          else
             isAccepted(i) = .false.
             rwb(i) = mRWSave(i)              ! go back to old walker
             call incrPersistency(rwb(i))
              if ( mTmove%mode > 0 .and. mTmove%tmoveddeloc .and. istmoved(rwb(i))) then
                call elocalSet(rwb(i),E_localSaveGet(rwb(i)),.false.) ! local energy for rwb(i) is last step tmoved and for mRWSave(i) is last difussion
                call resetPersistency(rwb(i))
              endif
             ! about the last condition on next line Tmove should not be done after rejection if the reson for rejection was tmove crossing
             if (mTmove%mode > 0 .and. mTMove%tmoverej == 0 .and. .not.noTmove(i)) then
                if (size(rwb) > 1) call abortp("propagateAndWeightAll: propagation with T move only with one parallel walker")
                tmove = mTmove%mode
                if (mWeight==1) then !for rey weight
                    tau=sTauAcc
                  else
                    tau = mTau
                endif
                call eConfigArray_new(ec,ne,1)
                call pos(rwb(i),x,y,z)
                call eConfigArray_set(ec,i,x,y,z)
                call eloc(0, ec, 'none', tmove=tmove, tau=tau, isMoved=isMoved)
                ! calculate Eloc/drift at T move position if moved
                if (isMoved(1)) call rw_setToBlock(mRWTMoves(1:1), ec)
                call eConfigArray_destroy(ec)
                mIsMoved(1) = isMoved(1)
             end if
          end if
          sR2all = sR2all + r2(i)
       end do
    else
       if (mode==2) then
          accepted = 1.d0
       endif
    endif

    !if(mMove == 3) call rw_recalculateElocBlock(rwb,3,isAccepted)

    ! T Moves weight with E_loc at T Move position
    if (mTmove%mode > 0 .and. mTMove%tmovew ==1 ) then
      call internal_rwtotmove()
    end if



    ! weighting step
    select case (mWeight)
    case (0)   ! no weighting
       continue
    case (1)   ! 'Reynolds'
       do i=1,walkerBlockSize
          s1 = E_local(mRWSave(i)) - ERef
          sn1 = E_local(rwb(i)) - ERef
          dw = exp(-0.5d0*sTauAcc*(s1+sn1))
          call multiplyWeight(rwb(i), dw)
       enddo
    case (2)   ! 'Umrigar'
       do i=1,walkerBlockSize
          mp = min(1.d0,accRatio(i))   ! Metropolis probability
          s1 =  (E_local(mRWSave(i)) - sERef0)*sqrt(vvr(i)) - (ERef-sERef0)
          sn1 = (E_local(rwb(i)) - sERef0)*sqrt(vvnr(i)) - (ERef-sERef0)
          dw = exp(-sTauAcc*(0.5d0*mp*(s1+sn1)+(1.d0-mp)*s1))
          call multiplyWeight(rwb(i), dw)
       enddo
    case (3)   ! 'Accept'
       ! reweight accepted moves only (with tau instead of taua)
       do i=1,walkerBlockSize
          if (isAccepted(i)) then
             s1 = E_local(mRWSave(i)) - ERef
             sn1 = E_local(rwb(i)) - ERef
             dw = exp(-0.5d0*mTau*(s1+sn1))
             call multiplyWeight(rwb(i), dw)
          endif
       enddo
    case (4)   ! 'Tmove non symetric'
       do i=1,walkerBlockSize
            dw = exp(-mtau*(E_local(rwb(i)) - ERef))
            call multiplyWeight(rwb(i), dw)
       enddo
    case (5)   ! 'Tmove symetric'
       do i=1,walkerBlockSize
          s1 = E_local(mRWSave(i)) - ERef
          sn1 = E_local(rwb(i)) - ERef
          dw = exp(-0.5d0*mTau*(s1+sn1))
          call multiplyWeight(rwb(i), dw)
       enddo
    end select

    !!if (logmode >= 3) then
    !!  call assert(walkerBlockSize==1,"propagator: statistics only for walker blocksize == 1")
    !!  call addData(mDBGWeights,dw)
    !!  call addData(mDBGELoc,E_local(rwb(1)))
    !!  call addData(mDBGERef,ERef)
    !!end if

    if (mTmove%mode > 0 .and. mTMove%tmovew==0 ) then
      call internal_rwtotmove()
    end if
    ! T Moves, with E_loc/drift at T Move position, weight with E_loc at drift-diffusion pos.

   contains
      subroutine internal_rwtotmove()

       select case (mTMove%tmoverej)
       case (1)
          do i = 1, walkerBlockSize
             if (isAccepted(i) .and. mIsMoved(i)) then
             call setWeight(mRWTMoves(i),wgt(rwb(i)))
                rwb(i) = mRWTMoves(i)
                call addData(mTMoveCount,1.d0)
             else
                call addData(mTMoveCount,0.d0)
             end if
          end do
       case (0,2)
          do i = 1, walkerBlockSize
             if (mIsMoved(i)) then
             call setWeight(mRWTMoves(i),wgt(rwb(i)))
             if(mTMove%tmoveddeloc) then
               call elocalSave(mRWTMoves(i),E_local(mRWTMoves(i))) ! save tmove local energy in randomwalker rw%elocsave
               call elocalSet(mRWTMoves(i),E_local(rwb(i)),.true.)       ! set the local enegy to drift difusion local energy
             endif ! at this point the local energy of the walker is drift difusion eloc and tmoved local energy is saved in elocsave
                rwb(i) = mRWTMoves(i)
                call addData(mTMoveCount,1.d0)
             else
                call addData(mTMoveCount,0.d0)
             end if
          end do
       case default
          call abortp("propagateAndWeightAll: illegal mTMoveRej")
       end select

       end subroutine internal_rwtotmove

  end subroutine propagateAndWeightAll



  subroutine twoLevelProp(rwb, accRatio, r2)
  !-----------------------------------------!
    type(RandomWalker), intent(inout) :: rwb(:)
    real(r8), intent(out) :: accRatio(:)
    real(r8), intent(out) :: r2(:)

    integer i,w
    real(r8) :: u
    real(r8) :: phiOld(size(rwb)), uOld(size(rwb))
    logical :: calcjs(size(rwb))
    real(r8) :: diffx(ne),diffy(ne),diffz(ne)
    real(r8) :: x(ne),y(ne),z(ne),x1(ne),y1(ne),z1(ne)
    real(r8) :: xNew(ne),yNew(ne),zNew(ne)
    real(r8) :: driftx(ne),drifty(ne),driftz(ne)
    type(eConfigArray) :: ecOld,ecNew

    call eConfigArray_new(ecOld,ne,size(rwb))
    call eConfigArray_new(ecNew,ne,size(rwb))

    if (.not. mARstep) then
      ! force acceptance
      accRatio = 1d0
    else
      do w=1,size(rwb)
        phiOld(w) = phi(rwb(w))
        uOld(w) = ju(rwb(w))
        do i=1,ne
           diffx(i) = mSqrtTau*mygran()
           diffy(i) = mSqrtTau*mygran()
           diffz(i) = mSqrtTau*mygran()
        enddo
        call pos(rwb(w),x,y,z)
        call eConfigArray_set(ecOld,w,x,y,z)

        x1 = x  + diffx
        y1 = y  + diffy
        z1 = z  + diffz
        r2(w)=dot_product(diffx,diffx) + dot_product(diffy,diffy) + dot_product(diffz,diffz)
        call eConfigArray_set(ecNew,w,x1,y1,z1)
      enddo
      call rw_setToBlock(rwb,ecNew,1)
      do w=1,size(rwb)
        accRatio(w) = (phi(rwb(w))/phiOld(w))**2 !* exp(-0.5d0*mp/mTau)
        calcjs(w) = (accRatio(w) > myran())
      enddo

      call rw_setToBlock(rwb,ecNew,3,calcjs)
      do w=1,size(rwb)
        if (calcjs(w)) then
          accRatio(w) = exp(2.d0*(ju(rwb(w))-uOld(w)) )
        else
          accRatio(w) = 0.00D+00
        endif
      enddo
    endif

    call eConfigArray_destroy(ecOld)
    call eConfigArray_destroy(ecNew)
  end subroutine twoLevelProp


  subroutine gaussProp(rwb, accRatio, r2)
  !-----------------------------------------!
    type(RandomWalker), intent(inout) :: rwb(:)
    real(r8), intent(out) :: accRatio(:)
    real(r8), intent(out) :: r2(:)

    integer i,w,tmove
    real(r8) :: u,tau
    real(r8) :: phiOld(size(rwb)), uOld(size(rwb))
    logical :: calcjs(size(rwb))
    real(r8) :: diffx(ne),diffy(ne),diffz(ne)
    real(r8) :: x(ne),y(ne),z(ne),x1(ne),y1(ne),z1(ne)
    real(r8) :: xNew(ne),yNew(ne),zNew(ne)
    real(r8) :: driftx(ne),drifty(ne),driftz(ne)
    type(eConfigArray) :: ecOld,ecNew
    logical :: isMoved(size(rwb))

    call eConfigArray_new(ecOld,ne,size(rwb))
    call eConfigArray_new(ecNew,ne,size(rwb))

    if (.not. mARstep) then
      ! force acceptance
      accRatio = 1d0
    else
      do w=1,size(rwb)
        phiOld(w) = phi(rwb(w))
        uOld(w) = ju(rwb(w))
        do i=1,ne
           diffx(i) = mSqrtTau*mygran()
           diffy(i) = mSqrtTau*mygran()
           diffz(i) = mSqrtTau*mygran()
        enddo
        call pos(rwb(w),x,y,z)
        call eConfigArray_set(ecOld,w,x,y,z)

        x1 = x  + diffx
        y1 = y  + diffy
        z1 = z  + diffz
        r2(w)=dot_product(diffx,diffx) + dot_product(diffy,diffy) + dot_product(diffz,diffz)
        call eConfigArray_set(ecNew,w,x1,y1,z1)
      enddo
      if (mTmove%mode > 0) then
         if (size(rwb) > 1) call abortp("propagation with T move only with one parallel walker")
         tmove = mTmove%mode
         tau = mTau
         call rw_setToBlock(rwb, ecNew, tmove=tmove, tau=tau, isMoved=isMoved)
         ! calculate Eloc/drift at T move position if moved
         if (isMoved(1)) call rw_setToBlock(mRWTMoves(1:1), ecNew)
         mIsMoved(1) = isMoved(1)
      else
         call rw_setToBlock(rwb,ecNew)
      end if
      do w=1,size(rwb)
        accRatio(w) = (phi(rwb(w))/phiOld(w))**2 * exp(2.d0*(ju(rwb(w))-uOld(w)) )
      enddo
    endif

    call eConfigArray_destroy(ecOld)
    call eConfigArray_destroy(ecNew)
  end subroutine gaussProp


  subroutine reyProp(rwb,accRatio,r2,noTmove)
  !---------------------------------!
    ! Reynolds propagator (JCP, 1982 paper)
    ! this version propagates a block of walkers rwb
    type(RandomWalker), intent(inout) :: rwb(:)        ! random walker to propagate
    real(r8), intent(out)               :: accRatio(:)   ! acceptance ratio
    real(r8), intent(out)               :: r2(:)         ! squared diffusion displacement
    logical,intent(out)               :: noTmove(:)

    integer i,w,tmove
    real(r8) :: mp,u,tau
    real(r8) :: phiOld(size(rwb)),uOld(size(rwb))
    real(r8) :: diffx(ne),diffy(ne),diffz(ne)
    real(r8) :: x(ne),y(ne),z(ne),x1(ne),y1(ne),z1(ne)
    real(r8) :: xNew(ne),yNew(ne),zNew(ne)
    real(r8) :: driftx(ne),drifty(ne),driftz(ne)
    type(eConfigArray) :: ecOld,ecNew
    logical :: isMoved(size(rwb)), driftDiffCross, TMoveCross, crossing

    noTmove=.false.
    TMoveCross = .false.
    call assert(size(rwb)>=1,'(reyProp): walker block array size not set')
    call assert(size(rwb)==size(accRatio) .and. size(rwb)==size(r2),'(reyProp): illegal array sizes')

    if (logmode > 5) then
       write(iul,*) 'in reyProp:'
       do i=1,size(rwb)
          call display(rwb(i),2,iul)
       enddo
    endif
    call eConfigArray_new(ecOld,ne,size(rwb))
    call eConfigArray_new(ecNew,ne,size(rwb))

    do w=1,size(rwb)
       phiOld(w) = phi(rwb(w))
       uOld(w) = ju(rwb(w))

       do i=1,ne
          diffx(i) = mSqrtTau*mygran()
          diffy(i) = mSqrtTau*mygran()
          diffz(i) = mSqrtTau*mygran()
       enddo

       call pos(rwb(w),x,y,z)
       call eConfigArray_set(ecOld,w,x,y,z)
       call drift(rwb(w),driftx,drifty,driftz)
       x1 = x + mDriftScal*mTau*driftx + diffx
       y1 = y + mDriftScal*mTau*drifty + diffy
       z1 = z + mDriftScal*mTau*driftz + diffz
       call eConfigArray_set(ecNew,w,x1,y1,z1)

       r2(w) = dot_product(diffx,diffx) + dot_product(diffy,diffy) + dot_product(diffz,diffz)
    enddo

    if (mTmove%mode > 0) then
      if (size(rwb) > 1) call abortp("propagation with T move only with one parallel walker")
      tmove = mTmove%mode
      if (mWeight==1) then !for rey weight
         tau=sTauAcc
       else
         tau = mTau
      endif
      call rw_setToBlock(rwb, ecNew, tmove=tmove, tau=tau, isMoved=isMoved)
      ! calculate Eloc/drift at T move position if moved
      if (isMoved(1)) call rw_setToBlock(mRWTMoves(1:1), ecNew)
      mIsMoved(1) = isMoved(1)
    else
      call rw_setToBlock(rwb,ecNew)
    end if

    ! calculate acceptance ratio 'accRatio' only if necessary
    !  force acceptance with accRatio = 1d0, rejection with = 0d0
    if (.not. mARstep) then                                ! force acceptance
       accRatio = 1d0
    else
       do w=1,size(rwb)
          driftDiffCross = (phi(rwb(w))*phiOld(w) < 0.d0)
          TMoveCross = mIsMoved(1) .and. (phi(rwb(w))*phi(mRWTMoves(1)) < 0.d0)
          select case (mTMove%tmovecross)
          case (0)
            crossing = driftDiffCross
          case (1)
            crossing = driftDiffCross .or. TMoveCross
          case (2)
            crossing = driftDiffCross .neqv. TMoveCross
          end select
          ! ! Tmove should not be done after rejection if the reson for rejection was tmove crossing
          if (crossing .and. .not.driftDiffCross) noTmove(w)=.true.

          !!!if (logmode >= 3) print*, crossing, driftDiffCross, TMoveCross, mIsMoved(1)
          !!!if (phi(rwb(w))*phiOld(w) < 0.d0 .and. mRejCross) then  ! reject node cross
          if (crossing .and. mRejCross) then
             accRatio(w) = 0d0
             if (logmode >= 5) write(iul,*) 'attempted xing from walker ',w,' of ',size(rwb)
          else                                                  ! calculate accRatio
             call pos(rwb(w),xNew,yNew,zNew)
             call eConfigArray_get(ecOld,w,x,y,z)
             call drift(rwb(w),driftx,drifty,driftz)
             x1 = x - xNew - mDriftScal*mTau*driftx
             y1 = y - yNew - mDriftScal*mTau*drifty
             z1 = z - zNew - mDriftScal*mTau*driftz
             mp = dot_product(x1,x1) + dot_product(y1,y1) + dot_product(z1,z1) - r2(w)
             !accRatio = (psi(rw)/psiOld)**2 * exp(-0.5d0*mp/mTau)
             accRatio(w) = (phi(rwb(w))/phiOld(w))**2 * exp(-0.5d0*mp/mTau + 2.d0*(ju(rwb(w))-uold(w)) )
          endif
       enddo
    endif

    call eConfigArray_destroy(ecOld)
    call eConfigArray_destroy(ecNew)
  end subroutine reyProp


  subroutine umrProp(rwb,accRatio,r20,vvr,vvnr)
  !-------------------------------------------!
    ! Move according to Umrigar/Nightingale/Runge paper (JCP, 1993)
    ! velocity ratios are returned for the weighting step
    type(RandomWalker), intent(inout) :: rwb(:)         ! random walker block to propagate
    real(r8), intent(out)               :: accRatio(:)    ! acceptance ratio
    real(r8), intent(out)               :: r20(:)         ! squared diffusion displacement
    real(r8), intent(out)               :: vvr(:)         ! velocity ratio
    real(r8), intent(out)               :: vvnr(:)        ! new velocity ratio

    integer                           :: i,nnu0
    real(r8)                            :: x(ne),y(ne),z(ne)
    real(r8)                            :: xOld(ne),yOld(ne),zOld(ne)
    real(r8)                            :: driftx(ne),drifty(ne),driftz(ne)
    real(r8)                            :: phiOld(size(rwb)),uOld(size(rwb))
    real(r8)                            :: rnnuc(ne)
    integer                           :: nnuc(ne)
    real(r8)                            :: vv,vv1,vvn,vvn1
    real(r8)                            :: tau,chi,dx,dy,dz,tmp,tmp1,t2
    real(r8)                            :: vf,v2,v20,vz,vxz,vrhox,vrhoy,vrhoz,zz1,rho,z11,rho11
    real(r8)                            :: ddx,ddy,ddz,vv1x,vv1y,vv1z,p1,q1,gf,gf1
    real(r8)                            :: zeta,rr,costh,sinth,phi1,chi1,chi2,chi3
    real(r8)                            :: erfc
    integer                           :: w,tmove
    type(eConfigArray) :: ecOld,ecNew
    logical :: isMoved(size(rwb))

    call assert(size(rwb)>=1,'(umrProp): walker block array size not set')
    call assert(size(rwb)==size(accRatio) .and. size(rwb)==size(r20),'(umrProp): illegal array sizes')
    call assert(size(rwb)==size(vvr) .and. size(rwb)==size(vvnr),'(umrProp): illegal vvr sizes')
    call eConfigArray_new(ecOld,ne,size(rwb))
    call eConfigArray_new(ecNew,ne,size(rwb))

    t2 = mSqrtTau
    do w=1,size(rwb)
       call pos(rwb(w),x,y,z)
       call eConfigArray_set(ecOld,w,x,y,z)
       phiOld(w) = phi(rwb(w))
       uOld(w)   = ju(rwb(w))
    enddo

    do w=1,size(rwb)
      call eConfigArray_get(ecOld,w,x,y,z)
      call getNearestNucleusDist(rwb(w),rnnuc,nnuc)
      call drift(rwb(w),driftx,drifty,driftz)

      r20(w) = 0d0                         ! to sum squared diffusion displacement
      vv = 0d0                          ! sum unchanged velocities (abs.value)
      vv1 = 0d0                         ! sum modified velocities
      gf  = 1d0                         ! Green's function
      gf1 = 1d0                         ! Green's function for reversed move

      ! move all electrons
      do i=1,ne
         if (logmode >= 5) write(iul,*) 'el. ',i
         nnu0 = nnuc(i)
         zz1 = rnnuc(i)

         zeta = sqrt(atoms(nnu0)%za**2 + 1d0/mTau)

         v2  = driftx(i)**2 + drifty(i)**2 + driftz(i)**2
         v20 = v2
         vv  = vv + v20

         ! Calculate drift
         tmp = (atoms(nnu0)%za*zz1)**2
          vxz = (driftx(i)*(x(i)-atoms(nnu0)%cx)+drifty(i)*(y(i)-atoms(nnu0)%cy) &
              +  driftz(i)*(z(i)-atoms(nnu0)%cz))/(zz1*sqrt(v2))


         v2  = v2 * (0.5d0*(1d0 + vxz) + tmp/((4d0 + tmp)*10d0))
         vf  = (-1d0 + sqrt(1d0 + 2d0*v2*mTau))/(v2*mTau)
         if (logmode .ge. 5) write(iul,'(A6,2G12.5)') &
             'vf,a=',vf,tmp/(4d0+tmp)
         vv1x = vf*driftx(i)                     ! modified velocity
         vv1y = vf*drifty(i)
         vv1z = vf*driftz(i)

         vz = (vv1x*(x(i)-atoms(nnu0)%cx) &              ! z-component
             + vv1y*(y(i)-atoms(nnu0)%cy) &
             + vv1z*(z(i)-atoms(nnu0)%cz))/zz1
         vrhox = vv1x - vz/zz1*(x(i)-atoms(nnu0)%cx)
         vrhoy = vv1y - vz/zz1*(y(i)-atoms(nnu0)%cy)
         vrhoz = vv1z - vz/zz1*(z(i)-atoms(nnu0)%cz)
         rho   = sqrt(vrhox**2 + vrhoy**2 + vrhoz**2)
         z11   = max(zz1 + vz*mTau,0d0)
         rho11 = 2d0*rho*mTau*z11 / (zz1 + z11)

         if (logmode >= 5) write(iul,'(A22,5G10.3)')  &
             'zz1,vz,z11,rho,rho11:',zz1,vz,z11,rho,rho11

         ! new position after drift
         ddx = atoms(nnu0)%cx + rho11/rho*vrhox  &
              + z11/zz1*(x(i)-atoms(nnu0)%cx)
         ddy = atoms(nnu0)%cy + rho11/rho*vrhoy  &
              + z11/zz1*(y(i)-atoms(nnu0)%cy)
         ddz = atoms(nnu0)%cz + rho11/rho*vrhoz  &
              + z11/zz1*(z(i)-atoms(nnu0)%cz)

         if (logmode >= 5) write(iul,'(a4,3g12.5)') 'dd:',ddx,ddy,ddz

         ! Calc prob. if sampling from Exponential or Gaussian
         if (mNoExp) then                         ! don't sample exponential
            q1  = 0d0
            p1  = 1d0
            chi = .5d0                            ! don't call mrg, force gaussian
         else
            q1 = 0.5d0*erfc((zz1+vz*mTau)/sqrt(2d0*mTau))
            p1 = 1d0 - q1
            chi = myran()
         endif
         if (chi > q1) then                       ! sample from Gaussian
            if (logmode >= 5) write(iul,'(A13,2G12.5)') 'gauss:q1,chi',q1,chi
            chi = mygran()
            dx = t2*chi
            x(i) = ddx + dx
            chi = mygran()
            dy = t2*chi
            y(i) = ddy + dy
            chi = mygran()
            dz = t2*chi
            z(i) = ddz + dz
            r20 = r20 + dx**2 + dy**2 + dz**2
         else                                      ! Sample from Exponential
            if (logmode >= 5) write(iul,'(A13,2G12.5)')  'expon:q1,chi',q1,chi
            chi1 = myran()
            chi2 = myran()
            chi3 = myran()
            rr = - log(chi1*chi2*chi3) / (2d0*zeta)
            chi = myran()
            costh = 1d0 - 2d0*chi
            sinth = sqrt(1d0 - costh*costh)
            chi = myran()
            phi1 = 2d0*pi*chi
            dx = rr*sinth*cos(phi1)
            x(i) = atoms(nnu0)%cx + dx
            dy = rr*sinth*sin(phi1)
            y(i) = atoms(nnu0)%cy + dy
            dz = rr*costh
            z(i) = atoms(nnu0)%cz + dz
            r20 = r20 + dx**2 + dy**2 + dz**2
         endif
         tmp = (x(i)-ddx)**2 + (y(i)-ddy)**2  &
              +    (z(i)-ddz)**2
         tmp = exp(-tmp/(2d0*mTau))/(2d0*pi*mTau)**1.5d0
         tmp1 = sqrt((x(i)-atoms(nnu0)%cx)**2       &
              + (y(i)-atoms(nnu0)%cy)**2 + (z(i)-atoms(nnu0)%cz)**2)
         tmp1 = zeta**3/pi*exp(-2d0*zeta*tmp1)
         gf = gf * (p1*tmp + q1*tmp1)
         if (logmode .ge. 5) write(iul,'(A8,5G24.16)') 'g,g1,g2', &
              p1*tmp + q1*tmp1,tmp,tmp1,p1,q1
         vv1 = vv1 + vf*v20
      enddo

      vvr(w) = vv1/vv
      call eConfigArray_set(ecNew,w,x,y,z)
    enddo

    ! calculate psi and E_local for all walkers in rwb at new positions

    if (mTmove%mode > 0) then
      if (size(rwb) > 1) call abortp("propagation with T move only with one parallel walker")
      tmove = mTmove%mode
      tau = mTau
      call rw_setToBlock(rwb, ecNew, tmove=tmove, tau=tau, isMoved=isMoved)
      ! calculate Eloc/drift at T move position if moved
      if (isMoved(1)) call rw_setToBlock(mRWTMoves(1:1), ecNew)
      mIsMoved(1) = isMoved(1)
    else
      call rw_setToBlock(rwb,ecNew)
    end if

    ! calculate terms for weighting and acceptance ratio

    do w=1,size(rwb)

      call eConfigArray_get(ecNew,w,x,y,z)
      call eConfigArray_get(ecOld,w,xOld,yOld,zOld)
      call getNearestNucleusDist(rwb(w),rnnuc,nnuc)
      call drift(rwb(w),driftx,drifty,driftz)

      ! Check node crossing, calculate Metropolis probability accratio

      ! Calculate mod. velocity for new position (as above)
      !!!   (Instead of Calculating twice, better use vf(i))
      vvn = 0d0              ! sum unchanged new velocities (abs.value)
      vvn1 = 0d0             ! sum modified new velocities
      do i=1,ne
         if (logmode >= 5) write(iul,*) 'el. ',i
         nnu0 = nnuc(i)
         zz1 = rnnuc(i)

         zeta = sqrt(atoms(nnu0)%za**2 + 1d0/mTau)

         v2 = driftx(i)**2 + drifty(i)**2 + driftz(i)**2
         v20 = v2
         vvn = vvn + v20

         ! Calculate drift (with modified a-factor)
         tmp = (atoms(nnu0)%za*zz1)**2
         vxz = (driftx(i)*(x(i)-atoms(nnu0)%cx)+drifty(i)*(y(i)-atoms(nnu0)%cy) &
             +  driftz(i)*(z(i)-atoms(nnu0)%cz))/(zz1*sqrt(v2))
         v2  = v2 * (0.5d0*(1d0+ vxz) + tmp/((4d0 + tmp)*10d0))
         vf  = (-1d0 + sqrt(1d0 + 2d0*v2*mTau))/(v2*mTau)
         if (logmode >= 5) write(iul,'(A6,2G12.5)') 'vf,a=',vf,tmp/(4d0+tmp)
         vv1x = vf*driftx(i)               ! modified velocity
         vv1y = vf*drifty(i)
         vv1z = vf*driftz(i)

         vz = (vv1x*(x(i)-atoms(nnu0)%cx)  &      ! z-component
              + vv1y*(y(i)-atoms(nnu0)%cy) &
              + vv1z*(z(i)-atoms(nnu0)%cz))/zz1
         vrhox = vv1x - vz/zz1*(x(i)-atoms(nnu0)%cx)
         vrhoy = vv1y - vz/zz1*(y(i)-atoms(nnu0)%cy)
         vrhoz = vv1z - vz/zz1*(z(i)-atoms(nnu0)%cz)
         rho   = sqrt(vrhox**2 + vrhoy**2 + vrhoz**2)
         z11   = max(zz1 + vz*mTau,0d0)
         rho11 = 2d0*rho*mTau*z11 / (zz1 + z11)

         ddx = atoms(nnu0)%cx + rho11/rho*vrhox  &
              + z11/zz1*(x(i)-atoms(nnu0)%cx)
         ddy = atoms(nnu0)%cy + rho11/rho*vrhoy  &
              + z11/zz1*(y(i)-atoms(nnu0)%cy)
         ddz = atoms(nnu0)%cz + rho11/rho*vrhoz  &
              + z11/zz1*(z(i)-atoms(nnu0)%cz)

         if (mNoExp) then                     ! no exponential sampling
            q1 = 0d0
            p1 = 1d0
         else
            q1 = 0.5d0*erfc((zz1+vz*mTau)/sqrt(2d0*mTau))
            p1 = 1d0 - q1
         endif
         tmp = (xOld(i)-ddx)**2 + (yOld(i)-ddy)**2  &
              +    (zOld(i)-ddz)**2
         tmp =  exp(-tmp/(2d0*mTau))/(2d0*pi*mTau)**1.5d0
         tmp1 = sqrt((xOld(i)-atoms(nnu0)%cx)**2  &
              + (yOld(i)-atoms(nnu0)%cy)**2 + (zOld(i)-atoms(nnu0)%cz)**2)
         tmp1 = zeta**3/pi*exp(-2d0*zeta*tmp1)
         gf1 = gf1 * (p1*tmp + q1*tmp1)
         if (logmode >= 5) write(iul,'(A8,3G12.5)') 'g,g1,g2',  &
              p1*tmp + q1*tmp1,tmp,tmp1

         vvn1 = vvn1 + vf*v20
      enddo

      if (.not. mARstep) then
         accratio(w) = 1d0
      else if (phi(rwb(w))*phiOld(w) < 0.d0 .and. mRejCross) then     ! node crossing
         accratio(w) = 0d0
         if (logmode >= 5) write(iul,*) 'attempted xing'
      else
         accratio(w) = (phi(rwb(w))/phiOld(w))**2 * exp( 2.d0*(ju(rwb(w))-uold(w)) )  * gf1/gf
      endif

      if (logmode >= 5) write(iul,*) ' Metrop.prob=',accratio(w),gf1,gf  &
                                                    ,phi(rwb(w)),phiOld(w)

      vvnr(w) = vvn1 / vvn
    enddo

    call eConfigArray_destroy(ecNew)
    call eConfigArray_destroy(ecOld)

  end subroutine umrProp


  !-----------------------------------------------------!
  subroutine propagateAndWeightOne(rw,ERef,mode,accepted)
  !-----------------------------------------------------!

    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: ERef
    integer, intent(in)               :: mode        ! 0: default (i.e. use mARStep)
                                                     ! 1,2: with,without accept step
    real(r8), intent(out)               :: accepted    ! ratio of accepted steps

    ! dummy routine
    accepted = 0

  end subroutine propagateAndWeightOne

end module propagator_m



