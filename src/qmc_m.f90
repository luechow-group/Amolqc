! Copyright (C) 1996, 1999, 2012-2015, 2018 Arne Luechow
! Copyright (C) 2000 Sebastian Manten
! Copyright (C) 2012-2013 Alexander Sturm
! Copyright (C) 2013, 2015-2016, 2018 Kaveh Haghighi Mood
! Copyright (C) 2018-2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

MODULE qmc_m

! containing qmc parameters and qmc run routines
! correction in Branching -> NSplit ps=0
! correction of bgte in qmc_m.f90 and propagator_m.f90 included

  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value, ieee_is_normal
  use kinds_m, only: r8, i8
  use OMP_LIB
  use global_m
  use wfData_m, only: do_epart, vpot0
  use properties_m
  use statistics
  use newStatistics_m, only: stat,vectorstat
  use randomWalker_m
  use rwSample_m
  use propagator_m
  use rwStatistics_m
  use parsing_m
  use utils_m, only: getTimeString
  use elocData_m, only: setElocCutoff, getElocCutoff, eloc_initialize, &
         mElocCut, mElocCutCount, mDriftCut, mDriftCutCount
  use eloc_m
  use aosData_m, only: aos_initialize
  use mos_m, only: mos_initialize
  use waveFunction_m, only: trajectory_head
  use reconfg_m
  use psiMax_m, only: psimax, RAWDATA, NOANALYSIS
  use rhoMax_m, only: rhoMax_t
  use energyPart_m
!  use assignment
  use assign_m
  use qmcSample_m, only: currentSampleSize
  use random_m, only: myran

  ! uses for on the fly blocking statistics
  use decorrelation_m, only: DecorrelationData_t, Decorrelation_t
#ifdef MPI
  use parallelDDDA_m, only: myMPIReduceDecorr
#endif
  use kinds_m, only: i8

!  use ce_new_m
  implicit none

  private
  public  :: qmc_init, initWalkerStat, initTrajectory, qmc_writeParams, qmc_run, qmc_Energy, qmc_StdDev, qmc_Variance

  integer, parameter :: WEIGHT_NONE=0, WEIGHT_REY=1, WEIGHT_UMR=2, WEIGHT_ACCEPT=3, WEIGHT_TMSM=4, WEIGHT_TMSC=5
  integer, parameter :: POP_NONE=0, POP_GLOBAL=1, POP_LOCAL=2

  ! parameters (constant during run)
  real(r8)      :: mWFac=1.d0                ! scale factor in Umrigar pop control
  real(r8)      :: mLThresh=0.5d0, mUThresh=2.0d0  ! thresholds for branching
  real(r8)      :: mTargetAR=0.5d0           ! target acceptance ratio (AR) for AR adaptation
  real(r8)      :: mEpsAR=0.01d0             ! required accuracy for AR
  real(r8)      :: mStdDev=0.05d0            ! target std deviation for result
  real(r8)      :: mMaxTau=1._r8             ! max value for time step adaptation
  integer(i8)   :: mSteps=1000               !
  integer(i8)   :: mBlockLen=100             ! blocking of steps
  integer(i8)   :: mStepsDiscard=0           ! steps to discard initially
  integer(i8)   :: mStepStride=10            ! stride for data collection
  integer     :: mTargetSampleSize=0             ! actual size for VMC, target size for DMC, this node
  integer     :: mTargetSampleSizeAllNodes=0     ! same for all nodes
  integer     :: mPerMax=10                ! max allowed persistency
  integer     :: mWeight=WEIGHT_NONE       ! weighting strategy (none=0,Rey=1,Umr,Accept)
  integer     :: mPopctrl=POP_GLOBAL       ! 0: noPopctrl,=1 : Global Popctrl in MPIJob, =2 : separated/local for each node
  integer     :: mRcf=1                    ! 1: Our Method 2: Implementation following Caffarel
  integer     :: mWalkerBlock=1            ! # of walkers to be moved and calculated "simultaneously"
  integer     :: mAutoCorrMax=500          ! calc autocorrelation up to mAutoCorrMax steps
  integer     :: mShowDetails=0            ! details in output
  logical     :: mCheckStdDev=.false.      ! check if target std dev is reached
  logical     :: mBranch=.false.
  logical     :: mReconf=.false.
  logical     :: mJoin=.true.
  logical     :: mKill=.false.
  logical     :: mShowSteps=.false.
  logical     :: mSplit=.true.
  logical     :: mKillPersist=.false.
  logical     :: mWalkerStatistics = .false.  ! do topological analysis for walkers?
  logical     :: mTrajectory = .false.        ! record trajectories of walkers?
  logical     :: mLoadBalance = .false.       ! balance sample sizes in parallel runs
  logical     :: mAdaptTau = .false.          ! adapt tau at beginning of run to some AR
  logical     :: mAutoCorrelation = .false.   ! calc autocorrelation of E_L data (only VMC)
  character(len=1) :: mStatType='c'
  character(len=3) :: mMethod='VMC'        ! QMC method

  !Future Walking
  logical    :: mFuture = .false.            ! switch for Future-Walking
  integer    :: mDprops = 0                  ! number of Generations (saved fathers managed)
  integer    :: mStprops = 0                 ! Offset between fathers taken (tau(corr) should be taken into account)
                                              ! stprops*dprops should equal Block length
  integer    :: mTau2   = 0                 ! Steps between output pro Generation (tau itself is like die usual tau)
  integer    :: mNotau2 = 0                  ! no. of outputs pro stprops*dprops

  ! status parameters
  integer :: sBlock=0
  integer :: sSampleSize=0
  integer :: sOldLogmode=0
  real(r8)  :: sERef=0
  integer :: sMode=0            ! propagator mode
  real(r8)  :: sEMean=0           ! QMC result: Energy
  real(r8)  :: sEMeanStdDev = 0   ! std deviation (of mean)
  real(r8)  :: sVar=0             ! variance (of local energy)
  real(r8)  :: sVarError=0._r8 ! error of the variance
  type(simpleStat),save :: ElocBlockStat
  type(simpleStat),save :: varBlockStat
  type(simpleStat),save :: VenBlockStat
  type(simpleStat),save :: VeeBlockStat
  type(simpleStat),save :: TBlockStat
  type(weightStat),save :: sTotalStat
  !Assignment parameters
  logical             :: mAssignment=.false.
  ! Parameters for Maxima Search
  logical            :: mMaxSearch = .false.
  logical            :: mMaxAnalysis = .false.
  integer            :: mBlockOffset = 1
  real(r8),allocatable :: x_max(:),y_max(:),z_max(:)
  real(r8)             :: best
  ! accumulate seen samples into a history
  logical            :: mAccumulate = .false.
  integer            :: mAccumulateSampleSize = 0

contains


  subroutine qmc_init(lines, nl, sample, psimax_obj)
  !------------------------------------------------!

    character(len=120), intent(in)        :: lines(:)
    integer, intent(in)                   :: nl
    type(RWSample), intent(inout)         :: sample
    type(psimax), intent(inout), optional :: psimax_obj


    integer iflag,i
    character(len=3)           :: s,mt
    character(len=12)          :: s1
    logical                    :: found,found1,found2
    real(r8) tau,ds,cf
    integer bgte,move,tauFlag,iflag1,vb,blockdiscard
    integer(i8) accSize
    logical ar,rc,eLocalCutOff
    integer :: int_rcv(nproc), ierr


    type(TMoveControl) :: tmove

    if (present(psimax_obj)) then
       mMaxSearch = psimax_obj%do_max_search()        !!! no longer necessary: data in psimax_obj
       mMaxAnalysis = psimax_obj%do_max_analysis()
    end if

    ! set default values for standard methods
    found1 = finda(lines,nl,'vmc')
    found2 = finda(lines,nl,'VMC')
    if (found1 .or. found2) then
       ! set VMC default values)
       mMethod='VMC';ar=.true.;rc=.false.; move=2; mWeight=WEIGHT_NONE; mBranch=.false.; mKillPersist=.false.
       mPerMax=0; mLoadBalance=.false.;eLocalCutOff=.false.
       mBlockLen=200; mStepsDiscard=0; mAutoCorrelation=.true.
       mTargetAR=0.5d0
    endif

    found1 = finda(lines,nl,'dmc')
    found2 = finda(lines,nl,'DMC')
    mFuture = finda(lines,nl,'FWDMC')
    if (found1 .or. found2 .or. mFuture) then
       ! DMC calculation
       mMethod='DMC';ar=.true.;rc=.true.; move=1; mWeight=WEIGHT_REY; mBranch=.true.; mJoin=.true.; mKill=.false.
       mSplit=.true.;mKillPersist=.true.;mPerMax=10;eLocalCutOff=.true.
       mLoadBalance=.false.; if (nproc>1) mLoadBalance=.true.
       mBlockLen=1000; mStepsDiscard=0; mAutoCorrelation=.false.
       mTargetAR=0.90d0
    endif

    ! QMC parameters
    mTargetSampleSize = 0   ! default value (0) is: "keep current sample size"
    call getinta(lines,nl,'walker=',mTargetSampleSize,iflag)
    mTargetSampleSizeAllNodes = mTargetSampleSize * nproc

    ! propagate walkers (parallel) in small blocks
    call getinta(lines,nl,'walker_block=',mWalkerBlock,iflag)
    if (iflag == 0) then
       call eloc_initialize(mWalkerBlock)
       call aos_initialize(mWalkerBlock)
       call mos_initialize(mWalkerBlock)
    endif

    ! total # of QMC steps
    call getint8a(lines,nl,'steps=',mSteps,iflag)
    mShowSteps = finda(lines,nl,'show_steps')

    ! block length for statistics of correlated data
    call getint8a(lines,nl,'block_len=',mBlockLen,iflag)

    ! allow finishing a run if a given std dev is reached
    call getdbla(lines,nl,'std_dev=',mStdDev,iflag)
    if (iflag == 0) then
       mCheckStdDev = .true.
    else
       mCheckStdDev = .false.
    endif

    ! step_stride for all calls using current walkers
    call getint8a(lines,nl,'step_stride=',mStepStride,iflag)

    ! discard first steps in all statistics
    call getint8a(lines,nl,'discard=',mStepsDiscard,iflag)
    found = finda(lines,nl,'discard_all')
    if (found) mStepsDiscard = mSteps

    if (mMethod == 'DMC' .and. iflag /= 0 .and. .not. found) then
       call error('$qmc: discard required for DMC runs')
    end if

    ! accumulation of large samples with size acc_size for optimization
    mAccumulate = finda(lines,nl,'accumulate')
    if (mAccumulate .and. mMethod /= 'VMC') call abortp('$qmc: accumulate requires a VMC run')
    call getint8a(lines,nl,'acc_size=',accSize,iflag)
    if (iflag==0) mSteps = (accSize - 1) * mStepStride / getSampleSize(sample) + mStepsDiscard

    ! allow to kill persistent walkers after mPerMax steps
    call getinta(lines,nl,'persist=',mPerMax,iflag)
    mKillPersist = .false.
    if (mPerMax > 0)  mKillPersist = .true.

    ! DMC branching
    call getloga(lines,nl,'branch=',mBranch,iflag)
    call getdbla(lines,nl,'lower_thresh=',mLThresh,iflag)
    call getdbla(lines,nl,'upper_thresh=',mUThresh,iflag)
    call getloga(lines,nl,'join=',mJoin,iflag)
    call getloga(lines,nl,'kill=',mKill,iflag)
    if (mKill .and. mJoin) call abortp('initQMCParams: join and kill exclude each other')
    call getloga(lines,nl,'split=',mSplit,iflag)

    ! DMC load balancing
    found = finda(lines,nl,'no_load_balance')
    if (found) then
       mLoadBalance = .false.
    else
       found = finda(lines,nl,'load_balance')
       if (found) mLoadBalance = .true.
    endif

    ! propagator parameters
    call getstra(lines,nl,'move=',s,iflag)
    if (iflag == 0) then
       select case (s)
          case ('Rey', 'rey')
             move = 1
          case ('Umr', 'umr')
             move = 2
          case ('Two', 'two')
             move = 3
          case ('Gss', 'gss')
             move = 4
          case default
             call abortp("$qmc: unknown move given")
       end select
    end if

    call getstra(lines,nl,'weight=',s,iflag)   ! weight overrides default values
    if (iflag == 0) then
       if (s=='none') then
          mWeight = WEIGHT_NONE
       else if (s=='Rey' .or. s=='rey') then
          mWeight = WEIGHT_REY
       else if (s=='Umr' .or. s=='umr') then
          mWeight = WEIGHT_UMR
       else if (s=='Acc' .or. s=='acc') then
          mWeight = WEIGHT_ACCEPT
       else if (s=='Tsm' .or. s=='tsm') then
          mWeight = WEIGHT_TMSM
       else if (s=='Tsc' .or. s=='tsc') then
          mWeight = WEIGHT_TMSC
       else
          call abortp('$qmc: unknown weight given')
       endif
    endif
    if (move==2 .and. .not.(mWeight==WEIGHT_UMR .or. mWeight==WEIGHT_NONE) .or. mWeight==WEIGHT_UMR .and. .not.move==2) then
       call abortp("$qmc: Umrigar move only with Umrigar weight allowed")
    end if

    ! T moves (only with ecps)
    tmove%mode = 0
    call getstra(lines,nl,'T_moves=',s1,iflag)
    if (iflag == 0) then
       if (s1=='simple') tmove%mode = 1
       if (s1=='SIMPLE') then !reasonable default valuess
         tmove%mode = 1
         tmove%tmoveddeloc=.true.
         mShowDetails = 1
       end if
       if (s1=='sc') tmove%mode = 2
       if (s1=='sc1') tmove%mode = 3
       if (s1=='SC') then !reasonable default valuess
         tmove%mode = 3
         tmove%tmoveddeloc=.true.
         mShowDetails = 1
       endif
    end if

    if (tmove%mode > 0) then
       tmove%tmoverej = 0
       call getinta(lines,nl,'T_move_reject=',tmove%tmoverej,iflag)
        if(tmove%tmoverej>2)    call abortp("illegal value for T_move_reject given")
       tmove%tmovew = 0
       call getinta(lines,nl,'T_move_wgt=',tmove%tmovew,iflag)
       if (iflag == 0 .and. tmove%tmoverej > 0 .and. tmove%tmovew>1 ) call abortp("$qmc: T_move_weight=1 only for T_move_reject=0")
       tmove%tmovecross=0
       call getinta(lines,nl,'T_move_cross=',tmove%tmovecross,iflag)
       found = finda(lines,nl,'T_move_dd_eloc')
       if(found) tmove%tmoveddeloc=.true.
    end if

    if ((mWeight==WEIGHT_TMSM .or. mWeight==WEIGHT_TMSC )  .and. tmove%mode < 1  ) then
       call abortp("$qmc: tmo is only valid for  tmove")
    end if

    call getdbla(lines,nl,'wfac=',mWFac,iflag)
    E_trial = 0.d0
    call getdbla(lines,nl,'eref=',E_trial,iflag)
    call getdbla(lines,nl,'E_ref=',E_trial,iflag1)
    if (iflag /= 0 .and. iflag1 /= 0 .and. mWeight /= WEIGHT_NONE) then
       E_trial = qmc_Energy()
       if (E_trial==0) call abortp('(qmc_init): E_ref or previous vmc calc for dmc required')
    end if
    mt = "all"
    call getstra(lines,nl,'move_typ=',mt,iflag)
    call getinta(lines,nl,'blocks_tau_eff=',bgte,iflag)
    if (iflag /= 0) bgte = -1

    ! set or guess time step tau
    mAdaptTau = .true.
    call getdbla(lines,nl,'initial_tau=',tau,iflag)
    if (iflag /= 0) then
       call getdbla(lines,nl,'tau=',tau,tauFlag)
       if (tauFlag == 0) then
          mAdaptTau = .false.
        else
          tau = propagator_timeStep()
       endif
    end if
    call getdbla(lines,nl,'accept_ratio=',mTargetAR,iflag)

    ds = 1.d0
    call getdbla(lines,nl,'drift_scal=',ds,iflag)
    found = finda(lines,nl,'no_accept_step')
    if (found) then
       ar = .false.
    else
       found = finda(lines,nl,'accept_step')
       if (found) ar = .true.
    endif
    found = finda(lines,nl,'allow_cross')
    if (found) rc = .false.
    found = finda(lines,nl,'reject_cross')
    if (found) rc = .true.

    ! Rothstein/Vrbik E_local/Drift cutoff
    if (finda(lines,nl,'no_elocal_cutoff')) eLocalCutOff=.false.
    cf = 1.d0
    call getdbla(lines,nl,'elocal_cutoff=',cf,iflag)
    call setElocCutOff(eLocalCutOff,tau,cf)

    ! turn on/off autocorrelation calculation
    found = finda(lines,nl,'no_auto_corr')
    if (found) then
       mAutoCorrelation = .false.
    else
       found = finda(lines,nl,'auto_corr')
       if (found) mAutoCorrelation = .true.
    endif

    call getinta(lines,nl,'auto_corr_max=',mAutoCorrMax,iflag)
    if (iflag == 0) mAutoCorrelation = .true.

    ! Initialization of Future Walking

    if(mFuture) then
     call getinta(lines,nl,'tau2=',mTau2,iflag)
     if(iflag /=0) call abortp("qmc_future: tau2 required")
     call getinta(lines,nl,'notau2=',mNotau2,iflag)
     if(iflag /=0) call abortp("qmc_future: notau2 required")
     call getinta(lines,nl,'stprops=',mStprops,iflag)
     if(iflag /=0) call abortp("qmc_future: stprops required")
    endif

    ! Popcontrol and Init Parallel. Old Scheme

    call getinta(lines,nl,'popctrl=',mPopctrl,iflag)
    if (iflag /=0) mPopctrl = POP_GLOBAL

    if (finda(lines,nl,'old_par')) then
      mLoadBalance = .false.
      mPopctrl = POP_LOCAL
    endif

    ! Init for Stochastic Reconfiguration

    if (finda(lines,nl,'rcf')) then
      call getinta(lines,nl,'rcf=',mRcf,iflag)
      mReconf = .true.
      mBranch = .false.
      !!!mKillPersist = .false.
      !!!mPopctrl = 0
      mLoadBalance = .false.
    endif

    sOldLogmode = logmode
    call getinta(lines,nl,'verbose=',vb,iflag)
    if (iflag==0) then
       logmode = vb
    endif

    found = finda(lines,nl,'show_details')
    if (found) mShowDetails = 1

    call propagator_reset()

    if (move == 2) then
        found = finda(lines,nl,'no_exp')
        if (found) call setNoExp()
    endif

    blockdiscard = mStepsDiscard / mBlockLen
    call propagator_init(mWeight,move,tmove,blockdiscard,bgte,tau,ds,ar,rc,mt,mWalkerBlock)

    !Max Analysis - Check for consistent sample size on each core
    !happens if outliers are removed and not replaced on any core
    if (mMaxAnalysis) then
      call myMPIGatherInteger(getSampleSize(sample),1,int_rcv,ierr)
      if (MASTER) then
        do i=1,nproc
          if (int_rcv(i) /= getSampleSize(sample) ) then
            call abortp("(qmc_init): Different sample size on each core for maxima calculation.")
          endif
        enddo
      endif
      call myMPIBarrier(ierr)
    endif

  end subroutine qmc_init


  ! access functions
  real(r8) function qmc_Energy()
     qmc_Energy = sEMean
  end function qmc_Energy

  real(r8) function qmc_StdDev()
     qmc_StdDev = sEMeanStdDev
  end function qmc_StdDev

  real(r8) function qmc_Variance()
     qmc_Variance = sVar
  end function qmc_Variance


  subroutine qmc_writeParams(iu, psimax_obj, rhoMax)
  !----------------------------------------!

    integer, intent(in)                   :: iu
    type(psimax), intent(inout), optional :: psimax_obj
    type(rhoMax_t), intent(inout), optional :: rhoMax
    character(len=10)   :: weight
    character(len=20)   :: s
    real(r8) cf,tau
    logical yn

    select case(mWeight)
    case(WEIGHT_NONE); weight="none"
    case(WEIGHT_REY); weight="Reynolds"
    case(WEIGHT_UMR); weight="Umrigar"
    case(WEIGHT_ACCEPT); weight="Accept"
    case(WEIGHT_TMSM); weight="TMSM"
    case(WEIGHT_TMSC); weight="TMSC"
    case default; call abortp("qmc_writeParams: illegal weight")
    end select

    s = mMethod
    write(iu,'(/3A/)') '   * * *  ',trim(s),' calculation  * * *'
    write(iu,'(A/)') '    QMC parameters:'

    write(iu,'(1X,A21,F12.5,3X,A21,L12)') 'tau =',propagator_timeStep(),' adapt tau =',mAdaptTau
    write(iu,'(1X,2(A21,I12,3X))') ' total walker =',mTargetSampleSizeAllNodes,    &
         ' local walker =',mTargetSampleSize
    write(iu,'(1X,A21,I12,3X,A21,I12)') ' steps =',mSteps,' discard =',mStepsDiscard
    write(iu,'(1X,2(A21,I12,3X))') ' block_len =',mBlockLen,' walker_block =',mWalkerBlock
    write(iu,'(1X,2(A21,I12,3X))') ' step_stride =',mStepStride
    if (mCheckStdDev) write(iu,'(1X,A21,F12.5)') ' target std dev =',mStdDev
    if (mAdaptTau) write(iu,'(1X,A21,F12.5)') 'target accept ratio =',mTargetAR
    write(iu,'(1X,2(A21,F12.5,3X))') 'E_ref =',E_trial,'wfac =',mWFac
    call getElocCutOff(yn,tau,cf)
    write(iu,'(1X,A21,L12,3X,A21,F12.5,3X)') 'E_loc_cutoff =',yn,'factor =',cf
    write(iu,'(1X,A21,L12,3X,A21,I12)') ' kill_persist =',mKillPersist,' max_persist =',mPerMax
    write(iu,'(2(1X,A21,L12,2X))') ' load balance =',mLoadBalance,' branch =',mBranch
    select case (mPopctrl)
    case (POP_NONE)
       write(iu,'(1X,A21,L12,3X,A21,I12)') ' future walking =',mFuture,' pop ctrl = none'
    case (POP_GLOBAL)
       write(iu,'(1X,A21,L12,3X,A21,I12)') ' future walking =',mFuture,' pop ctrl = global'
    case (POP_LOCAL)
       write(iu,'(1X,A21,L12,3X,A21,I12)') ' future walking =',mFuture,' pop ctrl = local'
    case default
       call abortp("qmc_writeParams: illegal mPopctrl value")
    end select
    write(iu,'(1X,A21,L12,3X,A21,I12)') ' Reconf =',mReconf,' RcfMethod =',mRcf
    write(iu,'(1X,A21,L12,3X,A21,I12)') 'accumulate =',mAccumulate
    write(iu,*)

    if (mFuture) call propOutput(iu)

    call propagator_writeParams(iu)

  end subroutine qmc_writeParams



  subroutine qmc_run(sample, psimax_obj, rhoMax)
  !------------------------------------!

    type(RWSample), intent(inout)         :: sample
    type(psimax), intent(inout), optional :: psimax_obj
    type(rhoMax_t), intent(inout), optional :: rhoMax

    integer                     :: block,bs,i,j,k,idx
    integer(i8)                   :: n,st,tauFoundStep
    real(r8)                      :: elocCounter,ERef,varAllNodes
    real(r8)                      :: EMeanAllNodes,sampleSizeAllNodes,eTotal,vTotal
    real(r8)                      :: E,var,stddev,stdDevAllNodes
    real(r8)                      :: EMeanlocal
    real(r8)                      :: VenMeanAllNodes, VeeMeanAllNodes, TMeanAllNodes
    real(r8)                      :: sampleWeightAllNodes
    real(r8)                      :: EL,EECPL,EECPNL
    real(r8)                      :: tcorr,ACvar,ACNcorr
    real(r8)                      :: accepted ! 0.0 for (all) rejected, 1.0 for (all) accepted
#ifdef WTIMER
    real(r8)                      :: wtimer1,wtimer2,ewtime,twtime
#endif
    real(r8), allocatable         :: Estep(:),AC(:,:),ACRingBuffer(:,:),ACresult(:)  ! autocorrelation date
    logical                     :: first,ACavail,tauFound
    type(RandomWalker), pointer :: rwp
    type(RandomWalker), pointer :: rwbp(:)   ! pointer to walker block==array
    type(weightStat) :: ElocStepStat          ! < E_local > (weighted)
    type(weightStat) :: VenStepStat          ! < E_local > (weighted)
    type(weightStat) :: VeeStepStat          ! < E_local > (weighted)
    type(weightStat) :: TStepStat          ! < E_local > (weighted)
    type(weightStat) :: testStat
!!!    type(simpleStat) :: totalAccStat      ! < acceptance ratio >
    type(simpleStat) :: blockAccStat      ! < acceptance ratio > of block
    type(simpleStat) :: adaptAccStat      ! < acceptance ratio > of some steps for adapting tau
    type(simpleStat) :: eRefStat          ! < E_ref >
    type(simpleStat) :: wgtStat           ! < total-Weight >
    type(simpleStat) :: blockStatlocal    ! for local (Procs) BlockStat
    type(stat)       :: autocorrEStat
    type(vectorstat) :: autocorrStat
    type(RWSample) :: history ! used to accumulate samples
    integer :: histSize
    integer(i8) :: dataAllNodes,tStep

    type(DecorrelationData_t), target :: decorrDtot, decorrD
    type(Decorrelation_t) :: decorr

#ifdef WTIMER
    if (MASTER) then
      wtimer1 = omp_get_wtime()
      wtimer = 0.d0
    endif
#endif

    ! init flyvbjerg-petersen
    call decorr%Create(decorrD)

    if (mAccumulate) call internal_accumulateInit()

    if (mWalkerStatistics) call rwStatInit(mStatType)

    if (mReconf .and. mRcf>2) call qmc_reconfInit(sample,mRcf)

    if (mFuture) call future_init(sample,mStprops,mTau2,mNotau2)

    if (mTargetSampleSize==0) mTargetSampleSize = getSampleSize(sample)
    if (mTargetSampleSizeAllNodes==0) mTargetSampleSizeAllNodes = getSampleSizeAllNodes(sample)
    if (mMethod == 'VMC' .and. mTargetSampleSize /= getSampleSize(sample))  &
         call abortp('qmc_run: inconsistent walker numbers in VMC')

    if (mWeight==WEIGHT_NONE) then
       call resetWeights(sample)
       call resetPersistencies(sample)
    endif

    if (mAutoCorrelation) then
       n = getSampleSize(sample)
       allocate(Estep(n),AC(0:mAutoCorrMax,n),ACRingBuffer(mAutoCorrMax,n),ACresult(0:mAutoCorrMax))
       call autocorrStat%create(mAutoCorrMax+1)
       call autocorrEStat%create()
    end if

    call reset(ElocBlockStat)
    call reset(VenBlockStat)
    call reset(VeeBlockStat)
    call reset(TBlockStat)
    call reset(varBlockStat)
    call reset(sTotalStat)
    call reset(blockStatlocal)
    call reset(eRefStat)
    call reset(wgtStat)
    call reset(adaptAccStat)
    call reset(testStat)
    mElocCut = 0; mElocCutCount = 0
    mDriftCut = 0; mDriftCutCount = 0

    sERef = E_trial
    sMode = 0           !
    tStep = 0           ! counter for total steps (no blocking) after discard (for autocorrelation)
    block = 0


    if (MASTER .and. logmode >= 2) then
       if (present(psimax_obj) .and. present(rhoMax)) then
          call qmc_writeParams(iul, psimax_obj, rhoMax)
       else if (present(psimax_obj)) then
          call qmc_writeParams(iul, psimax_obj)
       else if (present(rhoMax)) then
          call qmc_writeParams(iul, rhoMax=rhoMax)
       end if

       if (mMaxAnalysis) then
          write(iul,'(/A)') '             step       size              <E>                   <V>     AR   #max  fmin'
          write(iul,'(A)')  ' --------------------------------------------------------------------------------------'
       else
          write(iul,'(/A)') '             step       size              <E>                   <V>     AR '
          write(iul,'(A)')  ' --------------------------------------------------------------------------'
       end if
    endif

    tauFoundStep = 0
    elocCounter = 0
    ERef = sERef
    call reset(ElocStepStat)
    call reset(VenStepStat)
    call reset(VeeStepStat)
    call reset(TStepStat)
    call reset(blockAccStat)

    STEPS: do st=1,mSteps

       ! loop over walker (in blocks of bs)
       idx = 0
       bs = mWalkerBlock  ! # of walker for simultaneous propagation
       call getFirstWalkerBlock(sample,bs,rwbp)
       first = .true.
       WALKERLOOP: do

          call propagateAndWeight(rwbp,ERef,sMode,accepted)
          elocCounter = elocCounter + bs

          do i=1,bs
             if (st > mStepsDiscard) call decorr%Add_data(E_local(rwbp(i)))

             call addData(ElocStepStat,E_local(rwbp(i)),wgt(rwbp(i)))
             call addData(VenStepStat,E_pot(rwbp(i)) - E_ee(rwbp(i)) - vpot0,wgt(rwbp(i)))
             call addData(VeeStepStat,E_ee(rwbp(i)),wgt(rwbp(i)))
             call addData(TStepStat,E_local(rwbp(i)) - E_pot(rwbp(i)),wgt(rwbp(i)))
!              if (logmode >= 4) then
!                 if (wgt(rwbp(i)) > 3.d0) then
!                    call eloc_eloc_getCurrentElocData1(tphi,tu,teloc,tvpot,txdrift,tydrift,tzdrift)
!                    call eloc_eloc_getCurrentElocEpart1(tkin,tvee,tvne)
!                    call eloc_getECPContribs(tEECPL,tEECPNL)
!                    write(iull,'(2i7,g10.2,3f12.6)') st, idx+i, wgt(rwbp(i)), EL, EECPL, EECPNL
!                 end if
!              end if
             if (mAutoCorrelation) Estep(idx+i) = E_local(rwbp(i))     ! save E_L for autocorrelation
          enddo
          idx = idx + bs
          call addData(blockAccStat,accepted)
          call addData(adaptAccStat,accepted)

          if (mod(st,mStepStride)==0 .and. st > mStepsDiscard) then
             call internal_doStuff()
          endif

       if (.not. isNext(sample)) exit

          bs = mWalkerBlock
          call getNextWalkerBlock(sample,bs,rwbp)
          first = .false.

       enddo WALKERLOOP

       if (mWeight /= WEIGHT_NONE .and. mPopctrl /= POP_NONE) then
          call pop_control(sample,sampleWeightAllNodes,ERef,st)
          if (st > mStepsDiscard) then
             call addData(eRefStat,ERef)
             call addData(wgtStat,sampleWeightAllNodes)
          endif
       endif

       if (mKillPersist .and. .not. mReconf) then
          call kill_persist(sample)
          if (mMethod == 'VMC') call changeSampleSizeTo(sample,mTargetSampleSize)
       endif

       if (mShowSteps .or. logmode>=5) then
          call getSampleEnergyAndVarianceAllNodes(sample,E,var,stddev)
          if (MASTER) write(iul,'(a,2i6,g20.10,3g14.4,i8)') 'STEP:',block+1,st,E,var,stddev
       endif

       if (mBranch) call qmc_branch(sample)

       if (mReconf) call qmc_rcf(sample)

       if (mLoadBalance) call loadBalanceSamples(sample)

       if (mAdaptTau) then
          call qmc_adaptTau(adaptAccStat,tauFound)
          if (tauFound) tauFoundStep = st
       endif

       if (st > mStepsDiscard) then
          if (mFuture) call futurewlk(sample)
          tStep = tStep + 1
          if (mAutoCorrelation) call autocorrAdd()
       endif

       if (mod(st,mBlockLen)==0) then ! end of block

          EMeanlocal = mean(ElocStepStat)
          VenMeanAllNodes = meanAllNodes(VenStepStat)
          VeeMeanAllNodes = meanAllNodes(VeeStepStat)
          TMeanAllNodes = meanAllNodes(TStepStat)
          EMeanAllNodes = meanAllNodes(ElocStepStat)
          stdDevAllNodes = stdDevMeanAllNodes(ElocStepStat)
          varAllNodes = varianceAllNodes(ElocStepStat)

          if (st <= mStepsDiscard) then
             if(mPopctrl == POP_GLOBAL) then
                sERef = EMeanAllNodes
             else if (mPopctrl == POP_LOCAL) then
                sERef = EMeanlocal
             endif
          else
             sTotalStat = sTotalStat + ElocStepStat     ! local statistics
             call addData(ElocBlockStat,EMeanAllNodes)
             call addData(VenBlockStat,VenMeanAllNodes)
             call addData(VeeBlockStat,VeeMeanAllNodes)
             call addData(TBlockStat,TMeanAllNodes)
             call addData(varBlockStat,varAllNodes)
             call addData(blockStatlocal,Emeanlocal)
             if (mPopctrl == POP_GLOBAL) then
                sERef = mean(ElocBlockStat)
             else if (mPopctrl == POP_LOCAL) then
                sERef = mean(blockStatlocal)
             endif
          endif

          if (mWeight == WEIGHT_UMR) call propagator_setERef0(sERef)  ! for Umrigar weighting

          sampleSizeAllNodes = getSampleSizeAllNodes(sample)
          if (logmode>=2) then
             if (present(psimax_obj) .and. mMaxAnalysis .and. st>mStepsDiscard) then
                write(iul,'(i15,i10,f16.5,a4,f10.5,f10.3,f8.3,i4,f13.5)') st, nint(sampleSizeAllNodes), &
                     EMeanAllNodes, ' +/-', stdDevAllNodes, varAllNodes, mean(blockAccStat), &
                     psimax_obj%get_diff_max(), psimax_obj%get_first_f()
             else
                write(iul,'(i15,i10,f16.5,a4,f10.5,f10.3,f8.3)') st, nint(sampleSizeAllNodes), &
                     EMeanAllNodes, ' +/-', stdDevAllNodes, varAllNodes, mean(blockAccStat)
             end if
          end if

          block = INT(st/mBlockLen)
          call propagator_endOfBlock(block)  ! allow propagator to do stuff at end of block

          if (mWalkerStatistics) then
             call rwStatPrint(42,st)
             call rwStatReset()
          endif


          if (mCheckStdDev .and. st > mStepsDiscard+10*mBlockLen) then
             if (stdDevMean(ElocBlockStat) < mStdDev) exit STEPS
          endif

          elocCounter = 0
          ERef = sERef
          call reset(ElocStepStat)
          call reset(VenStepStat)
          call reset(VeeStepStat)
          call reset(TStepStat)
          call reset(blockAccStat)

       endif ! end of block

    enddo STEPS

    eTotal = meanAllNodes(sTotalStat)
    vTotal = varianceAllNodes(sTotalStat)
    dataAllNodes = dataCountAllNodes(sTotalStat)

    if (st > mStepsDiscard + mBlockLen) then

      if (mAutoCorrelation) then
         ACavail = autocorrStat%count() > 0
         if (ACavail) ACresult = autocorrStat%mean() - (autocorrEStat%mean())**2
      end if


#ifdef MPI
      call myMPIReduceDecorr(decorrD, decorrDtot)
#else
      decorrDtot = decorrD
#endif
      if (Master) then
        call decorr%Create(decorrDtot)

        if (logmode >= 3) then
            write(988,*) "Decorrelation object:"
            call decorr%Write(988)
        end if

        if (mAutoCorrelation .and. ACavail) then
           ACvar = ACresult(0)
           ACresult = ACresult / ACresult(0)
           do i=1,mAutoCorrMax
              if (ACresult(i) < 0.05) exit
           end do
           if (i==1) then
              ACNcorr = 1
           else if (i==mAutoCorrMax+1) then
              ACNcorr = 0     ! denotes larger than mAutoCorrMax
           else
              ! linear interpolation between R(i-1) and R(i)
              ACNcorr = i-1 + (0.05-ACresult(i-1))*1.d0/(ACresult(i)-ACresult(i-1))
           end if
           !! write preliminarily full autocorrelation data
           if (logmode >= 3) then
              write(987,*) ' full autocorrelation data:'
              write(987,*) i
              do i=0,size(ACresult)-1
                 write(987,'(i5,f12.5)') i,ACresult(i)
              end do
           end if
        end if

        ! calculation of correlation length
        ! note that the block statistic contains locally the block means of ALL nodes (see above)
        tcorr = variance(ElocBlockStat)/vTotal * &
                dble(dataAllNodes)/dble(datacount(ElocBlockStat))

        ! write final output
        write(iul,'(//2A/)') '  FINAL RESULT:'
        write(iul,'(a,f13.5,a,f8.5,a)') ' total energy (mean E_loc)    = ',eTotal,' +/-',  &
                stdDevMean(ElocBlockStat),' E_h'
        write(iul,'(a,f13.5,a,f8.5,a)') ' kinetic energy               = ',mean(TBlockStat),' +/-',  &
                stdDevMean(TBlockStat),' E_h'
        write(iul,'(a,f13.5,a,f8.5,a)') ' e-n potential energy         = ',mean(VenBlockStat),' +/-',  &
                stdDevMean(VenBlockStat),' E_h'
        write(iul,'(a,f13.5,a,f8.5,a)') ' e-e potential energy         = ',mean(VeeBlockStat),' +/-',  &
                stdDevMean(VeeBlockStat),' E_h'
        write(iul,'(a,f13.5,a)')        ' n-n potential energy         = ',vpot0,' E_h'
        write(iul,'(a,f13.5,a,f8.5,a)') ' variance (of E_loc)          = ',vTotal, ' +/-', &
                stdDevMean(varBlockStat), ' E_h^2'
        write(iul,'(a,f13.5,A)') ' block average variance       = ',mean(varBlockStat),' E_h^2'
        if (mWeight/=WEIGHT_NONE .and. mPopctrl /= POP_NONE) then
           write(iul,'(a,f13.5,a,f8.5,a)')       ' mean E_ref (sigma_i)         = ',mean(eRefStat),' +/-', &
                sqrt(variance(eRefStat)),' E_h'
           write(iul,'(a,f13.2,a,f8.2)')         ' mean weight (sigma_i)        = ',mean(wgtStat),' +/-', &
                sqrt(variance(wgtStat))
           write(iul,'(a,f13.2,a,f13.2)')        ' minimum weight               = ',minValue(wgtStat), &
                ' maximum weight = ',maxValue(wgtStat)
           write(iul,'(a,f13.4)')                ' tau_acc                      = ',propagator_timeStepAcc()
        endif
        if (tauFoundStep > 0) then
           write(iul,'(a,f13.4,a,i12)')          ' tau (adapted)                = ',propagator_timeStep(), &
                                                 ' fixed at step ',tauFoundStep
        endif
        if (mAutoCorrelation .and. ACavail) then
           if (ACNcorr > 0) then
              write(iul,'(a,f9.1)') ' N_corr (<5%)                 = ',ACNcorr
           else
              write(iul,'(a,i9)') ' N_corr (<5%)                 > ',mAutoCorrMax
           end if
        end if
        write(iul,'(a,f9.1)')  ' N_corr (global)              = ',tcorr
        write(iul,'(a)')  ''
        write(iul,'(a)')  '  FLYVBJERG-PETERSEN: '
        call decorr%Write_results(iul)
        call decorr%Destroy

        if (mShowDetails >= 1) then
           if (propagator_getTMove() > 0) then
              write(iul,'(a,f9.5)') ' T move ratio                 = ', propagator_getTMoveRatio()
           end if
           if (mElocCutCount > 0) then
              write(iul,'(a,f9.5)') ' E_local cutoff ratio         = ', dble(mElocCut) / dble(mElocCutCount)
              write(iul,'(a,f9.5)') ' drift cutoff ratio           = ', dble(mDriftCut) / dble(mDriftCutCount)
           end if
        end if
        if (datacount(ElocBlockStat) < 20) then
           write(iul,'(/a/a)') ' WARNING: stddev and global N_corr may be unreliable', &
             '  (number of blocks not discarded < 20)'
        end if
        if (20*tcorr > mBlockLen) then
           write(iul,'(/a/a)') ' WARNING: stddev and global N_corr may be unreliable', &
             '  (block_len < 20*N_corr)'
        end if

      end if

    else

      if (MASTER) then
        write(iul,'(/A,F15.5,A,F15.5)') ' qmc: Emean = ',EMeanAllNodes,' var = ',varAllNodes
      endif

    endif

    if (mFuture)then
       call propCalculate()
       if (MASTER) call propPrint()
       call deallocFWArrays()
    endif

    ! save result for access with qmc_Energy etc.
    if (st > mStepsDiscard + mBlockLen) then
       sEMean = eTotal
       sEMeanStdDev = stdDevMean(ElocBlockStat)
       sVar = vTotal
       sVarError = stdDevMean(varBlockStat)
    else
       sEMean = EMeanAllNodes
       sVar = varAllNodes
       sEMeanStdDev = ieee_value(sEMeanStdDev, ieee_quiet_nan)
       sVarError = ieee_value(sVarError, ieee_quiet_nan)
    endif
    ! save result also in global module
    call setCurrentResult(sEMean,sEMeanStdDev,sVar,sVarError)

    if (present(psimax_obj) .and. mMaxSearch) then
       call psimax_obj%write_results()
       call psimax_obj%destroy()
    end if
    if (present(rhoMax)) then
       if (rhoMax%isInitialized()) then
          call rhoMax%writeResults(iul)
          call rhoMax%destroy()
       end if
    end if

    logmode = sOldLogmode

    if (mAccumulate) then
      ! remember the actual sample size
      mAccumulateSampleSize = getSampleSize(sample)

      ! replace the current sample with the "history" sample
      call destroySample(sample)
      sample = history

      call destroySample(history)

      n = getSampleSizeAllNodes(sample)
      if (MASTER .and. logmode >= 2) then
        write(iul,'(/a,i12)') " sample accumulation: new total sample size is ", n
      endif
    endif

    if (mAutoCorrelation) then
      deallocate(Estep,AC,ACRingBuffer,ACresult)
      call autocorrStat%destroy()
      call autocorrEStat%destroy()
    endif

    if (mReconf .and. mRcf>2) then
      call qmc_reconfDestroy()
    endif

#ifdef WTIMER
    if (MASTER .and. logmode >= 2) then
      wtimer2 = omp_get_wtime()
      write(iul,'(/a)') " qmc wall clock times (master, in sec and relative to eloc):"
      twtime = wtimer2-wtimer1
      write(iul,'(a,f20.3)')       " total:    ",twtime
      write(iul,'(a,f20.3,f10.3)') " eloc:     ",wtimer(WTIMER_ELOC),wtimer(WTIMER_ELOC)/twtime
      ewtime = wtimer(WTIMER_ELOC)
      write(iul,'(a,f20.3,f10.3)') " phi:      ",wtimer(WTIMER_PHI),wtimer(WTIMER_PHI)/ewtime
      write(iul,'(a,f20.3,f10.3)') " jastrow:  ",wtimer(WTIMER_JAS),wtimer(WTIMER_JAS)/ewtime
      write(iul,'(a,f20.3,f10.3)') " pseudo:   ",wtimer(WTIMER_PP),wtimer(WTIMER_PP)/ewtime
      write(iul,'(a,f20.3,f10.3)') " AO:       ",wtimer(WTIMER_AO),wtimer(WTIMER_AO)/ewtime
      write(iul,'(a,f20.3,f10.3)') " MO:       ",wtimer(WTIMER_MO),wtimer(WTIMER_MO)/ewtime
      write(iul,'(a,f20.3,f10.3)') " AOMO:     ",wtimer(WTIMER_AOMO),wtimer(WTIMER_AOMO)/ewtime
      write(iul,'(a,f20.3,f10.3)') " mdet:     ",wtimer(WTIMER_MDET),wtimer(WTIMER_MDET)/ewtime
    endif
#endif

  contains

    subroutine internal_accumulateInit()
      ! if this is an accumulation run, the sample that was passed in is actual
      ! the history of previously seen samples. if this is the first accumulation
      ! run after creating a sample, store this size as the number of actual
      ! walkers used

      history = sample

      if (getSampleSize(sample) == currentSampleSize()) then
        mAccumulateSampleSize = getSampleSize(sample)
      endif

      if (MASTER .and. logmode >= 3) then
        write(iul,'(a,i8)') " accumulate run: sample copied to history. Initial (local) size: ", getSampleSize(history)
      endif

      ! the actual walkers used in this run are now only the last
      ! |mAccumulateSampleSize| samples in the history
      if (mAccumulateSampleSize < getSampleSize(sample)) then
         call reduceToLast(sample, mAccumulateSampleSize)
         if (MASTER .and. logmode >= 3) then
            write(iul,'(a,i8)') " accumulate: reduced actual sample to ", mAccumulateSampleSize
         endif
      endif

      ! pre-allocate the new history to avoid expensive reallocation in the
      ! qmc loop
      histSize = getSampleSize(history) + INT(mSteps / mStepStride) * mAccumulateSampleSize
      call reallocateSample(history, histSize)
    end subroutine internal_accumulateInit

    subroutine internal_doStuff()
      integer :: w, i

       if (present(psimax_obj) .and. mMaxSearch) then
          do w = 1, size(rwbp)
             if (mMethod == 'DMC') then
                if (psimax_obj%get_analyze_mode() /= RAWDATA &
                .and. psimax_obj%get_analyze_mode() /= NOANALYSIS) then
                   call abortp("internal_doStuff: on the fly dmc sample maximization"&
                           //" is atm only possible with rawdata mode (or no analysis)")
                else
                   do i = 1, INT(wgt(rwbp(w)) + myran())
                      ! this is of course unneccessarily costly
                      call psimax_obj%opt(rwbp(w), update_walker=.false.)
                   end do
                end if
             else if (mMethod == 'VMC') then
                call psimax_obj%opt(rwbp(w), update_walker=.false.)
             end if
          end do
       end if

      if (present(rhoMax)) then
         if (rhoMax%isInitialized()) then
            do w = 1, size(rwbp)
               if (mMethod == 'DMC') then
                  do i = 1, INT(wgt(rwbp(w)) + myran())
                     ! this is of course unneccessarily costly
                     call rhoMax%opt(rwbp(w), update_walker=.false.)
                  end do
               else if (mMethod == 'VMC') then
                  call rhoMax%opt(rwbp(w), update_walker=.false.)
               end if
            end do
         end if
      end if

       if (mWalkerStatistics) then
          do i = 1, bs
             rwp => rwbp(i)
             call rwAddToStat(rwp)
          enddo
       endif

       if (mAccumulate) then
          do i = 1, bs
            call appendWalker(history, rwbp(i))
          enddo
       endif

       ! collecting trajectory data
       if (mTrajectory .and. first) call display(rwbp(1),2,43)

    end subroutine internal_doStuff

    subroutine autocorrAdd()
       if (idx /= getSampleSize(sample)) call abortp("autocorrAdd: incorrect idx value")
       if (tStep > mAutoCorrMax) then
         do i=1,idx
            AC(0,i) = Estep(i)*Estep(i)
            do k=1,mAutoCorrMax
               j = INT(mod(tStep-k-1,INT(mAutoCorrMax, i8)) + 1)
               AC(k,i) = ACRingBuffer(j,i)*Estep(i)
            end do
            j = INT(mod(tStep-1,INT(mAutoCorrMax, i8)) + 1)
            ACRingBuffer(j,i) = Estep(i)
            call autocorrStat%add(AC(:,i))
            call autocorrEStat%add(Estep(i))
         end do
       else ! fill ring buffer
         j = INT(mod(tStep-1,INT(mAutoCorrMax, i8)) + 1)
         do i=1,idx
            ACRingBuffer(j,i) = Estep(i)
         end do
       end if

    end subroutine autocorrAdd

  end subroutine qmc_run



  subroutine qmc_branch(sample)
  !---------------------------!

  ! different branching schemes for current walker sample

    real(r8), parameter :: KILL_THRESH=1.d-12
    real(r8), parameter :: SAMPLE_SIZE_FACTOR=0.9d0

    type(RWSample), intent(inout) :: sample

    integer  ps,sampleSize,nsplit,i,maxSize
    real(r8)   wNew,xi,w,wJoin
    type(RandomWalker), pointer :: rwp,rwpJoin


    ! Kill walkers with weight < KILL_THRESH
    rwp => getFirst(sample)
    do
       if (wgt(rwp) < KILL_THRESH) then
          call deleteCurrentWalker(sample)
          if (.not.isValid(sample)) exit     ! true after deleting last element
       else
          if (.not.isNext(sample)) exit
          rwp => getNext(sample)
       endif
    enddo

    if (mJoin) then  ! join two walkers with weight < mLThresh (Umrigar, 1993)

       rwpJoin => null()                     ! index of 1st walker
       rwp => getFirst(sample)                ! oder "=>" ???
       do
          w = wgt(rwp)
          if (w < mLThresh) then
             if (.not.associated(rwpJoin)) then          ! 1st walker for joining
                wJoin   = w
                rwpJoin  => rwp
                if (.not.isNext(sample)) exit
                rwp => getNext(sample)
             else                            ! 2nd found
                xi = myran()
                wNew = w + wJoin
                if (xi < w/wNew) then
                   call setWeight(rwp,wNew)
                   rwpJoin = rwp
                else
                   call setWeight(rwpJoin,wNew)
                endif
                call deleteCurrentWalker(sample)
                if (.not.isValid(sample)) exit
                rwpJoin => null()
             endif
          else
             if (.not.isNext(sample)) exit
             rwp => getNext(sample)
          endif
       enddo

    else if (mKill) then

       rwp => getFirst(sample)
       do
          w = wgt(rwp)
          if (w < mLThresh) then
             xi = myran()
             wNew = 2.d0*mLThresh
             if (xi < w/wNew) then
                ! walker survives
                call setWeight(rwp,wNew)
                if (.not.isNext(sample)) exit
                rwp => getNext(sample)
             else
                call deleteCurrentWalker(sample)
                if (.not.isValid(sample)) exit
             endif
          else
             if (.not.isNext(sample)) exit
             rwp => getNext(sample)
          endif
       enddo

    endif

    if (mSplit) then

       sampleSize = getSampleSize(sample)
       maxSize    = INT(SAMPLE_SIZE_FACTOR*getMaxSampleSize(sample))
       rwp => getFirst(sample)
       ps = 0
       do
          w = wgt(rwp)
          if (w > mUThresh) then
             nsplit = int(w)
             ! do not split if maxSize reached
             if (logmode >= 3 .and. nsplit>3) write(iul,'(A,I6)') 'BRANCH: nsplit=',nsplit
             if (getSampleSize(sample)+nsplit-1 > maxSize) then
                write(999,*) "WARNING:qmc_branch: maxSize hit",getMyTaskId(),getSampleSize(sample)+nsplit-1,maxSize
                exit
             end if
             wNew = w/nsplit
             call setWeight(rwp,wNew)
             do i=1,nsplit-1
                call appendWalker(sample,rwp)
             enddo
          endif
          ps = ps +1
          if (ps == sampleSize) exit
          rwp => getNext(sample)
       enddo

    endif


  end subroutine qmc_branch



  subroutine qmc_rcf(sample)
  !-------------------------!
  ! reconfiguration of sample

  type(RWSample), intent(inout) :: sample
  type(RandomWalker), pointer :: rwp

  select case (mRcf)
    case (1)
      call qmc_reconf1(sample)
    case (2)
      call qmc_reconf2(sample)
    case (3:4)
      call qmc_reconfNew(sample)
    case (5:6)
      ! kill persistent walkers by setting weight to zero
      rwp => getFirst(sample)
      do
        if (persist(rwp) > mPerMax) call setWeight(rwp,0.d0)
        if (.not.isNext(sample)) exit
        rwp => getNext(sample)
      enddo
      call qmc_reconfNew(sample)
    case default
      call abortp("qmc_rcf: illegal value for mRcf")
  end select

  end subroutine qmc_rcf


  subroutine pop_control(sample,sampleWeightAllNodes,ERef,st)
  !---------------------------------------------------------!

  ! modify reference energy for population control
  ! Popcontrl global = 1 ; local = 2

  type(RWSample), intent(inout) :: sample
  real(r8),intent(out)            :: ERef
  integer(i8),intent(in)          :: st
  real(r8)                        :: Samplesize
  real(r8),intent(out)            :: sampleWeightAllNodes
  real(r8)                        :: sampleWeight
  integer                       :: ierr


    if(mPopctrl == POP_GLOBAL) then
       Sampleweight = getSampleWeightAllNodes(sample)
       Samplesize   = mTargetSampleSizeAllNodes
       sampleWeightAllNodes = getSampleWeightAllNodes(sample)
    elseif(mPopctrl == POP_LOCAL) then
       Sampleweight = getSampleWeight(sample)
       Samplesize   = mTargetSampleSize
       sampleWeightAllNodes = getSampleWeightAllNodes(sample)
    endif
    ERef = sERef - mWFac*log(Sampleweight/Samplesize) !
    if (logmode >= 3) then
       write(iul,'(A,2F15.4,I6,2F15.4,I15)') 'POP:',Sampleweight,getSampleWeight(sample),   &
             getSampleSize(sample),ERef,sERef,st
       flush(iul)
    endif
    if (sampleWeight > 2*Samplesize) then
       write(iul + mytid,'(A,G13.3)') ' sample overflow: all weights=',Sampleweight
       write(iul + mytid,'(A,2G18.6)') ' sample overflow: ERef,sERef=',ERef,sERef
       call displaySample(sample,iul + mytid)
       call myMPIBarrier(ierr)
       call abortp('QMC: walker overflow')
    endif

  end subroutine pop_control


  subroutine kill_persist(sample)
  !-----------------------------!
    ! kill all persistent walkers from sample
    type(RWSample), intent(inout) :: sample

    type(RandomWalker), pointer :: rwp

    rwp => getFirst(sample)
    do
       if (persist(rwp) > mPerMax) then
          call deleteCurrentWalker(sample)
          if (.not.isValid(sample)) exit
       else
          if (.not.isNext(sample)) exit
          rwp => getNext(sample)
       endif
    enddo

  end subroutine kill_persist



  subroutine qmc_adaptTau(accStat,tauFound)
  !---------------------------------------!

    ! adapting time step tau to the target acceptance ratio (AR)
    ! as accumulated in accStat
    ! note: adaptation is based on the model AR = exp(-a*tau)
    ! or tau = - log(AR)/a
    integer, parameter             :: minAdaptSteps = 500
    type(simpleStat),intent(inout) :: accStat      ! < acceptance ratio > of some steps for adapting tau
    logical, intent(out)           :: tauFound
    real(r8) accRatio,oldTau,newTau

    tauFound = .false.

    if (dataCountAllNodes(accStat) < minAdaptSteps) return

    accRatio = meanAllNodes(accStat)

    ! adapt tau even if convergence reached
    oldTau = propagator_timeStep()
    newTau = oldTau * log(mTargetAR) / log(accRatio)

    ! set tau to mMaxTau, if larger. Then stop adaptation
    if (newTau > mMaxTau .or. .not. ieee_is_normal(newTau)) then
       newTau = mMaxTau
    end if

    call propagator_setTimeStep(newTau)

    if ((abs(accRatio-mTargetAR) < mEpsAR) .or. (newTau == mMaxTau)) then
       mAdaptTau = .false.
       tauFound = .true.
    else
       call reset(accStat)
    endif

  end subroutine qmc_adaptTau



  subroutine initWalkerStat(lines,nl)
  !---------------------------------!

  ! read and initialize $walker_stat section

  integer                     :: nl
  character(len=120)          :: lines(nl)
  character(len=1)            :: statType

  if (nproc > 1) call abortp('walkerstat only in serial runs')
  statType = 'c'
  if (finda(lines,nl,'spherical')) statType = 's'

  mWalkerStatistics = .true.
  mStatType = statType

  open(42,file=trim(baseName)//'.reg')

  end subroutine initWalkerStat



  subroutine initTrajectory(lines,nl)
  !---------------------------------!

  ! read and initialize $gen section

  use global_m
  integer                     :: nl
  character(len=120)          :: lines(nl)

  if (nproc > 1) call abortp('trajectory only in serial runs')
  open(43,file=trim(baseName)//'.trj')

  mTrajectory = .true.

  call trajectory_head(43)

  end subroutine initTrajectory

end MODULE qmc_m
