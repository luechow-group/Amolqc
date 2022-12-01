! Copyright (C) 1996, 1999, 2012-2015 Arne Luechow
! Copyright (C) 2013, 2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module qmcSample_m

  use kinds_m, only: r8
  use global_m
  use error_m
  use parsing_m
  use statistics
  use randomWalker_m
  use rwSample_m
  use propagator_m
  use waveFunction_m
  use initialPositions_m
  use wfData_m, only: atoms
  use mpiInterface_m, only: myMPIGatherDouble

  implicit none

  private
  public  :: sample, initInitialWalker, initSampleSize, currentSampleSize, displaceInitialWalker

  integer :: mSampleType=1         ! chooses method to generate initial sample
  integer :: mSampleMode=1         ! 'create','read','change_size'
  integer :: mWalkerGenMode=1      ! simple gaussian (1), density (2), or lmocreate (3)
  integer :: mInitSampleSize=0     ! initial sample size
  integer :: mCurrentSampleSize=0  ! current sample size (set by size= or new_size=)
  ! the following only for VMC generation of sample with one walker (sample1)
  integer :: mBlockSize=100        ! steps in block
  integer :: mBlockDiscard=5       ! blocks to discard initially
  integer :: mMove=1               ! -> mMove in propagator
  integer :: mPrint=0              ! print level
  integer :: mDWGen=0              ! Descenden Weighting Gen
  real(r8)  :: mSampleTau=-1.0d0     ! (initial) tau for sample building (neg: calc tau from nuc charges)
  real(r8)  :: mTargetAccRatio=0.7d0 ! desired acceptance ratio for VMC sample building
  real(r8)  :: mEmin,mEmax           ! min/max values for gaussian samples
  logical :: mRejCross=.false.     ! -> mRejCross in propagator
  character(len=3)  :: mMoveType='all'  ! 'all' | 'one'
  real(r8), allocatable :: mX(:), mY(:), mZ(:)

  type(TMoveControl) :: mTMove ! -> mTMove in propagator

contains

  !------------------------------!
  subroutine sample(lines,nl,smpl)
  !------------------------------!

    integer, intent(in)           :: nl
    character(len=120), intent(in) :: lines(nl)
    type(RWSample), intent(inout) :: smpl
    integer iflag,newSize,nNew,nOld,oldLogmode,v,n,ierr
    integer riter,nbins,gsize,removeMode,maxsize
    real(r8)  :: mpiV1(5),mpiV2(5*MaxMPINodes)
    real(r8)  :: Emean,EmeanAll,var,varAll,sigma,sigmaAll,fac,rfac,tol,varscale
    character(len=15)           :: typ,start
    character(len=3)            :: s
    character(len=80)           :: posfilename
    logical                     :: found,replace,global,keepLast

    call assert(ne > 0,"initSampleParams: wf is not yet read")
    oldLogmode = logmode

    ! read $sample section

    if (finda(lines,nl,'create')) then
       mSampleMode = 1
       if (logmode>=2) write(iul,'(A)') ' creating new sample'
    else if (finda(lines,nl,'read')) then
       mSampleMode = 2
       if (logmode>=2) write(iul,'(A)') ' reading sample'
    else if (finda(lines,nl,'change_size')) then
       mSampleMode = 3
       if (logmode>=2) write(iul,'(A)') ' changing size of sample'
    else if (finda(lines,nl,'histogram')) then
       mSampleMode = 4
       if (logmode>=2) write(iul,'(A)') ' creating histogram of sample'
    else if (finda(lines,nl,'remove_outliers')) then
       mSampleMode = 5
       if (logmode>=2) write(iul,'(A)') ' removing outliers'
    else if (finda(lines,nl,'unweight')) then
       mSampleMode = 6
       if (logmode>=2) write(iul,'(A)') ' unweighting sample'
    else
       mSampleMode = 1
       if (logmode>=2) write(iul,'(A)') ' creating new sample'
    end if
    call getinta(lines,nl,'verbose=',v,iflag)
    if (iflag==0) then
       oldLogmode = logmode
       logmode = v
    endif

    typ = 'random'
    call getstra(lines,nl,'generate=',typ,iflag)
    if (typ=='vmc') then
       mSampleType=1
    else if (typ=='random') then
       mSampleType=2
    else if (typ=='single_point') then
       mSampleType=4
    else if (typ=='pos_file') then
       mSampleType=5
    else
       call abortp('sample: generate mode not implemented')
    endif

    select case (mSampleMode)
    case (2)
       call getstra(lines,nl,'pos_file=',posfilename,iflag)
       call assert(iflag==0,'sample: read requires file argument')

       call assert(.not. (finda(lines,nl,"size=") .and. finda(lines,nl,"size_total=")) &
               ,'sample: give size or size_total, not both.')

       if (finda(lines,nl,"size=")) then
          call getinta(lines,nl,'size=',mInitSampleSize,iflag)
       else if (finda(lines,nl,"size_total=")) then
          call getinta(lines,nl,'size_total=',mInitSampleSize,iflag)
          ! note the integer division in the next line, some walkers may be lost
          mInitSampleSize = mInitSampleSize / nproc
       else
          call error('sample: read requires size or size_total argument')
       end if
    case (4:5)
       removeMode = 2
       call getinta(lines,nl,'algo=',removeMode,iflag)
       tol = 3.0
       call getdbla(lines,nl,'tol=',tol,iflag)
       replace = .true.
       found = finda(lines,nl,'no_replace')
       if (found) replace = .false.
       if (PARALLEL_RUN) then
          global = .true.
       else
          global = .false.
       end if
       found = finda(lines,nl,'global')
       if (found) global = .true.
       found = finda(lines,nl,'local')
       if (found) global = .false.
       rfac = 10.d0
       call getdbla(lines,nl,'remove_factor=',rfac,iflag)
       riter = 100
       call getinta(lines,nl,'iterations=',riter,iflag)
       fac = 5.d0
       call getdbla(lines,nl,'histogram_width=',fac,iflag)
       nbins = 20
       call getinta(lines,nl,'bins=',nbins,iflag)
    end select

    start = 'density'
    call getstra(lines,nl,'start=',start,iflag)
    if (start=='gaussian') then
       mWalkerGenMode = 1
    else if (start=='density') then
       mWalkerGenMode = 2
       call getdbla(lines,nl,'rspread=',fac,iflag)
       if (iflag==0) call wf_setRadiusSpread(fac)
    else if (start=='lmo') then
       mWalkerGenMode = 3
       gsize = 100
       call getinta(lines,nl,'grid_size=',gsize,iflag)
       varscale = 1.0d0
       call getdbla(lines,nl,'scale=',varscale,iflag)
       call init_lmoSampling(gsize,varscale)
    else
       call abortp('sample: start mode not implemented')
    endif

    mEmin = 0.d0
    call getdbla(lines,nl,'E_min=',mEmin,iflag)
    call getdbla(lines,nl,'E_max=',mEmax,iflag)
    if (iflag==0) call assert(mEmin<mEmax,'sample: E_min must be smaller than E_max')
    if (iflag==0 .and. typ=='random') mSampleType = 3
    select case (mSampleMode)
    case (1)
       call getinta(lines,nl,'size=',mInitSampleSize,iflag)
       mCurrentSampleSize = mInitSampleSize
       if (iflag /= 0) call abortp('$sample: create: size required')
    case (3)
       if (finda(lines,nl,'new_size=init_size')) then
          newSize = mInitSampleSize
       else
          call getinta(lines,nl,'new_size=',newSize,iflag)
          if (iflag /= 0) call abortp('$sample: change_size: new_size required')
       end if
       mCurrentSampleSize = newSize

       keepLast = .false.
       if(finda(lines,nl,'last')) keepLast = .true.
    end select

    call getinta(lines,nl,'ndwgen=',mDWGen,iflag)
    if (mDWGen /= 0) call setDescendentWeightingGenerations(mDWGen)

    call getinta(lines,nl,'size_max=',maxsize,iflag)
    if (iflag /= 0) maxsize = INT(1.3_r8 * mInitSampleSize)

    call getinta(lines,nl,'print=',mPrint,iflag)

    if (typ == 'vmc') then
       call getinta(lines,nl,'steps=',mBlockSize,iflag)
       call getinta(lines,nl,'discard=',mBlockDiscard,iflag)
       call getdbla(lines,nl,'tau=',mSampleTau,iflag)
       call getdbla(lines,nl,'acc_ratio=',mTargetAccRatio,iflag)
       call getstra(lines,nl,'move_typ=',mMoveType,iflag)
       call getstra(lines,nl,'move=',s,iflag)
       if (iflag == 0) then
          select case(s)
          case('Rey'); mMove=1
          case('Umr'); mMove=2
          case default; call abortp("initSampleParams: illegal move")
          end select
       endif
       call getstra(lines,nl,'T_moves=',s,iflag)
       if (iflag == 0) then
          select case(s)
          case('simple'); mTMove%mode=1
          case('all_electrons'); mTMove%mode=2
          case default; call abortp("initSampleParams: illegal T_move")
          end select
       endif
       call getloga(lines,nl,'reject_cross',mRejCross,iflag)
    endif


    ! do sample command

    select case (mSampleMode)
    case (1:2)
       if (mSampleMode==1) then
          call getInitialSample(smpl,maxsize)
       else
          call readInitialSample(smpl,maxsize,posfilename)
          if (MASTER .and. logmode >= 2) then
             write(iul,'(/2A/)') ' created initial sample from file ',trim(posfilename)
          endif
       end if
       call getSampleEnergyAndVariance(smpl,Emean,var,sigma)
       if (PARALLEL_RUN) then
          call getSampleEnergyAndVarianceAllNodes(smpl,EmeanAll,varAll,sigmaAll)
          n = getSampleSizeAllNodes(smpl)
       else
          n = getSampleSize(smpl)
       endif

       if (MASTER .and. logmode >= 1) then
          if (PARALLEL_RUN) then
             write(iul,'(A,I12)')          ' sample size (all nodes): ',n
             write(iul,'(3(A,F12.4))') ' <E>_allnodes = ',EmeanAll,' +/- ',sigmaAll,' V_allnodes = ',varAll
          else
             write(iul,'(A,I12)')       ' sample size: ',n
             if (n > 1) then
                write(iul,'(3(A,F12.4))') ' <E>  = ',Emean,' +/- ',sigma,      ' V   = ',var
             else
                write(iul,'(3(A,F12.4))') '  E   = ',Emean
             endif
          endif
          write(iul,*)
       endif

       ! collect sample info from all nodes
       if (PARALLEL_RUN) then
          mpiV1(1)=mytid; mpiV1(2)=getSampleSize(smpl); mpiV1(3)=Emean; mpiV1(4)=sigma; mpiV1(5)=var
          call myMPIGatherDouble(mpiV1,5,mpiV2,ierr)
          if (MASTER .and. logmode >= 3) then
             write(iul,*) ' mytid   size                     <E>              <V>   '
             write(iul,*) '---------------------------------------------------------'
             do n=0,nproc-1
                write(iul,'(I5,I10,3X,F14.4,A,2F14.4)') nint(mpiV2(5*n+1)),nint(mpiV2(5*n+2)),  &
                                                mpiV2(5*n+3),' +/-',mpiV2(5*n+4),mpiV2(5*n+5)
             enddo
             write(iul,*)
          endif
       endif

    case (3)  ! change size
       nOld = getSampleSize(smpl)
       if(keepLast) then
          call reduceToLast(smpl,newSize)
       else
          call changeSampleSizeTo(smpl,newSize)
       endif
       nNew = getSampleSize(smpl)
       call assert(newSize==nNew,"changeSampleSizeTo failed")
       if (logmode >= 2) then
          if (nNew < nOld) then
             write(iul,'(a,i7,a,i7,a)') ' sample size has been reduced by ',nOld-nNew,' to ',nNew,' walkers per node'
             write(iul,*)
          else if (nNew > nOld) then
             write(iul,'(a,i7,a)') ' sample size has been increased by ',nNew-nOld,' (random duplication) '
             write(iul,'(a,i7,a)') ' to ', nNew,' walkers per node'
             write(iul,*)
          else
             write(iul,'(/a/)') ' sample size is not changed'
          endif
       endif
    case (4)
       call assert(getSampleSize(smpl)>40,'sample: sample size too small for histogram')
       if (logmode >= 2) then
       endif
       call getSampleEnergyAndVariance(smpl,Emean,var,sigma)
       call writeSampleHistogram(smpl,Emean,var,fac,nbins)
    case (5)
       if (removeMode==1) then
          call assert(getSampleSize(smpl)>2*nbins,'sample: sample size too small for histogram')
       end if
       if (logmode >= 2) then
          if (removeMode==2) then
             if (global) then
                write(iul,'(/A,I3)') '  using trimmed mean (parallel version)'
             else
                write(iul,'(A,I3)') '  using trimmed mean (local version)'
             end if
          else
             write(iul,'(A,I3)') ' using local histograms'
          end if
       endif
       if (removeMode == 1) then
          call getSampleEnergyAndVariance(smpl,Emean,var,sigma)
          call removeSampleOutliers(smpl,Emean,var,fac,rfac,riter,nbins)
       else if (removeMode == 2) then
          call removeSampleOutliers2(smpl,tol,global,replace)
       else
          call abortp("illegal removeMode in sample")
       end if
    case (6)
       call unWeight(smpl)
    case default
       call abortp('sample: illegal mSampleMode')
    end select

    ! finalize
    if (start == 'lmo') then
       call destroy_LMOSampling()
    endif

    logmode = oldLogmode

  end subroutine sample


  !-------------------------------------!
  subroutine initInitialWalker(lines,nl)
  !-------------------------------------!

    integer, intent(in) :: nl
    character(len=120), intent(in) :: lines(nl)
    integer i,j,a
    real(r8) :: scale
    character(len=4) :: s
    character(len=:), allocatable :: words(:)
    real(r8) :: r(3)

    call assert(ne > 0,"initInitialWalker: wf is not yet read")
    if (ALLOCATED(mX)) deallocate(mX)
    if (ALLOCATED(mY)) deallocate(mY)
    if (ALLOCATED(mZ)) deallocate(mZ)

    ! read $init_walker section: read initial walker
    allocate(mX(ne),mY(ne),mZ(ne))
    read(lines(2),'(A4)') s
    if (TRIM(s) == 'col') then
       do i=1,ne
          read(lines(i+2),*) mX(i),mY(i),mZ(i)
       enddo
    else if (s == 'free') then
       do i=1,ne
          words = Py_split(lines(i+2))
          if (TRIM(words(1)) == 'between') then
             r = 0._r8
             do j=2,SIZE(words)/LEN(words)
                read(words(j),*) a
                r = r + atoms(a)%Get_position()
             end do
             r = r / (SIZE(words)/LEN(words) - 1)
          else if (TRIM(words(1)) == 'scaled') then
             read(words(4),*) scale
             if (asserts) call assert(scale >= 0._r8)
             if (asserts) call assert(scale <= 1._r8)
             read(words(2),*) a
             r = (1._r8 - scale)*atoms(a)%Get_position()
             read(words(3),*) a
             r = r + scale * atoms(a)%Get_position()
          else if (TRIM(words(1)) == 'at') then
             read(words(2),*) a
             r = atoms(a)%Get_position()
          else
             ! read like col format
             if (asserts) call assert(SIZE(words)/LEN(words) == 3)
             do j=1,SIZE(words)/LEN(words)
                read(words(j),*) r(j)
             end do
          end if
          mX(i) = r(1); mY(i) = r(2); mZ(i) = r(3)
       enddo
    else
       read(lines(2),*) (mX(i),i=1,ne)
       read(lines(3),*) (mY(i),i=1,ne)
       read(lines(4),*) (mZ(i),i=1,ne)
    endif
    if (logmode>=2) then
       write(iul,*) ' setting initial walker position:'
       do i=1,ne
          write(iul,'(3F9.5)') mX(i), mY(i), mZ(i)
       enddo
    endif

  end subroutine initInitialWalker


  !-------------------------------------!
  subroutine displaceInitialWalker(lines,nl)
  !-------------------------------------!

    integer, intent(in) :: nl
    character(len=120), intent(in) :: lines(nl)
    character(len=1) :: s
    real(r8) :: d
    integer :: i, iflag

    call assert(ALLOCATED(mX), "displaceInitialWalker: initial walker is not initialized, call init_walker first.")

    call getdbla(lines, nl, "displacement=", d, iflag)
    call assert(iflag==0, "displaceInitialWalker: displacement has to be given.")

    call getstra(lines, nl, "dimension=", s, iflag)
    call assert(iflag==0, "displaceInitialWalker: dimension has to be given.")

    if (s=='x') then
      do i=1,ne
        mX(i) = mX(i) + d
      enddo
    else if (s=='y') then
      do i=1,ne
        mY(i) = mY(i) + d
      enddo
    else if (s=='z') then
      do i=1,ne
        mZ(i) = mZ(i) + d
      enddo
    else
      call abortp('displaceInitialWalker: dimension has to be x,y or z')
    end if

    if (logmode>=2) then
       write(iul,*) ' displacement: ', d
       write(iul,*) ' dimension: ', s
       write(iul,*) ' displacing initial walker position to:'
       do i=1,ne
          write(iul,'(3F9.5)') mX(i), mY(i), mZ(i)
       enddo
    endif

  end subroutine displaceInitialWalker


  integer function initSampleSize()
  !-------------------------------!
    initSampleSize = mInitSampleSize
  end function initSampleSize


  integer function currentSampleSize()
  !----------------------------------!
    ! returns the number that the current sample is supposed to be (i.e. the
    ! size that was set with size= or new_size=)
    currentSampleSize = mCurrentSampleSize
  end function currentSampleSize


  subroutine getInitialSample(sample,maxsize)
  !-----------------------------------------!
    type(RWSample), intent(inout) :: sample
    integer, intent(in)           :: maxsize

    call assert(mInitSampleSize>0 .and. maxsize>0,'getInitialSample: sample size not set')
    if (.not. isEmpty(sample)) call destroySample(sample)
    call initSample(sample,maxsize)

    select case(mSampleType)
    case (1)
       call sample1(sample)
    case (2)
       call sample2(sample)
    case (3)
       call sample3(sample)
    case (4)
       call sample4(sample)
    case default
       call abortp('getInitialSample: readChkPoint not yet implemented')
       !call rwReadChkpoint(sample)
       !mytid = 0
       !!MPI call mpi_comm_rank(MPI_COMM_WORLD, mytid, ierr)
       !seed = seed0
       !if (seed < 0) call abortp (' seed must greater than 0')
       !seed = seed + mytid ! add task ID for parallel runs
       !seed = -seed        ! i.e. initialize at first call
       !dummy = init_mrg(seed)  ! initizalize mrg
       !if (discd1/=0) then
       !   taueff = -1.d0
       !endif
    end select

  end subroutine getInitialSample



  subroutine readInitialSample(sample,maxsize,posfilename)
  !------------------------------------------------------!
    type(RWSample), intent(inout) :: sample
    integer, intent(in)           :: maxsize
    character(len=80), intent(in) :: posfilename
    integer io,ssize
    integer, parameter :: iu=19

    call assert(maxsize>0,'readInitialSample: max sample size not set')
    if (.not. isEmpty(sample)) call destroySample(sample)
    call initSample(sample,maxsize)
    if (MASTER) then
       open(iu,file=posfilename,form='unformatted',status='old',iostat=io)
       if (io /= 0) call abortp('$sample: position file '//trim(posfilename)//' could not be opened')
    end if
    if (mInitSampleSize==0) then
       ssize = 0
       call readSamplePositions(sample,iu,ssize)
       mInitSampleSize = ssize
    else
       ssize = mInitSampleSize
       call readSamplePositions(sample,iu,ssize)
    end if

    if (MASTER) then
       close(iu)
    end if
  end subroutine readInitialSample


  !------------------------!
  subroutine sample1(sample)
  !------------------------!

! VMC calculation with one walker using time step mSampleTau
! and block size mBlockSize. Discard mSampleDiscard blocks
! add walker at block end to sample until mSampleSize is reached

    integer, parameter :: BufSize = 20            ! ring buffer size for current acc.ratio
    integer, parameter :: MaxSteps = 10000
    integer, parameter :: MaxBlocks = 100000
    real(r8), parameter  :: MaxTau = 1.0d0
    real(r8), parameter  :: MinTau = 0.000005d0
    real(r8), parameter  :: TauScal = 0.3d0         ! scaling factor for tau adaptation

    type(RWSample), intent(inout) :: sample
    type(RandomWalker) :: rw
    type(RandomWalker),pointer :: rwpp
    type(RandomWalker),pointer :: rwp(:)
    integer :: step,bufStep,block,mode,n,bs,iull
    real(r8)  :: accbuf(BufSize)      ! ring buffer for recent acceptance ratio
    real(r8)  :: blockTime            ! actual total time in a block
    real(r8)  :: targetBlockTime      ! desired blockTime
    real(r8)  :: tau                  ! actual time step
    real(r8)  :: eloc,currentAccRatio,eRef,accepted,lastAcceptedSteps
    real(r8)  :: accRatio
    integer :: ctpersist

    call propagator_reset()
    if (mSampleTau <= 0.0d0) then
       mSampleTau = propagator_getInitialTau(mTargetAccRatio)
    endif
    call propagator_init(0,mMove,mTmove,mBlockDiscard,0,mSampleTau,1.d0,.true.,mRejCross,mMoveType,1)

    mode = mWalkerGenMode
    ! use electron position of init_walker section if available
    if (.not. allocated(mX)) then
       allocate(mX(ne),mY(ne),mZ(ne))
       call createRandomElectronPositions(mode,mX,mY,mZ)
    end if

    call rw_new(rw)
    call setTo(rw,mX,mY,mZ)

    call appendWalker(sample,rw)
    bs = 1
    call getFirstWalkerBlock(sample,bs,rwp)    ! random walk with 1st walker of sample (not a copy!)
    if (logmode >= 5) call display(rwp(1),2,iul)
    call resetTo(rwp(1),mX,mY,mZ)          ! update 1st walker
    if (logmode >= 5) call display(rwp(1),2,iul)

    ! do some initial steps without accept/reject step
    ! to avoid initial persistence
    eRef = 0.d0
    accRatio = 0.0d0
    mode = 2     ! no reject
    do n=1,10
       do step=1,5
          call propagateAndWeight(rwp,eRef,mode,accepted)
       enddo
       call getAccRatio(accRatio)
       if (accRatio > 0.01d0) exit
    enddo
    if (n==11) call abortp('sample1: generation of starting vmc walker failed')

    if (logmode >= 5) then
       write(iul,*) "start of sample VMC:"
       call display(rwp(1),2,iul)
    endif

    targetBlockTime = mSampleTau * mBlockSize / mTargetAccRatio
    tau = propagator_timeStep()
    mode = 1    ! with reject step

    ! do BufSize steps to fill ring-buffer until acc ratio > 0.2
    do n=1,10
       lastAcceptedSteps = 0
       do step = 1,BufSize
          call propagateAndWeight(rwp,eRef,mode,accepted)
          accbuf(step) = accepted
          lastAcceptedSteps = lastAcceptedSteps + accepted
       enddo
    if (dble(lastAcceptedSteps)/BufSize >  0.2) exit
       tau = max(tau*0.2d0,minTau)
       write(123,*) 'ringbuffer:',tau,dble(lastAcceptedSteps)/BufSize
       call propagator_setTimeStep(tau)
    enddo
    if (n==11) call abortp('sample1: generation of accept ratio > 0.2 failed')

    ctpersist = 0

    do block=1,MaxBlocks

       blockTime = 0
       do step=1,MaxSteps

          call propagateAndWeight(rwp,eRef,mode,accepted)
          eloc = E_local(rwp(1))
          blockTime = blockTime + tau

          if (logmode>=5) then
             write(iul,'(2(A,I5),4(A,F10.6))') 'sample1: block=',block,' step=',step, &
                                               ' tau=,',tau,' blockTime=',blockTime, &
                                               ' E_L=',eloc,' acc=',accepted
             call display(rwp(1),2,iul)
          endif

          bufStep = mod(step,BufSize)+1
          lastAcceptedSteps = lastAcceptedSteps - accbuf(bufStep) + accepted
          accbuf(bufStep) = accepted
          currentAccRatio = lastAcceptedSteps / BufSize
!
!          ! adapt tau to reach target acceptance ratio
!          if (currentAccRatio > mTargetAccRatio) then
!             tau = tau * (1 + TauScal*log(currentAccRatio/mTargetAccRatio))
!          else
!             tau = tau / (1 - TauScal*log(currentAccRatio/mTargetAccRatio))
!          endif
!          tau = max(tau,MinTau)
!          tau = min(tau,MaxTau)
!          write(123,*) 'setting tau:',tau,currentAccRatio,mTargetAccRatio,TauScal*log(currentAccRatio/mTargetAccRatio)
!          call propagator_setTimeStep(tau)

          if (blockTime > targetBlockTime) exit
       enddo

       if (block > mBlockDiscard .and. persist(rwp(1)) < 5) call appendWalker(sample,rwp(1))
       if (persist(rwp(1)) >= 5) ctpersist = ctpersist + 1
       if (getSampleSize(sample) >= mInitSampleSize) exit

    enddo

    if (logmode >= 5) then
       write(103,*) "Check Persistency of the Sample"
       rwpp => getFirst(sample)
       do
          write(103,*) "Persistency",persist(rwp(1))
       if (.not.isNext(sample)) exit
          rwpp => getNext(sample)
       enddo
       write(103,*)
    endif

    if (MASTER .and. logmode >= 2) then
       write(iul,*) ' initial sample created with VMC'
       if (mWalkerGenMode==1) then
          write(iul,*) ' with start=gaussian'
       else if (mWalkerGenMode==2) then
          write(iul,*) ' with start=density'
       endif
       write(iul,'(A,I8)')   ' approx steps per block:      ',int(mBlockSize / mTargetAccRatio)
       write(iul,'(A,I8)')   ' # of blocks:                 ',block - 1
       write(iul,'(A,F8.4)') ' initial time step:           ',mSampleTau
       write(iul,'(A,F8.4)') ' final time step:             ',propagator_timeStep()
       write(iul,'(A,F8.2)') ' final acceptance ratio:      ',currentAccRatio
       write(iul,'(A,I8)')   ' discarded persistent walker: ',ctpersist
       write(iul,*)
    endif
    if (logmode >= 3) then
       iull = 80 + mytid
       write(iull,*)
       write(iull,*) ' initial sample created with VMC'
       write(iull,'(A,I8)')   ' sample size now:             ',getSampleSize(sample)
       write(iull,'(A,I8)')   ' approx steps per block:      ',int(mBlockSize / mTargetAccRatio)
       write(iull,'(A,I8)')   ' # of blocks:                 ',block - 1
       write(iull,'(A,F8.4)') ' initial time step:           ',mSampleTau
       write(iull,'(A,F8.4)') ' final time step:             ',propagator_timeStep()
       write(iull,'(A,F8.2)') ' final acceptance ratio:      ',currentAccRatio
       write(iull,'(A,I8)')   ' discarded persistent walker: ',ctpersist
       write(iull,*)
    endif

  end subroutine sample1



  !-------------------------
  subroutine sample2(sample)
  !-------------------------

    ! create sample from random electron positions

    type(RWSample), intent(inout) :: sample
    real(r8), allocatable :: x(:), y(:), z(:)
    type(RandomWalker) :: rw
    integer mode

    allocate(x(ne),y(ne),z(ne))

    mode = mWalkerGenMode
    if (logmode>=4) write(iul,*) 'sample2 with mode = ',mode

    ! use initial walker if it exists
    call rw_new(rw)
    if (allocated(mX)) then
       call setTo(rw,mX,mY,mZ)
    else
       call createRandomElectronPositions(mode,x,y,z)
       call setTo(rw,x,y,z)
    endif

    do
       call appendWalker(sample,rw)
       if (logmode>=4) write(iul,*) 'sample2 ',getSampleSize(sample),E_local(rw)
    if (getSampleSize(sample) == mInitSampleSize) exit
       call createRandomElectronPositions(mode,x,y,z)
       call resetTo(rw,x,y,z)
    enddo

    if (MASTER .and. logmode >= 2) then
       select case (mWalkerGenMode)
       case (1)
          write(iul,'(/A)') ' created initial gaussian distributed random sample'
       case (2)
          write(iul,'(/A)') ' created initial random sample using atom densities'
       case (3)
          write(iul,'(/A)') ' created initial random sample using approximate lmo sampling'
       case default
          call abortp('sample2: illegal WalkerGenMode')
       end select
    endif

  end subroutine sample2



  !-------------------------
  subroutine sample3(sample)
  !-------------------------

    ! create sample from random electron positions
    ! with additional energy treshold for accepted positions

    integer, parameter :: maxSteps=100000
    type(RWSample), intent(inout) :: sample
    integer             :: i,mode
    real(r8), allocatable :: x(:), y(:), z(:)
    real(r8)              :: E
    type(simpleStat)    :: stat
    type(RandomWalker)  :: rw

    allocate(x(ne),y(ne),z(ne))
    call reset(stat)

    mode = mWalkerGenMode
    call createRandomElectronPositions(mode,x,y,z)

    call rw_new(rw)
    call setTo(rw,x,y,z)

    do i=1,maxSteps
       E = E_local(rw)
       if (E > mEmin .and. E < mEmax) then
          call appendWalker(sample,rw)
          call addData(stat,E)
       endif
       if (getSampleSize(sample) == mInitSampleSize) exit
       call createRandomElectronPositions(mode,x,y,z)
       call resetTo(rw,x,y,z)
    enddo

    if (MASTER .and. logmode >= 2) then
       select case (mWalkerGenMode)
       case (1)
          write(iul,'(/A)') ' created initial gaussian distributed random sample'
       case (2)
          write(iul,'(/A)') ' created initial random sample using atom densities'
       case (3)
          write(iul,'(/A)') ' created initial random sample using approximate lmo sampling'
       case default
          call abortp('sample3: illegal WalkerGenMode')
       end select
       write(iul,'(2(A,F13.5))') ' Emin=',minValue(stat),' Emax=',maxValue(stat)
       write(iul,'(a,i5,a)') ' with ',i-getSampleSize(sample),' failed attempts using thresholds'
       write(iul,'(2(a,f13.3))') ' Emin=',mEmin,' Emax=',mEmax
    endif

  end subroutine sample3



  !-------------------------
  subroutine sample4(sample)
  !-------------------------
    ! create sample with identical electron positions, either from $init_walker input
    ! or from random electron positions (using gaussian distributions)

    type(RWSample), intent(inout) :: sample
    integer             :: i,mode
    type(RandomWalker)  :: rw

    ! use electron position of init_walker section if available
    mode = mWalkerGenMode
    if (.not.allocated(mX)) then
       allocate(mX(ne),mY(ne),mZ(ne))
       call createRandomElectronPositions(mode,mX,mY,mZ)
    end if

    call rw_new(rw)
    call setTo(rw,mX,mY,mZ)

    do i=1,mInitSampleSize
       call appendWalker(sample,rw)
    enddo

    if (MASTER .and. logmode >= 2) then
       write(iul,'(/A,I5,A)') ' created initial identical sample (generate=single_point)'
       if (logmode >= 3) then
          write(iul,*) ' at position:'
          call display(rw,2,iul)
       endif
       write(iul,*)
    endif

  end subroutine sample4




end module qmcSample_m
