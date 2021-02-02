! Copyright (C) 2006-2007, 2012, 2015-2016 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

! simple statistics module to provide data types
! simplifiying statistical analysis

! Contains the modules/data types
!  simple_statistics/simpleStat
!  weight_statistics/weightStat
!  block_statistics/blockStat
!  int_histogram_statistics/IntHistogram
!  double_histogram_statistics/IndexedHistogram

! all modules are accessible thru generic interface module 'statistics'
! (at the end of file), allowing generic names for all types
! e.g. addData(stat,data) and mean(stat) works with all stat/histogram types

! this version has min/max added to simpleStat. AL 17.2.06
! this version has histogram types added AL 1.2.07


!=======================
module simple_statistics
!=======================

  use kinds_m, only: r8, i8
  use error_m
  implicit none

  type simpleStat
    private
    real(r8)    :: averageSum = 0
    real(r8)    :: varianceSum = 0
    integer(i8) :: dataCounter = 0
    real(r8)    :: minValue = 0
    real(r8)    :: maxValue = 0
  end type simpleStat

  interface assignment(=)
     module procedure assign_simpleStat
  end interface

  interface operator(+)
     module procedure add_simpleStat
  end interface

contains

! -----------------------------
  subroutine ss_addData(self,data)
! -----------------------------

    type(simpleStat), intent(inout) :: self
    real(r8), intent(in)            :: data

    self%averageSum = self%averageSum + data
    self%varianceSum = self%varianceSum + data**2
    self%dataCounter = self%dataCounter + 1
    if (self%dataCounter==1) then
       self%minValue = data
       self%maxValue = data
    else
       self%minValue = min(data,self%minValue)
       self%maxValue = max(data,self%maxValue)
    end if

  end subroutine ss_addData


! ----------------------
  subroutine ss_reset(self)
! ----------------------

    type(simpleStat), intent(inout) :: self

    self%averageSum = 0
    self%varianceSum = 0
    self%dataCounter = 0
    self%minValue = 0
    self%maxValue = 0

  end subroutine ss_reset


! ---------------------------------------
  subroutine assign_simpleStat(self,stat)
! ---------------------------------------

    type(simpleStat), intent(out) :: self
    type(simpleStat), intent(in)  :: stat

    self%averageSum = stat%averageSum
    self%varianceSum = stat%varianceSum
    self%dataCounter = stat%dataCounter
    self%minValue = stat%minValue
    self%maxValue = stat%maxValue

  end subroutine assign_simpleStat


! --------------------------------------------------
  type(simpleStat) function add_simpleStat(self,add)
! --------------------------------------------------

    type(simpleStat), intent(in)  :: self
    type(simpleStat), intent(in)  :: add

    add_simpleStat%averageSum = self%averageSum + add%averageSum
    add_simpleStat%varianceSum = self%varianceSum + add%varianceSum
    add_simpleStat%dataCounter = self%dataCounter + add%dataCounter
    add_simpleStat%minValue = min(add%minValue,self%minValue)
    add_simpleStat%maxValue = max(add%maxValue,self%maxValue)

  end function add_simpleStat


! access functions


! ----------------------------
  real(r8) function ss_mean(self)
! ----------------------------

    type(simpleStat), intent(in)  :: self

    call assert(self%dataCounter > 0,"simpleStat:mean")

    ss_mean = self%averageSum / self%dataCounter

  end function ss_mean


! -------------------------------------
  real(r8) function ss_meanAllNodes(self)
! -------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(simpleStat), intent(in)  :: self
    real(r8) sendbuf(2),recvbuf(2)
#ifdef MPI
    integer ierr
#endif

    sendbuf(1) = self%averageSum; sendbuf(2) = self%dataCounter
    recvbuf = sendbuf  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(sendbuf,recvbuf,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    if (recvbuf(2) > 0) then
       ss_meanAllNodes = recvbuf(1) / recvbuf(2)
    else
       ss_meanAllNodes = 0
    endif

  end function ss_meanAllNodes


! --------------------------------
  real(r8) function ss_variance(self)
! --------------------------------

    type(simpleStat), intent(in)  :: self

    call assert(self%dataCounter > 0,"simpleStat:variance")

    ss_variance = self%varianceSum/self%dataCounter - ss_mean(self)**2

  end function ss_variance



! -----------------------------------------
  real(r8) function ss_varianceAllNodes(self)
! -----------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(simpleStat), intent(in)  :: self
    real(r8) sendbuf(3),recvbuf(3)
#ifdef MPI
    integer ierr
#endif

    sendbuf(1) = self%averageSum; sendbuf(2)=self%varianceSum; sendbuf(3) = self%dataCounter
    recvbuf = sendbuf  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(sendbuf,recvbuf,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    if (recvbuf(3) > 0) then
       ss_varianceAllNodes = recvbuf(2)/recvbuf(3)-(recvbuf(1)/recvbuf(3))**2
    else
       ss_varianceAllNodes = 0
    endif

  end function ss_varianceAllNodes


! ----------------------------------
  real(r8) function ss_stdDevMean(self)
! ----------------------------------

    ! returns the standard deviation of the mean

    type(simpleStat), intent(in)  :: self

    call assert(self%dataCounter > 1,"simpleStat:stdDevMean")

    ss_stdDevMean = sqrt(ss_variance(self) / (self%dataCounter - 1))

  end function ss_stdDevMean


! -------------------------------------------
  real(r8) function ss_stdDevMeanAllNodes(self)
! -------------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(simpleStat), intent(in)  :: self
    real(r8) sendbuf(3),recvbuf(3)
#ifdef MPI
    integer ierr
#endif

    sendbuf(1) = self%averageSum; sendbuf(2)=self%varianceSum; sendbuf(3) = self%dataCounter
    recvbuf = sendbuf  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(sendbuf,recvbuf,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    if (recvbuf(3) > 1) then
       ss_stdDevMeanAllNodes = sqrt((recvbuf(2)/recvbuf(3)-(recvbuf(1)/recvbuf(3))**2)/(recvbuf(3)-1.d0))
    else
       ss_stdDevMeanAllNodes = 0
    endif

  end function ss_stdDevMeanAllNodes


! ------------------------------------
  integer(i8) function ss_dataCount(self)
! ------------------------------------

    type(simpleStat), intent(in)  :: self

    ss_dataCount = self%dataCounter

  end function ss_dataCount


! ------------------------------------------
  integer(i8) function ss_dataCountAllNodes(self)
! ------------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(simpleStat), intent(in)  :: self
    integer(i8) send,recv
#ifdef MPI
    integer ierr
#endif

    send = self%dataCounter
    recv = send  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(send,recv,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    ss_dataCountAllNodes = recv

  end function ss_dataCountAllNodes


! ------------------------------
  real(r8) function ss_minValue(self)
! ------------------------------

    type(simpleStat), intent(in)  :: self

    call assert(self%dataCounter > 0,"simpleStat:minValue")

    ss_minValue = self%minValue

  end function ss_minValue

! ------------------------------
  real(r8) function ss_maxValue(self)
! ------------------------------

    type(simpleStat), intent(in)  :: self

    call assert(self%dataCounter > 0,"simpleStat:maxValue")

    ss_maxValue = self%maxValue

  end function ss_maxValue



  !----------------------------
  subroutine ss_serialize(self,iu)
  !----------------------------

    ! write object formatted to open file unit iu

    type(simpleStat), intent(in)  :: self
    integer, intent(in)           :: iu    ! open file unit

    write(iu,*) self%averageSum,self%varianceSum,self%dataCounter, &
                self%minValue,self%maxValue

  end subroutine ss_serialize


  !------------------------------
  subroutine ss_deSerialize(self,iu)
  !------------------------------

    ! write object formatted to open file unit iu

    type(simpleStat), intent(inout)  :: self
    integer, intent(in)              :: iu    ! open file unit

    read(iu,*) self%averageSum,self%varianceSum,self%dataCounter, &
               self%minValue,self%maxValue

  end subroutine ss_deSerialize


  !--------------------------
  subroutine ss_debugPrint(self)
  !--------------------------

    type(simpleStat), intent(in)  :: self

    print*,"debugPrint:averageSum=",self%averageSum
    print*,"debugPrint:varianceSum=",self%varianceSum
    print*,"debugPrint:dataCounter=",self%dataCounter

  end subroutine ss_debugPrint

end module simple_statistics




!=============================================================================


!=======================
module weight_statistics
!=======================

  ! weight_statistics is compatible with simple_statistics
  ! add_data uses optional weight argument
  use kinds_m, only: r8, i8
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
  use error_m
  implicit none

  type weightStat
    private
    real(r8)    :: averageSum = 0
    real(r8)    :: varianceSum = 0
    real(r8)    :: weightSum = 0
    integer(i8) :: dataCounter = 0
    real(r8)    :: minValue = 0
    real(r8)    :: maxValue = 0
  end type weightStat

  interface assignment(=)
     module procedure assign_weightStat
  end interface

  interface operator(+)
     module procedure add_weightStat
  end interface

contains

! ------------------------------------
  subroutine ws_addData(self,data,weight)
! ------------------------------------

    type(weightStat), intent(inout) :: self
    real(r8), intent(in)            :: data
    real(r8), optional,intent(in)   :: weight

    real(r8) :: w

    if (.not.present(weight)) then
       w = 1
    else
       w = weight
    endif

    self%averageSum = self%averageSum + w*data
    self%varianceSum = self%varianceSum + w*data**2
    self%weightSum = self%weightSum + w
    self%dataCounter = self%dataCounter + 1

    if (self%dataCounter==1) then
       self%minValue = data
       self%maxValue = data
    else
       self%minValue = min(data,self%minValue)
       self%maxValue = max(data,self%maxValue)
    end if
  end subroutine ws_addData


! ----------------------
  subroutine ws_reset(self)
! ----------------------

    type(weightStat), intent(inout) :: self

    self%averageSum = 0
    self%varianceSum = 0
    self%weightSum = 0
    self%dataCounter = 0
    self%minValue = 0
    self%maxValue = 0

  end subroutine ws_reset


! ---------------------------------------
  subroutine assign_weightStat(self,stat)
! ---------------------------------------

    type(weightStat), intent(out) :: self
    type(weightStat), intent(in)  :: stat

    self%averageSum = stat%averageSum
    self%varianceSum = stat%varianceSum
    self%weightSum = stat%weightSum
    self%dataCounter = stat%dataCounter
    self%minValue = stat%minValue
    self%maxValue = stat%maxValue

  end subroutine assign_weightStat


! --------------------------------------------------
  type(weightStat) function add_weightStat(self,add)
! --------------------------------------------------

    type(weightStat), intent(in)  :: self
    type(weightStat), intent(in)  :: add

    add_weightStat%averageSum = self%averageSum + add%averageSum
    add_weightStat%varianceSum = self%varianceSum + add%varianceSum
    add_weightStat%weightSum = self%weightSum + add%weightSum
    add_weightStat%dataCounter = self%dataCounter + add%dataCounter
    add_weightStat%minValue = min(add%minValue,self%minValue)
    add_weightStat%maxValue = max(add%maxValue,self%maxValue)

  end function add_weightStat


! access functions


! ----------------------------
  real(r8) function ws_mean(self)
! ----------------------------

    type(weightStat), intent(in)  :: self

    call assert(self%dataCounter > 0,"weightStat:mean")

    ws_mean = self%averageSum / self%weightSum

  end function ws_mean


! -------------------------------------
  real(r8) function ws_meanAllNodes(self)
! -------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(weightStat), intent(in)  :: self
    real(r8) sendbuf(2),recvbuf(2)
#ifdef MPI
    integer ierr
#endif

    sendbuf(1) = self%averageSum; sendbuf(2) = self%weightSum
    recvbuf = sendbuf  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(sendbuf,recvbuf,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    if (recvbuf(2) > 0) then
       ws_meanAllNodes = recvbuf(1) / recvbuf(2)
    else
       ws_meanAllNodes = 0
    endif

  end function ws_meanAllNodes


! --------------------------------
  real(r8) function ws_variance(self)
! --------------------------------

    type(weightStat), intent(in)  :: self

    call assert(self%weightSum > 0,"weightStat:variance")

    ws_variance = self%varianceSum/self%weightSum - ws_mean(self)**2

  end function ws_variance


! -----------------------------------------
  real(r8) function ws_varianceAllNodes(self)
! -----------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(weightStat), intent(in)  :: self
    real(r8) sendbuf(3),recvbuf(3)
#ifdef MPI
    integer ierr
#endif

    sendbuf(1) = self%averageSum; sendbuf(2)=self%varianceSum; sendbuf(3) = self%weightSum
    recvbuf = sendbuf  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(sendbuf,recvbuf,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    if (recvbuf(3) > 0) then
       ws_varianceAllNodes = recvbuf(2)/recvbuf(3)-(recvbuf(1)/recvbuf(3))**2
    else
       ws_varianceAllNodes = 0
    endif

  end function ws_varianceAllNodes


! ----------------------------------
  real(r8) function ws_stdDevMean(self)
! ----------------------------------

    ! returns the standard deviation of the mean
    ! careful: returns sensible values only if independent data !!!

    type(weightStat), intent(in)  :: self

    if (self%dataCounter == 1) then
      ws_stdDevMean = ieee_value(ws_stdDevMean, ieee_positive_inf)
    else
      ws_stdDevMean = sqrt(ws_variance(self) / (self%dataCounter - 1))
    end if

  end function ws_stdDevMean


! -------------------------------------------
  real(r8) function ws_stdDevMeanAllNodes(self)
! -------------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(weightStat), intent(in)  :: self
    real(r8) sendbuf(4),recvbuf(4)
#ifdef MPI
    integer ierr
#endif

    sendbuf(1) = self%averageSum; sendbuf(2)=self%varianceSum
    sendbuf(3) = self%weightSum; sendbuf(4)=self%dataCounter
    recvbuf = sendbuf  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(sendbuf,recvbuf,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    if (recvbuf(4) > 1) then
       ws_stdDevMeanAllNodes = sqrt( (recvbuf(2)/recvbuf(3)-(recvbuf(1)/recvbuf(3))**2)/(recvbuf(4)-1.d0) )
    else
       ws_stdDevMeanAllNodes = 0
    endif

  end function ws_stdDevMeanAllNodes


! ------------------------------------
  integer(i8) function ws_dataCount(self)
! ------------------------------------

    type(weightStat), intent(in)  :: self

    ws_dataCount = self%dataCounter

  end function ws_dataCount


! ------------------------------------------
  integer(i8) function ws_dataCountAllNodes(self)
! ------------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(weightStat), intent(in)  :: self
    integer(i8) send,recv
#ifdef MPI
    integer ierr
#endif

    send = self%dataCounter
    recv = send  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(send,recv,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    ws_dataCountAllNodes = recv

  end function ws_dataCountAllNodes


! ------------------------------------
  real(r8) function ws_totalWeight(self)
! ------------------------------------

    type(weightStat), intent(in)  :: self

    ws_totalWeight = self%weightSum

  end function ws_totalWeight


! --------------------------------------------
  real(r8) function ws_totalWeightAllNodes(self)
! --------------------------------------------

#ifdef MPI
    use MPI_F08
#endif

    type(weightStat), intent(in)  :: self
    real(r8) send,recv
#ifdef MPI
    integer ierr
#endif

    send = self%dataCounter
    recv = send  ! default for serial run

#ifdef MPI
    call MPI_ALLREDUCE(send,recv,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    ws_totalWeightAllNodes = recv

  end function ws_totalWeightAllNodes


! ------------------------------
  real(r8) function ws_minValue(self)
! ------------------------------

    type(weightStat), intent(in)  :: self

    call assert(self%dataCounter > 0,"weightStat:minValue")

    ws_minValue = self%minValue

  end function ws_minValue

! ------------------------------
  real(r8) function ws_maxValue(self)
! ------------------------------

    type(weightStat), intent(in)  :: self

    call assert(self%dataCounter > 0,"weightStat:maxValue")

    ws_maxValue = self%maxValue

  end function ws_maxValue



  !----------------------------
  subroutine ws_serialize(self,iu)
  !----------------------------

    ! write object formatted to open file unit iu

    type(weightStat), intent(in)  :: self
    integer, intent(in)           :: iu    ! open file unit

    write(iu,*) self%averageSum,self%varianceSum,self%weightSum, &
                self%dataCounter,self%minValue,self%maxValue

  end subroutine ws_serialize


  !------------------------------
  subroutine ws_deSerialize(self,iu)
  !------------------------------

    ! write object formatted to open file unit iu

    type(weightStat), intent(inout)  :: self
    integer, intent(in)              :: iu    ! open file unit

    read(iu,*) self%averageSum,self%varianceSum,self%weightSum, &
               self%dataCounter,self%minValue,self%maxValue

  end subroutine ws_deSerialize


  !---------------------------
  subroutine ws_debugPrint(self)
  !---------------------------

    type(weightStat), intent(in)  :: self

    print*,"debugPrint:averageSum=",self%averageSum
    print*,"debugPrint:varianceSum=",self%varianceSum
    print*,"debugPrint:dataCounter=",self%dataCounter

  end subroutine ws_debugPrint

end module weight_statistics



!=========================================================


!======================
module block_statistics
!======================

  use error_m
  use simple_statistics
  use weight_statistics
  implicit none

  type blockStat
     private
     type(simpleStat) :: allBlockStat
     type(weightStat) :: thisBlockStat
     type(weightStat) :: totalWeightStat
     integer          :: blockLength = 1
     integer          :: blockCounter = 1
     integer          :: discard = 0
     logical          :: blockEnd = .false.
  end type blockStat

  interface assignment(=)
     module procedure assign_blockStat
  end interface

  interface operator(+)
     module procedure add_blockStat
  end interface

contains

  !----------------------------------------
  subroutine bs_init(self,blockLength,discard)
  !----------------------------------------

    type(blockStat), intent(inout) :: self
    integer, optional, intent(in)  :: blockLength
    integer, optional, intent(in)  :: discard

    call ss_reset(self%allBlockStat)
    call ws_reset(self%thisBlockStat)
    call ws_reset(self%totalWeightStat)

    if (present(blockLength)) then
       call assert(blockLength > 0,"blockStat:init::blockLength>0")
       self%blockLength = blockLength
    else
       self%blockLength = 1
    endif

    self%blockCounter = 1

    if (present(discard)) then
       call assert(discard > 0,"blockStat:init::discard>0")
       self%discard = discard
    else
       self%discard = 0
    endif

    self%blockEnd = .false.

  end subroutine bs_init


  !-----------------------------------
  subroutine bs_addData(self,data,weight)
  !-----------------------------------

    type(blockStat), intent(inout) :: self
    real(r8), intent(in)            :: data
    real(r8), optional,intent(in)   :: weight

    real(r8) :: w

    if (present(weight)) then
       w = weight
    else
       w = 1
    endif

    self%blockEnd = .false.

    call ws_addData(self%thisBlockStat,data,w)

    if (self%blockCounter >= self%discard) then
       call ws_addData(self%totalWeightStat,data,w)
    endif

    if (ws_dataCount(self%thisBlockStat) == self%blockLength) then
       if (self%blockCounter > self%discard) then
          call ss_addData(self%allBlockStat,ws_mean(self%thisBlockStat))
       endif
       call ws_reset(self%thisBlockStat)
       self%blockEnd = .true.
       self%blockCounter = self%blockCounter + 1
    endif

  end subroutine bs_addData

  !---------------------
  subroutine bs_reset(self)
  !---------------------

    type(blockStat), intent(inout) :: self

    call ss_reset(self%allBlockStat)
    call ws_reset(self%thisBlockStat)
    call ws_reset(self%totalWeightStat)

    self%blockCounter = 1

    self%blockEnd = .false.

  end subroutine bs_reset


! --------------------------------------
  subroutine assign_blockStat(self,stat)
! --------------------------------------

    type(blockStat), intent(out) :: self
    type(blockStat), intent(in)  :: stat

    self%allBlockStat = stat%allBlockStat
    self%thisBlockStat = stat%thisBlockStat
    self%totalWeightStat = stat%totalWeightStat
    self%blockLength = stat%blockLength
    self%blockCounter = stat%blockCounter
    self%discard = stat%discard
    self%blockEnd = stat%blockEnd

  end subroutine assign_blockStat


! ------------------------------------------------
  type(blockStat) function add_blockStat(self,add)
! ------------------------------------------------

    type(blockStat), intent(in)  :: self
    type(blockStat), intent(in)  :: add

    call assert(self%blockLength == add%blockLength, &
         "blockStat:add:: block lengths must be equal")
    call assert(self%discard == add%discard, &
         "blockStat:add:: discard must be equal")


    ! TODO

    call bs_reset(add_blockStat)

  end function add_blockStat


! access functions


! ----------------------------
  real(r8) function bs_mean(self)
! ----------------------------

    type(blockStat), intent(in)  :: self

    bs_mean = ws_mean(self%totalWeightStat)

  end function bs_mean


! ---------------------------------
  real(r8) function bs_blockMean(self)
! ---------------------------------

    type(blockStat), intent(in)  :: self

    bs_blockMean = ss_mean(self%allBlockStat)

  end function bs_blockMean


! --------------------------------
  real(r8) function bs_variance(self)
! --------------------------------

    type(blockStat), intent(in)  :: self

    bs_variance = ws_variance(self%totalWeightStat)

  end function bs_variance

! ----------------------------------
  real(r8) function bs_stdDevMean(self)
! ----------------------------------

    ! returns the standard deviation of the mean

    type(blockStat), intent(in)  :: self

    bs_stdDevMean = ss_stdDevMean(self%allBlockStat)

  end function bs_stdDevMean


! ------------------------------------
  integer(i8) function bs_dataCount(self)
! ------------------------------------

    type(blockStat), intent(in)  :: self

    bs_dataCount = ws_dataCount(self%totalWeightStat)

  end function bs_dataCount

! ---------------------------------
  logical function bs_isBlockEnd(self)
! ---------------------------------

    type(blockStat), intent(in)  :: self

    bs_isBlockEnd = self%blockEnd

  end function bs_isBlockEnd

  !--------------------------
  subroutine bs_debugPrint(self)
  !--------------------------

    type(blockStat), intent(in)  :: self

    print*,"blockStat:debugPrint:"
    print*,"blockLength=",self%blockLength
    print*,"blockCounter=",self%blockCounter
    print*,"discard=",self%discard
    print*,"blockEnd=",self%blockEnd
    print*,"allBlockStat:"
    call ss_debugPrint(self%allBlockStat)
    print*,"thisBlockStat:"
    call ws_debugPrint(self%thisBlockStat)
    print*,"totalWeightStat:"
    call ws_debugPrint(self%totalWeightStat)
    print*

  end subroutine bs_debugPrint



end module block_statistics

!=================================================


!=============================!
module int_histogram_statistics
!=============================!
  use kinds_m, only: i8
  use error_m
  implicit none

  type IntHistogram
     private
     integer, pointer  :: bin(:) => null()
     integer           :: xMin = 0
     integer           :: xMax = 0
     integer           :: h = 1         ! step size
     integer           :: n = 0         ! # of bins
     integer           :: averageSum = 0
     integer           :: varianceSum = 0
     integer(i8)         :: dataCounter = 0
     integer           :: minValue = 0
     integer           :: maxValue = 0
  end type IntHistogram

  interface assignment(=)
     module procedure assign_IntHistogram
  end interface

  interface operator(+)
     module procedure add_IntHistogram
  end interface

contains

! ---------------------------------------------------
  subroutine ih_createHistogram(self,xMin,xMax,nBins)
! ---------------------------------------------------

    type(IntHistogram), intent(inout) :: self
    integer, intent(in)    :: xMin
    integer, intent(inout) :: xMax   ! possibly adapted to generate
                                     ! bins with constant width
    integer, intent(inout) :: nBins  ! possibly adapted to generate
                                     ! bins with constant width

    call assert(.not.associated(self%bin),"IntHistogram: already initialized")
    call assert(xMax>xMin.and.nBins>1,"IntHistogram.init: wrong input")

    self%h = nint(dble(xMax-xMin)/nBins)
    nBins = nint(dble(xMax-xMin)/self%h)
    self%n = nBins
    xMax = xMin + nBins*self%h
    self%xMin = xMin
    self%xMax = xMax
    allocate(self%bin(0:self%n-1))
    self%bin = 0
  end subroutine ih_createHistogram


! -------------------------
  subroutine ih_reset(self)
! -------------------------
    ! resets all data, keeps bin structure

    type(IntHistogram), intent(inout) :: self

    self%bin = 0
    self%averageSum = 0
    self%varianceSum = 0
    self%dataCounter = 0
    self%minValue = 0
    self%maxValue = 0

  end subroutine ih_reset


! ---------------------------
  subroutine ih_destroy(self)
! ---------------------------

    type(IntHistogram), intent(inout) :: self

    call ih_reset(self)
    deallocate(self%bin)

  end subroutine ih_destroy


! --------------------------------
  subroutine ih_addData(self,data)
! --------------------------------

    type(IntHistogram), intent(inout) :: self
    integer, intent(in)               :: data
    integer kx

    call assert(associated(self%bin),"IntHistogram: object not allocated")
    kx = max(0,(data - self%xMin)/self%h)
    kx = min(kx,self%n-1)
    self%bin(kx) = self%bin(kx) + 1

    self%averageSum = self%averageSum + data
    self%varianceSum = self%varianceSum + data**2
    self%dataCounter = self%dataCounter + 1
    if (self%dataCounter==1) then
       self%minValue = data
       self%maxValue = data
    else
       self%minValue = min(data,self%minValue)
       self%maxValue = max(data,self%maxValue)
    end if

  end subroutine ih_addData


! ---------------------------------------
  subroutine assign_IntHistogram(self,stat)
! ---------------------------------------

    type(IntHistogram), intent(out) :: self
    type(IntHistogram), intent(in)  :: stat

    call assert(associated(stat%bin),"IntHistogram: assigning from uninitialized object")
    if (.not.associated(self%bin)) then
       allocate(self%bin(0:stat%n-1))
       self%bin = 0
    else
       call assert(size(self%bin)==size(stat%bin),"IntHistogram: assigning incompatable objects")
    endif
    self%bin = stat%bin
    self%n   = stat%n
    self%h   = stat%h
    self%xMin = stat%xMin
    self%xMax  = stat%xMax
    self%averageSum = stat%averageSum
    self%varianceSum = stat%varianceSum
    self%dataCounter = stat%dataCounter
    self%minValue = stat%minValue
    self%maxValue = stat%maxValue

  end subroutine assign_IntHistogram


! --------------------------------------------------
  type(IntHistogram) function add_IntHistogram(self,add)
! --------------------------------------------------

    type(IntHistogram), intent(in)  :: self
    type(IntHistogram), intent(in)  :: add

    call assert(associated(self%bin).and.associated(add%bin),"IntHistogram: adding uninitialized object")
    call assert(self%xMin==add%xMin.and.   &
                self%n==add%n .and. self%h==add%h,   &
                "IntHistogram: adding incompatable objects")
    if (.not.associated(add_IntHistogram%bin)) then
       allocate(add_IntHistogram%bin(0:self%n-1))
       add_IntHistogram%bin = 0
    else
       call assert(size(add_IntHistogram%bin)==size(self%bin),  &
                "IntHistogram: adding incompatable objects")
    endif
    add_IntHistogram%bin = self%bin + add%bin
    add_IntHistogram%n   = self%n
    add_IntHistogram%h   = self%h
    add_IntHistogram%xMin = self%xMin
    add_IntHistogram%xMax  = self%xMax
    add_IntHistogram%averageSum = self%averageSum + add%averageSum
    add_IntHistogram%varianceSum = self%varianceSum + add%varianceSum
    add_IntHistogram%dataCounter = self%dataCounter + add%dataCounter
    add_IntHistogram%minValue = min(add%minValue,self%minValue)
    add_IntHistogram%maxValue = max(add%maxValue,self%maxValue)

  end function add_IntHistogram


! access functions

! ------------------------------
  function ih_getHistogram(self)
! ------------------------------
    type(IntHistogram), intent(in)  :: self
    integer,pointer                 :: ih_getHistogram(:,:)
    integer i

    if (associated(ih_getHistogram)) deallocate(ih_getHistogram)
    allocate(ih_getHistogram(2,self%n))

    do i=0,self%n-1
       ih_getHistogram(1,i) = self%xMin + i*self%h
       ih_getHistogram(2,i) = self%bin(i)
    enddo

  end function ih_getHistogram


! -----------------------------
  real(r8) function ih_mean(self)
! -----------------------------

    type(IntHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"IntHistogram:mean")

    ih_mean = self%averageSum / self%dataCounter

  end function ih_mean


! --------------------------------
  real(r8) function ih_variance(self)
! --------------------------------

    type(IntHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"IntHistogram:variance")

    ih_variance = self%varianceSum/self%dataCounter - ih_mean(self)**2

  end function ih_variance

! ----------------------------------
  real(r8) function ih_stdDevMean(self)
! ----------------------------------

    ! returns the standard deviation of the mean

    type(IntHistogram), intent(in)  :: self

    call assert(self%dataCounter > 1,"IntHistogram:stdDevMean")

    ih_stdDevMean = sqrt(ih_variance(self) / (self%dataCounter - 1))

  end function ih_stdDevMean


! ------------------------------------
  integer(i8) function ih_dataCount(self)
! ------------------------------------

    type(IntHistogram), intent(in)  :: self

    ih_dataCount = self%dataCounter

  end function ih_dataCount

! ------------------------------
  real(r8) function ih_minValue(self)
! ------------------------------

    type(IntHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"IntHistogram:minValue")

    ih_minValue = self%minValue

  end function ih_minValue

! ------------------------------
  real(r8) function ih_maxValue(self)
! ------------------------------

    type(IntHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"IntHistogram:maxValue")

    ih_maxValue = self%maxValue

  end function ih_maxValue


  !--------------------------
  subroutine ih_debugPrint(self)
  !--------------------------

    type(IntHistogram), intent(in)  :: self

    print*,"debugPrint:averageSum=",self%averageSum
    print*,"debugPrint:varianceSum=",self%varianceSum
    print*,"debugPrint:dataCounter=",self%dataCounter

  end subroutine ih_debugPrint


! -------------------------------------
  subroutine ih_writeHistogram(self,iu)
! -------------------------------------

    type(IntHistogram), intent(in)  :: self
    integer, intent(in)             :: iu    ! file unit to write to
    integer i

    do i=0,self%n-1
       write(iu,*) self%xMin + i*self%h,self%bin(i)
    enddo

  end subroutine ih_writeHistogram


end module int_histogram_statistics

!=================================================


!================================!
module double_histogram_statistics
!================================!

    ! bin 0 contains all data smaller xMin and all data up to xMin + h
    ! the center of the bin is therefor xMin + 1/2*h
    ! the last (n-th) bin n-1 contains all data from xMax-h up to xMax (and beyond)
  use kinds_m, only: r8, i8
  use error_m
  implicit none

  type DoubleHistogram
     private
     real(r8), pointer  :: bin(:) => null()
     real(r8)           :: xMin = 0
     real(r8)           :: xMax = 0
     real(r8)           :: h = 1         ! step size
     integer          :: n = 0         ! # of bins
     real(r8)           :: averageSum = 0
     real(r8)           :: varianceSum = 0
     real(r8)           :: weightSum = 0
     integer(i8)        :: dataCounter = 0
     real(r8)           :: minValue = 0
     real(r8)           :: maxValue = 0
  end type DoubleHistogram

  interface assignment(=)
     module procedure assign_DoubleHistogram
  end interface

  interface operator(+)
     module procedure add_DoubleHistogram
  end interface

contains

! ---------------------------------------------------
  subroutine dh_createHistogram(self,xMin,xMax,nBins)
! ---------------------------------------------------

    type(DoubleHistogram), intent(inout) :: self
    real(r8), intent(in)    :: xMin
    real(r8), intent(in)    :: xMax
    integer, intent(in)   :: nBins

    call assert(.not.associated(self%bin),"DoubleHistogram: already initialized")
    call assert(xMax>xMin.and.nBins>1,"DoubleHistogram.init: wrong input")

    self%h = (xMax-xMin)/nBins
    self%n = nBins
    self%xMin = xMin
    self%xMax = xMax
    allocate(self%bin(0:self%n-1))
    self%bin = 0
  end subroutine dh_createHistogram


! -------------------------
  subroutine dh_reset(self)
! -------------------------
    ! resets all data, keeps bin structure

    type(DoubleHistogram), intent(inout) :: self

    self%bin = 0
    self%averageSum = 0
    self%varianceSum = 0
    self%weightSum = 0
    self%dataCounter = 0
    self%minValue = 0
    self%maxValue = 0

  end subroutine dh_reset


! ---------------------------
  subroutine dh_destroy(self)
! ---------------------------

    type(DoubleHistogram), intent(inout) :: self

    call dh_reset(self)
    deallocate(self%bin)

  end subroutine dh_destroy


! ---------------------------------------
  subroutine dh_addData(self,data,weight)
! ---------------------------------------

    type(DoubleHistogram), intent(inout) :: self
    real(r8), intent(in)               :: data
    real(r8), optional,intent(in)      :: weight
    integer kx
    real(r8) :: w

    if (.not.present(weight)) then
       w = 1
    else
       w = weight
    endif


    call assert(associated(self%bin),"DoubleHistogram: object not allocated")
    kx = max(0,int(dble(data - self%xMin)/self%h))
    kx = min(kx,self%n-1)
    self%bin(kx) = self%bin(kx) + w

    self%averageSum = self%averageSum + w*data
    self%varianceSum = self%varianceSum + w*data**2
    self%weightSum = self%weightSum + w
    self%dataCounter = self%dataCounter + 1
    if (self%dataCounter==1) then
       self%minValue = data
       self%maxValue = data
    else
       self%minValue = min(data,self%minValue)
       self%maxValue = max(data,self%maxValue)
    end if

  end subroutine dh_addData


! ---------------------------------------
  subroutine assign_DoubleHistogram(self,stat)
! ---------------------------------------

    type(DoubleHistogram), intent(out) :: self
    type(DoubleHistogram), intent(in)  :: stat

    call assert(associated(stat%bin),"DoubleHistogram: assigning from uninitialized object")
    if (.not.associated(self%bin)) then
       allocate(self%bin(0:stat%n-1))
       self%bin = 0
    else
       call assert(size(self%bin)==size(stat%bin),"DoubleHistogram: assigning incompatable objects")
    endif
    self%bin = stat%bin
    self%n   = stat%n
    self%h   = stat%h
    self%xMin = stat%xMin
    self%xMax  = stat%xMax
    self%averageSum = stat%averageSum
    self%varianceSum = stat%varianceSum
    self%dataCounter = stat%dataCounter
    self%minValue = stat%minValue
    self%maxValue = stat%maxValue

  end subroutine assign_DoubleHistogram


! --------------------------------------------------
  type(DoubleHistogram) function add_DoubleHistogram(self,add)
! --------------------------------------------------

    type(DoubleHistogram), intent(in)  :: self
    type(DoubleHistogram), intent(in)  :: add

    call assert(associated(self%bin).and.associated(add%bin),"DoubleHistogram: adding uninitialized object")
    call assert(self%xMin==add%xMin.and.   &
                self%n==add%n .and. self%h==add%h,   &
                "DoubleHistogram: adding incompatable objects")
    if (.not.associated(add_DoubleHistogram%bin)) then
       allocate(add_DoubleHistogram%bin(0:self%n-1))
    else
       call assert(size(add_DoubleHistogram%bin)==size(self%bin),  &
                "DoubleHistogram: adding incompatable objects")
    endif
    add_DoubleHistogram%bin = self%bin + add%bin
    add_DoubleHistogram%n   = self%n
    add_DoubleHistogram%h   = self%h
    add_DoubleHistogram%xMin = self%xMin
    add_DoubleHistogram%xMax  = self%xMax
    add_DoubleHistogram%averageSum = self%averageSum + add%averageSum
    add_DoubleHistogram%varianceSum = self%varianceSum + add%varianceSum
    add_DoubleHistogram%dataCounter = self%dataCounter + add%dataCounter
    add_DoubleHistogram%minValue = min(add%minValue,self%minValue)
    add_DoubleHistogram%maxValue = max(add%maxValue,self%maxValue)

  end function add_DoubleHistogram


! access functions

! ------------------------------
  function dh_getHistogram(self)
! ------------------------------
    ! returns in getHistogram(1,:) start value of bin and
    ! in getHistogram(2,:) size of bin
    type(DoubleHistogram), intent(in)  :: self
    real(r8),pointer                     :: dh_getHistogram(:,:)
    integer i

    if (associated(dh_getHistogram)) deallocate(dh_getHistogram)
    allocate(dh_getHistogram(2,self%n))

    do i=0,self%n-1
       dh_getHistogram(1,i) = self%xMin + i*self%h
       dh_getHistogram(2,i) = self%bin(i)
    enddo

  end function dh_getHistogram


! ----------------------------------
  function dh_getMeanHistogram(self)
! ----------------------------------
    ! returns in getHistogram(1,:) mean x value of bin and
    ! in getHistogram(2,:) size of bin
    type(DoubleHistogram), intent(in)  :: self
    real(r8),pointer                 :: dh_getMeanHistogram(:,:)
    integer i

    if (associated(dh_getMeanHistogram)) deallocate(dh_getMeanHistogram)
    allocate(dh_getMeanHistogram(2,self%n))

    do i=0,self%n-1
       dh_getMeanHistogram(1,i) = self%xMin + (i+0.5d0)*self%h
       dh_getMeanHistogram(2,i) = self%bin(i)
    enddo

  end function dh_getMeanHistogram


! ----------------------------
  real(r8) function dh_mean(self)
! ----------------------------

    type(DoubleHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"DoubleHistogram:mean")

    dh_mean = self%averageSum / self%weightSum

  end function dh_mean


! --------------------------------
  real(r8) function dh_variance(self)
! --------------------------------

    type(DoubleHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"DoubleHistogram:variance")

    dh_variance = self%varianceSum/self%weightSum - dh_mean(self)**2

  end function dh_variance

! ----------------------------------
  real(r8) function dh_stdDevMean(self)
! ----------------------------------

    ! returns the standard deviation of the mean

    type(DoubleHistogram), intent(in)  :: self

    call assert(self%dataCounter > 1,"DoubleHistogram:stdDevMean")

    dh_stdDevMean = sqrt(dh_variance(self) / (self%dataCounter - 1))

  end function dh_stdDevMean


! ------------------------------------
  integer(i8) function dh_dataCount(self)
! ------------------------------------

    type(DoubleHistogram), intent(in)  :: self

    dh_dataCount = self%dataCounter

  end function dh_dataCount

! ------------------------------------
  real(r8) function dh_totalWeight(self)
! ------------------------------------

    type(DoubleHistogram), intent(in)  :: self

    dh_totalWeight = self%weightSum

  end function dh_totalWeight


! ---------------------------------
  real(r8) function dh_minValue(self)
! ---------------------------------

    type(DoubleHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"DoubleHistogram:minValue")

    dh_minValue = self%minValue

  end function dh_minValue

! ---------------------------------
  real(r8) function dh_maxValue(self)
! ---------------------------------

    type(DoubleHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"DoubleHistogram:maxValue")

    dh_maxValue = self%maxValue

  end function dh_maxValue


  !-----------------------------
  subroutine dh_debugPrint(self)
  !-----------------------------

    type(DoubleHistogram), intent(in)  :: self

    print*,"debugPrint:averageSum=",self%averageSum
    print*,"debugPrint:varianceSum=",self%varianceSum
    print*,"debugPrint:weightSum=",self%weightSum
    print*,"debugPrint:dataCounter=",self%dataCounter

  end subroutine dh_debugPrint


! -------------------------------------
  subroutine dh_writeHistogram(self,iu)
! -------------------------------------
    ! writing Histogram as x,y pairs where x is the left value
    ! of the bin
    type(DoubleHistogram), intent(in)  :: self
    integer, intent(in)             :: iu    ! file unit to write to
    integer i
    real(r8) mean,denom,x,normaldist(0:self%n-1),factor

    mean = dh_mean(self)
    denom = sqrt(2.d0*dh_variance(self))
    do i=0,self%n-1
       x = self%xMin + i*self%h
       normaldist(i) = 0.5d0*(erf((x+self%h - mean)/denom) - erf((x - mean)/denom))
    enddo
    factor = self%weightSum / sum(normaldist)

    do i=0,self%n-1
       x = self%xMin + i*self%h
       write(iu,*) x, self%bin(i), factor*normaldist(i)
    enddo

  end subroutine dh_writeHistogram


! -----------------------------------------
  subroutine dh_writeMeanHistogram(self,iu)
! -----------------------------------------
    ! writing Histogram as x,y pairs where x is the mean value
    ! of the bin
    type(DoubleHistogram), intent(in)  :: self
    integer, intent(in)             :: iu    ! file unit to write to
    integer i
    real(r8) mean,denom,x,normaldist(0:self%n-1),factor

    mean = dh_mean(self)
    denom = sqrt(2.d0*dh_variance(self))
    do i=0,self%n-1
       x = self%xMin + i*self%h
       normaldist(i) = 0.5d0*(erf((x+self%h - mean)/denom) - erf((x - mean)/denom))
    enddo
    factor = self%weightSum / sum(normaldist)

    do i=0,self%n-1
       x = self%xMin + i*self%h
       write(iu,*) x + 0.5d0*self%h, self%bin(i), factor*normaldist(i)
    enddo

  end subroutine dh_writeMeanHistogram

end module double_histogram_statistics



!=================================================


!=================================!
module indexed_histogram_statistics
!=================================!

    ! bin 0 contains all data smaller xMin and all data up to xMin + h
    ! the center of the bin is therefor xMin + 1/2*h
    ! the last (n-th) bin n-1 contains all data from xMax-h up to xMax (and beyond)
  use kinds_m, only: r8, i8
  use error_m
  use intList_m
  implicit none

  type IndexedHistogram
     private
     real(r8), pointer    :: bin(:) => null()
     real(r8)             :: xMin = 0
     real(r8)             :: xMax = 0
     real(r8)             :: h = 1         ! step size
     integer            :: n = 0         ! # of bins
     real(r8)             :: averageSum = 0
     real(r8)             :: varianceSum = 0
     real(r8)             :: weightSum = 0
     integer            :: dataCounter = 0
     integer, pointer   :: binCounter(:) => null()
     integer, pointer   :: indexList(:) => null()
     real(r8)             :: minValue = 0
     real(r8)             :: maxValue = 0
     integer            :: maxData = 0
  end type IndexedHistogram

  interface assignment(=)
     module procedure assign_IndexedHistogram
  end interface

  interface operator(+)
     module procedure add_IndexedHistogram
  end interface

  private
  public :: IndexedHistogram,idh_createHistogram,idh_reset,idh_destroy,idh_addData, assign_IndexedHistogram, &
            add_indexedHistogram, idh_getHistogram, idh_getMeanHistogram,idh_mean,idh_variance, idh_stdDevMean, &
            idh_dataCount,idh_minValue,idh_maxValue,idh_debugPrint,idh_writeHistogram,idh_writeMeanHistogram, &
            idh_getOutliers

contains

! ------------------------------------------------------------
  subroutine idh_createHistogram(self,xMin,xMax,nBins,maxData)
! ------------------------------------------------------------

    type(IndexedHistogram), intent(inout) :: self
    real(r8), intent(in)    :: xMin
    real(r8), intent(in)    :: xMax
    integer, intent(in)   :: nBins
    integer, intent(in)   :: maxData   ! expected maximum data count
    integer alstat

    write(*,*) "DBG:idx_createHistogram:",associated(self%bin)
    call assert(.not.associated(self%bin),"IndexedHistogram: already initialized")
    call assert(xMax>xMin.and.nBins>1 .and. maxData>1,"IndexedHistogram.init: illegal input")

    self%h = (xMax-xMin)/nBins
    self%n = nBins
    self%xMin = xMin
    self%xMax = xMax
    self%maxData = maxData
    allocate(self%bin(0:self%n-1),self%binCounter(0:self%n-1),self%indexList(maxData),stat=alstat)
    self%bin = 0; self%binCounter = 0; self%indexList = -1
    call assert(alstat==0,'indexedHistogram: allocation failed')
  end subroutine idh_createHistogram


! ------------------------------------
  subroutine idh_reset(self,xMin,xMax)
! ------------------------------------
    ! resets all data, keeps bin structure

    type(IndexedHistogram), intent(inout) :: self
    real(r8), optional, intent(in)  :: xMin
    real(r8), optional, intent(in)  :: xMax

    if (present(xMin) .and. present(xMax)) then
       self%h = (xMax-xMin)/self%n
       self%xMin = xMin
       self%xMax = xMax
    endif
    self%bin = 0
    self%averageSum = 0
    self%varianceSum = 0
    self%weightSum = 0
    self%dataCounter = 0
    self%binCounter = 0
    self%minValue = 0
    self%maxValue = 0
    self%indexList = -1
  end subroutine idh_reset


! ---------------------------
  subroutine idh_destroy(self)
! ---------------------------

    type(IndexedHistogram), intent(inout) :: self

    call idh_reset(self)
    deallocate(self%bin,self%binCounter,self%indexList)

  end subroutine idh_destroy


! ----------------------------------------
  subroutine idh_addData(self,data,weight)
! ----------------------------------------

    type(IndexedHistogram), intent(inout) :: self
    real(r8), intent(in)               :: data
    real(r8), optional,intent(in)      :: weight
    integer kx
    real(r8) :: w

    if (.not.present(weight)) then
       w = 1
    else
       w = weight
    endif

    call assert(associated(self%bin),"IndexedHistogram: object not allocated")
    call assert(self%dataCounter < self%maxData,"indexedHistogram: too many data")
    kx = max(0,int(dble(data - self%xMin)/self%h))
    kx = min(kx,self%n-1)
    self%bin(kx) = self%bin(kx) + w
    self%binCounter(kx) = self%binCounter(kx) + 1
    self%averageSum = self%averageSum + w*data
    self%varianceSum = self%varianceSum + w*data**2
    self%weightSum = self%weightSum + w
    self%dataCounter = self%dataCounter + 1
    self%indexList(self%dataCounter) = kx
    if (self%dataCounter==1) then
       self%minValue = data
       self%maxValue = data
    else
       self%minValue = min(data,self%minValue)
       self%maxValue = max(data,self%maxValue)
    end if

  end subroutine idh_addData


! ---------------------------------------------
  subroutine assign_IndexedHistogram(self,stat)
! ---------------------------------------------

    type(IndexedHistogram), intent(out) :: self
    type(IndexedHistogram), intent(in)  :: stat
    integer alstat

    call assert(associated(stat%bin),"IndexedHistogram: assigning from uninitialized object")
    if (.not.associated(self%bin)) then
       allocate(self%bin(0:stat%n-1),self%binCounter(0:stat%n-1),self%indexList(0:stat%n-1),stat=alstat)
       call assert(alstat==0,'assign_IndexedHistogram: allocation failed')
    else
       call assert(size(self%bin)==size(stat%bin),"IndexedHistogram: assigning incompatible objects")
    endif
    self%bin = stat%bin
    self%n   = stat%n
    self%h   = stat%h
    self%xMin = stat%xMin
    self%xMax  = stat%xMax
    self%indexList = stat%indexList
    self%averageSum = stat%averageSum
    self%varianceSum = stat%varianceSum
    self%weightSum = stat%weightSum
    self%dataCounter = stat%dataCounter
    self%binCounter = stat%binCounter
    self%minValue = stat%minValue
    self%maxValue = stat%maxValue

  end subroutine assign_IndexedHistogram


! --------------------------------------------------------------
  type(IndexedHistogram) function add_IndexedHistogram(self,add)
! --------------------------------------------------------------

    type(IndexedHistogram), intent(in)  :: self
    type(IndexedHistogram), intent(in)  :: add
    integer alstat,n1,n2

    call assert(associated(self%bin).and.associated(add%bin),"IndexedHistogram: adding uninitialized object")
    call assert(self%xMin==add%xMin.and.   &
                self%n==add%n .and. self%h==add%h,   &
                "IndexedHistogram: adding incompatable objects")
    if (.not.associated(add_IndexedHistogram%bin)) then
       allocate(add_IndexedHistogram%bin(0:self%n-1), add_IndexedHistogram%indexList(0:self%n-1),stat=alstat)
       call assert(alstat==0,'add_IndexedHistogram: allocation failed')
    else
       call assert(size(add_IndexedHistogram%bin)==size(self%bin),  &
                "IndexedHistogram: adding incompatable objects")
    endif
    add_IndexedHistogram%bin = self%bin + add%bin
    add_IndexedHistogram%n   = self%n
    add_IndexedHistogram%h   = self%h
    add_IndexedHistogram%xMin = self%xMin
    add_IndexedHistogram%xMax  = self%xMax
    n1 = self%dataCounter; n2 = add%dataCounter
    add_IndexedHistogram%indexList(1:n1) = self%indexList(1:n1)
    add_IndexedHistogram%indexList(n1+1:n1+n2) = add%indexList(1:n2)
    add_IndexedHistogram%averageSum = self%averageSum + add%averageSum
    add_IndexedHistogram%varianceSum = self%varianceSum + add%varianceSum
    add_IndexedHistogram%weightSum = self%weightSum + add%weightSum
    add_IndexedHistogram%dataCounter = self%dataCounter + add%dataCounter
    add_IndexedHistogram%binCounter = self%binCounter + add%binCounter
    add_IndexedHistogram%minValue = min(add%minValue,self%minValue)
    add_IndexedHistogram%maxValue = max(add%maxValue,self%maxValue)

  end function add_IndexedHistogram


! access functions

! ------------------------------
  function idh_getHistogram(self)
! ------------------------------
    ! returns in getHistogram(1,:) start value of bin and
    ! in getHistogram(2,:) size of bin
    type(IndexedHistogram), intent(in)  :: self
    real(r8),pointer                     :: idh_getHistogram(:,:)
    integer i

    if (associated(idh_getHistogram)) deallocate(idh_getHistogram)
    allocate(idh_getHistogram(2,self%n))

    do i=0,self%n-1
       idh_getHistogram(1,i) = self%xMin + i*self%h
       idh_getHistogram(2,i) = self%bin(i)
    enddo

  end function idh_getHistogram


! ----------------------------------
  function idh_getMeanHistogram(self)
! ----------------------------------
    ! returns in getHistogram(1,:) mean x value of bin and
    ! in getHistogram(2,:) size of bin
    type(IndexedHistogram), intent(in)  :: self
    real(r8),pointer                 :: idh_getMeanHistogram(:,:)
    integer i

    if (associated(idh_getMeanHistogram)) deallocate(idh_getMeanHistogram)
    allocate(idh_getMeanHistogram(2,self%n))

    do i=0,self%n-1
       idh_getMeanHistogram(1,i) = self%xMin + (i+0.5d0)*self%h
       idh_getMeanHistogram(2,i) = self%bin(i)
    enddo

  end function idh_getMeanHistogram

! ------------------------------
  real(r8) function idh_mean(self)
! -----------------------------

    type(IndexedHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"IndexedHistogram:mean")

    idh_mean = self%averageSum / self%weightSum

  end function idh_mean


! --------------------------------
  real(r8) function idh_variance(self)
! --------------------------------

    type(IndexedHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"IndexedHistogram:variance")

    idh_variance = self%varianceSum/self%weightSum - idh_mean(self)**2

  end function idh_variance

! ----------------------------------
  real(r8) function idh_stdDevMean(self)
! ----------------------------------

    ! returns the standard deviation of the mean

    type(IndexedHistogram), intent(in)  :: self

    call assert(self%dataCounter > 1,"IndexedHistogram:stdDevMean")

    idh_stdDevMean = sqrt(idh_variance(self) / (self%dataCounter - 1))

  end function idh_stdDevMean

! ------------------------------------
  integer(i8) function idh_dataCount(self)
! ------------------------------------

    type(IndexedHistogram), intent(in)  :: self

    idh_dataCount = self%dataCounter

  end function idh_dataCount

! ---------------------------------
  real(r8) function idh_minValue(self)
! ---------------------------------

    type(IndexedHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"IndexedHistogram:minValue")

    idh_minValue = self%minValue

  end function idh_minValue

! ---------------------------------
  real(r8) function idh_maxValue(self)
! ---------------------------------

    type(IndexedHistogram), intent(in)  :: self

    call assert(self%dataCounter > 0,"IndexedHistogram:maxValue")

    idh_maxValue = self%maxValue

  end function idh_maxValue


  !-----------------------------
  subroutine idh_debugPrint(self)
  !-----------------------------

    type(IndexedHistogram), intent(in)  :: self

    print*,"debugPrint:averageSum=",self%averageSum
    print*,"debugPrint:varianceSum=",self%varianceSum
    print*,"debugPrint:dataCounter=",self%dataCounter

  end subroutine idh_debugPrint


! -------------------------------------
  subroutine idh_writeHistogram(self,iu)
! -------------------------------------
    ! writing Histogram as x,y pairs where x is the left value
    ! of the bin
    type(IndexedHistogram), intent(in)  :: self
    integer, intent(in)             :: iu    ! file unit to write to
    integer i
    real(r8) mean,denom,x,normaldist(0:self%n-1),factor

    mean = idh_mean(self)
    denom = sqrt(2.d0*idh_variance(self))
    do i=0,self%n-1
       x = self%xMin + i*self%h
       normaldist(i) = 0.5d0*(erf((x+self%h - mean)/denom) - erf((x - mean)/denom))
    enddo
    factor = self%weightSum / sum(normaldist)

    do i=0,self%n-1
       x = self%xMin + i*self%h
       write(iu,*) x, self%bin(i), factor*normaldist(i)
    enddo
  end subroutine idh_writeHistogram


! -----------------------------------------
  subroutine idh_writeMeanHistogram(self,iu)
! -----------------------------------------
    ! writing Histogram as x,y pairs where x is the mean value
    ! of the bin
    type(IndexedHistogram), intent(in)  :: self
    integer, intent(in)             :: iu    ! file unit to write to
    integer i
    real(r8) mean,denom,x,normaldist(0:self%n-1),factor

    mean = idh_mean(self)
    denom = sqrt(2.d0*idh_variance(self))
    do i=0,self%n-1
       x = self%xMin + i*self%h
       normaldist(i) = 0.5d0*(erf((x+self%h - mean)/denom) - erf((x - mean)/denom))
    enddo
    factor = self%weightSum / sum(normaldist)

    do i=0,self%n-1
       x = self%xMin + i*self%h
       write(iu,*) x + 0.5d0*self%h, self%bin(i), factor*normaldist(i)
    enddo
  end subroutine idh_writeMeanHistogram

  !--------------------------------------------------------------------!
  subroutine idh_getOutliers(self,leftList,leftProb,rightList,rightProb)
  !--------------------------------------------------------------------!
     ! this returns the indices of the leftmost and rightmost non-empty bins
     ! with corresponding probable number (for a normal distribution) of entries
     ! note: the lists are pointers that are reallocated
    type(IndexedHistogram), intent(in)  :: self
    type(IntList), intent(inout)        :: leftList, rightList
    real(r8), intent(out)                 :: leftProb,rightProb
    real(r8) norm,mean,denom,x,normaldist
    integer kl,kr,j

    call assert(self%dataCounter > 2*self%n,'idh_getOutliers: too few data for analysis')
    if (isCreatedList(leftList)) then
       call clear(leftList)
    else
       call createList(leftList)
    endif
    if (isCreatedList(rightList)) then
       call clear(rightList)
    else
       call createList(rightList)
    endif

    ! get leftmost bin
    do kl=0,self%n/2
       if (self%binCounter(kl) > 0) exit
    enddo
    call assert(kl<=self%n/2,'idh_getOutliers: empty bins')
    do j=1,self%dataCounter
       if (self%indexList(j)==kl) then
          call pushBack(leftList,j)
       endif
    enddo
    call assert(getSize(leftList)==self%binCounter(kl),'idh_getOutliers: corrupted indexList')
    ! get rightmost bin
    do kr=self%n-1,self%n/2,-1
       if (self%binCounter(kr) > 0) exit
    enddo
    call assert(kr>=self%n/2,'idh_getOutliers: empty bins')
    do j=1,self%dataCounter
       if (self%indexList(j)==kr) then
          call pushBack(rightList,j)
       endif
    enddo
    call assert(getSize(rightList)==self%binCounter(kr),'idh_getOutliers: corrupted indexList')

    ! calculate corresponding expected counts (from normal distribution)
    norm = self%weightSum
    mean = idh_mean(self)
    denom = sqrt(2.d0*idh_variance(self))

    x = self%xMin + kl*self%h
    normaldist = 0.5d0*(erf((x+self%h - mean)/denom) - erf((x - mean)/denom))
    leftProb = normaldist*norm
    x = self%xMin + kr*self%h
    normaldist = 0.5d0*(erf((x+self%h - mean)/denom) - erf((x - mean)/denom))
    rightProb = normaldist*norm
  end subroutine idh_getOutliers

end module indexed_histogram_statistics



!==========================================================

! Generic module to allow use of same names for different
!
!================
module statistics
!================

  use simple_statistics
  use weight_statistics
  use int_histogram_statistics
  use double_histogram_statistics

  implicit none

  interface addData
     module procedure ss_addData,ws_addData,ih_addData,dh_addData
  end interface addData

  interface reset
     module procedure ss_reset,ws_reset,ih_reset,dh_reset
  end interface reset

  interface mean
     module procedure ss_mean,ws_mean,ih_mean,dh_mean
  end interface mean

  interface meanAllNodes
     module procedure ss_meanAllNodes,ws_meanAllNodes
  end interface meanAllNodes

  interface variance
     module procedure ss_variance,ws_variance,ih_variance,dh_variance
  end interface variance

  interface varianceAllNodes
     module procedure ss_varianceAllNodes,ws_varianceAllNodes
  end interface varianceAllNodes

  interface stdDevMean
     module procedure ss_stdDevMean, ws_stdDevMean,ih_stdDevMean,dh_stdDevMean
  end interface stdDevMean

  interface stdDevMeanAllNodes
     module procedure ss_stdDevMeanAllNodes, ws_stdDevMeanAllNodes
  end interface

  interface dataCount
     module procedure ss_dataCount,ws_dataCount,ih_dataCount,dh_dataCount
  end interface dataCount

  interface dataCountAllNodes
     module procedure ss_dataCountAllNodes,ws_dataCountAllNodes
  end interface dataCountAllNodes

  interface totalWeight
     module procedure ws_totalWeight
  end interface totalWeight

  interface totalWeightAllNodes
     module procedure ws_totalWeightAllNodes
  end interface totalWeightAllNodes

  interface minValue
     module procedure ss_minValue,ws_minValue,ih_minValue,dh_minValue
  end interface minValue

  interface maxValue
     module procedure ss_maxValue,ws_maxValue,ih_maxValue,dh_maxValue
  end interface maxValue

  interface debugPrint
     module procedure ss_debugPrint,ws_debugPrint
  end interface debugPrint

  interface createHistogram
     module procedure ih_createHistogram,dh_createHistogram
  end interface

  interface getHistogram
     module procedure ih_getHistogram,dh_getHistogram
  end interface

  interface getMeanHistogram
     module procedure dh_getMeanHistogram
  end interface

  interface writeHistogram
     module procedure ih_writeHistogram,dh_writeHistogram
  end interface

  interface writeMeanHistogram
     module procedure dh_writeMeanHistogram
  end interface

end module statistics

