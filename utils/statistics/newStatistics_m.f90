! Copyright (C) 2013 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module newStatistics_m

   ! new version of statistics objects
   ! designed for parallel efficiency
   ! mean returns mean for all nodes; var, stddev similarly
   ! mean, var, stddev, and weight return correct value ONLY FOR MASTER
   ! count returns value on ALL NODES
   ! NOTE: weighted version has ONE weight for each vector/matrix
   ! use array of stat/bstat for individual weights (or add corresponding types)
   !

   ! TODO: blocking and/or Flyvbjerg/Petersen method (for correct stddev)
   ! TODO: get current block mean (for logging, trends!)

   use kinds_m, only: r8, i8
   use error_m
   implicit none
   private

   integer, parameter :: MASTER=0

   type, public :: stat
      private
      real(r8)    :: averageSum = 0
      real(r8)    :: varianceSum = 0
      real(r8)    :: weightSum = 0
      integer(i8) :: dataCounter = 0
   contains
      procedure :: create => stat_create
      procedure :: destroy => stat_destroy
      procedure :: reset => stat_reset
      procedure :: add => stat_addData
      procedure :: mean => stat_mean
      procedure :: var => stat_var
      procedure :: stddev => stat_stddev
      procedure :: weight => stat_weight
      procedure :: count => stat_count
      procedure :: localmean => stat_localmean
      procedure :: localvar => stat_localvar
      procedure :: localstddev => stat_localstddev
      procedure :: localweight => stat_localweight
      procedure :: localcount => stat_localcount
   end type stat

   type, public :: blockstat
      private
      type(stat) :: allBlockStat
      type(stat) :: thisBlockStat
      type(stat) :: totalStat
      real(r8)     :: lastBlockMeanValue = 0
      integer    :: blockLen = 0
   contains
      procedure :: create => bstat_create
      procedure :: destroy => bstat_destroy
      procedure :: reset => bstat_reset
      procedure :: blocklength => bstat_blocklen
      procedure :: add => bstat_addData
      procedure :: mean => bstat_mean
      procedure :: var => bstat_var
      procedure :: stddev => bstat_stddev
      procedure :: weight => bstat_weight
      procedure :: count => bstat_count
      procedure :: lastblockmean => bstat_lastblockmean
      procedure :: blockcount => bstat_blockcount
      procedure :: localmean => bstat_localmean
      procedure :: localvar => bstat_localvar
      procedure :: localstddev => bstat_localstddev
      procedure :: localweight => bstat_localweight
      procedure :: localcount => bstat_localcount
      procedure :: locallastblockmean => bstat_locallastblockmean
      procedure :: localblockcount => bstat_localblockcount
   end type blockstat

   type, public :: vectorstat
      private
      real(r8),pointer    :: averageSum(:) => null()
      real(r8),pointer    :: varianceSum(:) => null()
      real(r8)            :: weightSum = 0
      integer(i8)         :: dataCounter = 0
      integer           :: vsize = 0
   contains
      procedure :: create => vstat_create
      procedure :: iscreated => vstat_iscreated
      procedure :: destroy => vstat_destroy
      procedure :: reset => vstat_reset
      procedure :: add => vstat_addData
      procedure :: mean => vstat_mean
      procedure :: var => vstat_var
      procedure :: stddev => vstat_stddev
      procedure :: weight => vstat_weight
      procedure :: count => vstat_count
      procedure :: localmean => vstat_localmean
      procedure :: localvar => vstat_localvar
      procedure :: localstddev => vstat_localstddev
      procedure :: localweight => vstat_localweight
      procedure :: localcount => vstat_localcount
   end type vectorstat

   type, public :: blockvectorstat
      private
      type(vectorstat) :: allBlockStat
      type(vectorstat) :: thisBlockStat
      type(vectorstat) :: totalStat
      real(r8),pointer   :: lastBlockMeanValue(:) => null()
      integer    :: blockLen = 0
   contains
      procedure :: create => bvstat_create
      procedure :: iscreated => bvstat_iscreated
      procedure :: destroy => bvstat_destroy
      procedure :: reset => bvstat_reset
      procedure :: blocklength => bvstat_blocklen
      procedure :: add => bvstat_addData
      procedure :: mean => bvstat_mean
      procedure :: var => bvstat_var
      procedure :: stddev => bvstat_stddev
      procedure :: weight => bvstat_weight
      procedure :: count => bvstat_count
      procedure :: lastblockmean => bvstat_lastblockmean
      procedure :: blockcount => bvstat_blockcount
      procedure :: localmean => bvstat_localmean
      procedure :: localvar => bvstat_localvar
      procedure :: localstddev => bvstat_localstddev
      procedure :: localweight => bvstat_localweight
      procedure :: localcount => bvstat_localcount
      procedure :: locallastblockmean => bvstat_locallastblockmean
      procedure :: localblockcount => bvstat_localblockcount
   end type blockvectorstat

   type, public :: matrixstat
      private
      real(r8),pointer    :: averageSum(:,:) => null()
      real(r8),pointer    :: varianceSum(:,:) => null()
      real(r8)            :: weightSum = 0
      integer(i8)         :: dataCounter = 0
      integer           :: vsize = 0
      integer           :: m=0
      integer           :: n=0
   contains
      procedure :: create => mstat_create
      procedure :: iscreated => mstat_iscreated
      procedure :: destroy => mstat_destroy
      procedure :: reset => mstat_reset
      procedure :: add => mstat_addData
      procedure :: mean => mstat_mean
      procedure :: var => mstat_var
      procedure :: stddev => mstat_stddev
      procedure :: weight => mstat_weight
      procedure :: count => mstat_count
      procedure :: localmean => mstat_localmean
      procedure :: localvar => mstat_localvar
      procedure :: localstddev => mstat_localstddev
      procedure :: localweight => mstat_localweight
      procedure :: localcount => mstat_localcount
   end type matrixstat

   type, public :: blockmatrixstat
      private
      type(matrixstat) :: allBlockStat
      type(matrixstat) :: thisBlockStat
      type(matrixstat) :: totalStat
      real(r8),pointer   :: lastBlockMeanValue(:,:) => null()
      integer    :: blockLen = 0
   contains
      procedure :: create => bmstat_create
      procedure :: iscreated => bmstat_iscreated
      procedure :: destroy => bmstat_destroy
      procedure :: reset => bmstat_reset
      procedure :: blocklength => bmstat_blocklen
      procedure :: add => bmstat_addData
      procedure :: mean => bmstat_mean
      procedure :: var => bmstat_var
      procedure :: stddev => bmstat_stddev
      procedure :: weight => bmstat_weight
      procedure :: count => bmstat_count
      procedure :: lastblockmean => bmstat_lastblockmean
      procedure :: blockcount => bmstat_blockcount
      procedure :: localmean => bmstat_localmean
      procedure :: localvar => bmstat_localvar
      procedure :: localstddev => bmstat_localstddev
      procedure :: localweight => bmstat_localweight
      procedure :: localcount => bmstat_localcount
      procedure :: locallastblockmean => bmstat_locallastblockmean
      procedure :: localblockcount => bmstat_localblockcount
   end type blockmatrixstat

contains

   subroutine stat_create(self)
      class(stat), intent(inout) :: self
      call stat_reset(self)
   end subroutine

   subroutine stat_destroy(self)
      class(stat), intent(inout) :: self
      self%averageSum = 0
      self%varianceSum = 0
      self%weightSum = 0
      self%dataCounter = 0
   end subroutine

   subroutine stat_reset(self)
      class(stat), intent(inout) :: self
      self%averageSum = 0
      self%varianceSum = 0
      self%weightSum = 0
      self%dataCounter = 0
   end subroutine

   subroutine stat_addData(self,data,wgt)
      class(stat), intent(inout) :: self
      real(r8), intent(in)         :: data
      real(r8), intent(in), optional :: wgt
      real(r8) w
      if (.not.present(wgt)) then
         w = 1
      else
         w = wgt
      endif

      self%averageSum = self%averageSum + w*data
      self%varianceSum = self%varianceSum + w*data**2
      self%weightSum = self%weightSum + w
      self%dataCounter = self%dataCounter + 1
   end subroutine stat_addData


! access functions


   real(r8) function stat_localmean(self)
      class(stat), intent(in)  :: self
      call assert(self%weightSum > 0,"stat:mean no data or no weight")
      stat_localmean = self%averageSum / self%weightSum
   end function stat_localmean


   real(r8) function stat_mean(self)
#ifdef MPI
      use MPI_F08
#endif
      class(stat), intent(in)  :: self
      real(r8) sendbuf(2),recvbuf(2)
      integer taskid
#ifdef MPI
      integer ierr
#endif
      sendbuf(1) = self%averageSum; sendbuf(2) = self%weightSum
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,2,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0) then
         call assert(recvbuf(2) > 0,"stat:mean no data or no weight")
         stat_mean = recvbuf(1) / recvbuf(2)
      else
         stat_mean = 0.d0
      end if
   end function stat_mean


   real(r8) function stat_localvar(self)
      class(stat), intent(in)  :: self
      call assert(self%weightSum > 0,"stat:var no data or no weight")
      stat_localvar = max(0.d0,self%varianceSum/self%weightSum - stat_localmean(self)**2)
   end function stat_localvar


   real(r8) function stat_var(self)
#ifdef MPI
      use MPI_F08
#endif
      class(stat), intent(in)  :: self
      real(r8) sendbuf(3),recvbuf(3)
      integer taskid
#ifdef MPI
      integer ierr
#endif
      sendbuf(1) = self%averageSum; sendbuf(2)=self%varianceSum; sendbuf(3) = self%weightSum
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,3,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0) then
         call assert(recvbuf(3) > 0,"stat:var no data or no weight")
         stat_var = max(0.d0,recvbuf(2)/recvbuf(3)-(recvbuf(1)/recvbuf(3))**2)
      else
         stat_var = 0.d0
      end if
   end function stat_var


   real(r8) function stat_localstddev(self)
      class(stat), intent(in)  :: self
      if (self%dataCounter > 1) then
         stat_localstddev = sqrt(stat_localvar(self) / (self%dataCounter - 1))
      else
         stat_localstddev = 0.d0
      end if
   end function stat_localstddev


   real(r8) function stat_stddev(self)
#ifdef MPI
      use MPI_F08
#endif
      class(stat), intent(in)  :: self
      real(r8) sendbuf(4),recvbuf(4)
      integer taskid
#ifdef MPI
      integer ierr
#endif
      sendbuf(1) = self%averageSum; sendbuf(2)=self%varianceSum
      sendbuf(3) = self%weightSum; sendbuf(4)=self%dataCounter
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,4,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0 .and. recvbuf(4) > 1.d0) then
         stat_stddev = sqrt( max(0.d0,(recvbuf(2)/recvbuf(3)-(recvbuf(1)/recvbuf(3))**2)/(recvbuf(4)-1.d0)) )
      else
         stat_stddev = 0.d0
      end if
   end function stat_stddev


   integer(i8) function stat_localcount(self)
      class(stat), intent(in)  :: self
      stat_localcount = self%dataCounter
   end function stat_localcount


   integer(i8) function stat_count(self)
#ifdef MPI
      use MPI_F08
#endif
      class(stat), intent(in)  :: self
      integer(i8) send,recv
#ifdef MPI
      integer ierr
#endif
      send = self%dataCounter
      recv = send  ! default for serial run
#ifdef MPI
      call MPI_REDUCE(send,recv,1,MPI_INTEGER8,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(recv,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,ierr)
#endif
      stat_count = recv
   end function stat_count


   real(r8) function stat_localweight(self)
      class(stat), intent(in)  :: self
      stat_localweight = self%weightSum
   end function stat_localweight


   real(r8) function stat_weight(self)
#ifdef MPI
      use MPI_F08
#endif
      class(stat), intent(in)  :: self
      real(r8) send,recv
#ifdef MPI
      integer ierr
#endif
      send = self%dataCounter
      recv = send  ! default for serial run
#ifdef MPI
      call MPI_REDUCE(send,recv,1,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
#endif
      stat_weight = recv
   end function stat_weight


 !!!!!!!!! block stat methods

   subroutine bstat_create(self,blockLen)
      class(blockstat), intent(inout) :: self
      integer, intent(in) :: blockLen
      call bstat_reset(self)
      self%blockLen = blockLen
   end subroutine

   subroutine bstat_destroy(self)
      class(blockstat), intent(inout) :: self
      self%blockLen = 0
   end subroutine

   subroutine bstat_reset(self)
      class(blockstat), intent(inout) :: self
      call self%allBlockStat%reset()
      call self%thisBlockStat%reset()
      call self%totalStat%reset()
      self%lastBlockMeanValue = 0
   end subroutine

   integer function bstat_blocklen(self)
      class(blockstat), intent(in) :: self
      bstat_blocklen = self%blockLen
   end function

   subroutine bstat_addData(self,data,wgt)
      class(blockstat), intent(inout) :: self
      real(r8), intent(in)         :: data
      real(r8), intent(in), optional :: wgt
      real(r8) mean
      if (.not.present(wgt)) then
         call self%thisBlockStat%add(data)
         call self%totalStat%add(data)
         if (self%thisBlockStat%localcount()==self%blockLen) then
            mean = self%thisBlockStat%localmean()
            call self%allBlockStat%add(mean)
            call self%thisBlockStat%reset()
            self%lastBlockMeanValue = mean
         end if
      else
         call self%thisBlockStat%add(data,wgt)
         call self%totalStat%add(data,wgt)
         if (self%thisBlockStat%localcount()==self%blockLen) then
            mean = self%thisBlockStat%localmean()
            call self%allBlockStat%add(mean)
            call self%thisBlockStat%reset()
            self%lastBlockMeanValue = mean
         end if
      endif
   end subroutine bstat_addData


! access functions

   real(r8) function bstat_localmean(self)
      class(blockstat), intent(in)  :: self
      bstat_localmean = self%totalStat%localmean()
   end function bstat_localmean

   real(r8) function bstat_locallastblockmean(self)
      class(blockstat), intent(in)  :: self
      bstat_locallastblockmean = self%lastBlockMeanValue
   end function bstat_locallastblockmean

   real(r8) function bstat_mean(self)
      class(blockstat), intent(in)  :: self
      bstat_mean = self%totalStat%mean()
   end function bstat_mean

   real(r8) function bstat_lastblockmean(self)
#ifdef MPI
      use MPI_F08
#endif
      class(blockstat), intent(in)  :: self
      real(r8) sendbuf,recvbuf
      integer nproc
#ifdef MPI
      integer ierr
#endif
      sendbuf = self%lastBlockMeanValue
      recvbuf = sendbuf  ! default for serial run
      nproc = 1
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
#endif
      bstat_lastblockmean = recvbuf / nproc
   end function bstat_lastblockmean

   real(r8) function bstat_localvar(self)
      class(blockstat), intent(in)  :: self
      bstat_localvar = self%totalStat%localvar()
   end function bstat_localvar


   real(r8) function bstat_var(self)
      class(blockstat), intent(in)  :: self
      bstat_var = self%totalStat%var()
   end function bstat_var


   real(r8) function bstat_localstddev(self)
      class(blockstat), intent(in)  :: self
      bstat_localstddev = self%allBlockStat%localstddev()
   end function bstat_localstddev


   real(r8) function bstat_stddev(self)
      class(blockstat), intent(in)  :: self
      bstat_stddev = self%allBlockStat%stddev()
   end function bstat_stddev


   integer(i8) function bstat_localcount(self)
      class(blockstat), intent(in)  :: self
      bstat_localcount = self%totalStat%localcount()
   end function bstat_localcount


   integer(i8) function bstat_count(self)
      class(blockstat), intent(in)  :: self
      bstat_count = self%totalStat%count()
   end function bstat_count


   integer(i8) function bstat_localblockcount(self)
      class(blockstat), intent(in)  :: self
      bstat_localblockcount = self%allBlockStat%localcount()
   end function bstat_localblockcount


   integer(i8) function bstat_blockcount(self)
      class(blockstat), intent(in)  :: self
      bstat_blockcount = self%allBlockStat%count()
   end function bstat_blockcount


   real(r8) function bstat_localweight(self)
      class(blockstat), intent(in)  :: self
      bstat_localweight = self%totalStat%localweight()
   end function bstat_localweight


   real(r8) function bstat_weight(self)
      class(blockstat), intent(in)  :: self
      bstat_weight = self%totalStat%weight()
   end function bstat_weight


!!!!!!!!! vstat methods

   subroutine vstat_create(self,n)
      class(vectorstat), intent(inout) :: self
      integer, intent(in)         :: n
      allocate(self%averageSum(n),self%varianceSum(n))
      self%vsize = n
      call vstat_reset(self)
   end subroutine

   logical function vstat_iscreated(self)
      class(vectorstat), intent(in) :: self
      vstat_iscreated = self%vsize > 0
   end function

   subroutine vstat_destroy(self)
      class(vectorstat), intent(inout) :: self
      deallocate(self%averageSum,self%varianceSum)
      self%vsize = 0
   end subroutine

   subroutine vstat_reset(self)
      class(vectorstat), intent(inout) :: self
      self%averageSum = 0
      self%varianceSum = 0
      self%weightSum = 0
      self%dataCounter = 0
   end subroutine

   subroutine vstat_addData(self,vdata,vwgt)
      class(vectorstat), intent(inout) :: self
      real(r8), intent(in)         :: vdata(:)
      real(r8), intent(in), optional :: vwgt
      if (.not.present(vwgt)) then
         self%averageSum = self%averageSum + vdata
         self%varianceSum = self%varianceSum + vdata**2
         self%weightSum = self%weightSum + 1.d0
         self%dataCounter = self%dataCounter + 1
      else
         self%averageSum = self%averageSum + vwgt*vdata
         self%varianceSum = self%varianceSum + vwgt*vdata**2
         self%weightSum = self%weightSum + vwgt
         self%dataCounter = self%dataCounter + 1
      endif
   end subroutine vstat_addData


! access functions


   function vstat_localmean(self) result(localmean)
      class(vectorstat), intent(in)  :: self
      real(r8)                    :: localmean(self%vsize)
      call assert(self%weightSum > 0d0,"vstat:mean no data or no weight")
      localmean = self%averageSum / self%weightSum
   end function vstat_localmean


   function vstat_mean(self) result(mean)
#ifdef MPI
      use MPI_F08
#endif
      class(vectorstat), intent(in)  :: self
      real(r8)                    :: mean(self%vsize)
      real(r8) sendbuf(self%vsize+1),recvbuf(self%vsize+1)
      integer taskid
#ifdef MPI
      integer ierr
#endif
      sendbuf(1:self%vsize) = self%averageSum; sendbuf(self%vsize+1) = self%weightSum
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,self%vsize+1,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0) then
         call assert(recvbuf(self%vsize+1) > 0,"vstat:mean no data")
         mean = recvbuf(1:self%vsize) / recvbuf(self%vsize+1)
      else
         mean = 0.d0
      end if
   end function vstat_mean


   function vstat_localvar(self) result(var)
      class(vectorstat), intent(in)  :: self
      real(r8)                    :: var(self%vsize)
      call assert(self%weightSum > 0,"vstat:variance no data or no weight")
      var = max(0.d0, self%varianceSum/self%weightSum - (self%averageSum / self%weightSum)**2)
   end function vstat_localvar


   function vstat_var(self) result(var)
#ifdef MPI
      use MPI_F08
#endif
      class(vectorstat), intent(in)  :: self
      real(r8)                    :: var(self%vsize)
      real(r8) sendbuf(2*self%vsize+1),recvbuf(2*self%vsize+1)
      integer taskid
#ifdef MPI
      integer ierr
#endif
      sendbuf(1:self%vsize) = self%averageSum
      sendbuf(self%vsize+1:2*self%vsize)=self%varianceSum
      sendbuf(2*self%vsize+1) = self%weightSum
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,2*self%vsize+1,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0) then
         call assert(recvbuf(2*self%vsize+1) > 0,"vstat:variance no data or no weight")
         var = max(0.d0, recvbuf(self%vsize+1:2*self%vsize) / recvbuf(2*self%vsize+1) &
               - ( recvbuf(1:self%vsize) / recvbuf(2*self%vsize+1) )**2 )
      else
         var = 0.d0
      end if
   end function vstat_var


   function vstat_localstddev(self) result(stddev)
      class(vectorstat), intent(in)  :: self
      real(r8)                    :: stddev(self%vsize)
      if (self%dataCounter > 1) then
         stddev = sqrt(vstat_localvar(self) / (self%dataCounter - 1))
      else
         stddev = 0.d0
      end if
   end function vstat_localstddev


   function vstat_stddev(self) result(stddev)
#ifdef MPI
      use MPI_F08
#endif
      class(vectorstat), intent(in)  :: self
      real(r8)                    :: stddev(self%vsize)
      real(r8) sendbuf(2*self%vsize+2),recvbuf(2*self%vsize+2)
      integer taskid
#ifdef MPI
      integer ierr
#endif
      sendbuf(1:self%vsize) = self%averageSum
      sendbuf(self%vsize+1:2*self%vsize)=self%varianceSum
      sendbuf(2*self%vsize+1) = self%weightSum
      sendbuf(2*self%vsize+2) = self%dataCounter
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,2*self%vsize+2,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0 .and. recvbuf(2*self%vsize+1) > 0 .and. recvbuf(2*self%vsize+2) > 1) then
         stddev = sqrt ( max(0.d0, ( recvbuf(self%vsize+1:2*self%vsize) / recvbuf(2*self%vsize+1)   &
                  - ( recvbuf(1:self%vsize) / recvbuf(2*self%vsize+1) )**2 ) / (recvbuf(2*self%vsize+2) - 1.d0) ) )
      else
         stddev = 0.d0
      end if
   end function vstat_stddev


   integer(i8) function vstat_localcount(self)
      class(vectorstat), intent(in)  :: self
      vstat_localcount = self%dataCounter
   end function vstat_localcount


   integer(i8) function vstat_count(self)
#ifdef MPI
      use MPI_F08
#endif
      class(vectorstat), intent(in)  :: self
      integer(i8) send,recv
#ifdef MPI
      integer ierr
#endif
      send = self%dataCounter
      recv = send  ! default for serial run
#ifdef MPI
      call MPI_REDUCE(send,recv,1,MPI_INTEGER8,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(recv,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,ierr)
#endif
      vstat_count = recv
   end function vstat_count


   function vstat_localweight(self) result(wgt)
      class(vectorstat), intent(in)  :: self
      real(r8)                    :: wgt
      wgt = self%weightSum
   end function vstat_localweight


   function vstat_weight(self) result(wgt)
#ifdef MPI
      use MPI_F08
#endif
      class(vectorstat), intent(in)  :: self
      real(r8)                    :: wgt
      real(r8) sendbuf,recvbuf
#ifdef MPI
      integer ierr
#endif
      sendbuf = self%weightSum
      recvbuf = sendbuf  ! default for serial run
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
#endif
      wgt = recvbuf
   end function vstat_weight


 !!!!!!!!! block vector stat methods

   subroutine bvstat_create(self,n,blockLen)
      class(blockvectorstat), intent(inout) :: self
      integer, intent(in) :: n
      integer, intent(in) :: blockLen
      call self%allBlockStat%create(n)
      call self%thisBlockStat%create(n)
      call self%totalStat%create(n)
      allocate(self%lastBlockMeanValue(n))
      self%blockLen = blockLen
   end subroutine

   logical function bvstat_iscreated(self)
      class(blockvectorstat), intent(in) :: self
      bvstat_iscreated = self%blockLen > 0
   end function

   subroutine bvstat_destroy(self)
      class(blockvectorstat), intent(inout) :: self
      call self%allBlockStat%destroy()
      call self%thisBlockStat%destroy()
      call self%totalStat%destroy()
      deallocate(self%lastBlockMeanValue)
      self%blockLen = 0
   end subroutine

   subroutine bvstat_reset(self)
      class(blockvectorstat), intent(inout) :: self
      call self%allBlockStat%reset()
      call self%thisBlockStat%reset()
      call self%totalStat%reset()
      self%lastBlockMeanValue = 0
   end subroutine

   integer function bvstat_blocklen(self)
      class(blockvectorstat), intent(in) :: self
      bvstat_blocklen = self%blockLen
   end function

   subroutine bvstat_addData(self,vdata,vwgt)
      class(blockvectorstat), intent(inout) :: self
      real(r8), intent(in)           :: vdata(:)
      real(r8), intent(in), optional :: vwgt
      real(r8) mean(size(vdata))
      if (.not.present(vwgt)) then
         call self%thisBlockStat%add(vdata)
         call self%totalStat%add(vdata)
         if (self%thisBlockStat%localcount()==self%blockLen) then
            mean = self%thisBlockStat%localmean()
            call self%allBlockStat%add(mean)
            call self%thisBlockStat%reset()
            self%lastBlockMeanValue = mean
         end if
      else
         call self%thisBlockStat%add(vdata,vwgt)
         call self%totalStat%add(vdata,vwgt)
         if (self%thisBlockStat%localcount()==self%blockLen) then
            mean = self%thisBlockStat%localmean()
            call self%allBlockStat%add(mean)
            call self%thisBlockStat%reset()
            self%lastBlockMeanValue = mean
         end if
      endif
   end subroutine bvstat_addData


! access functions

   function bvstat_localmean(self) result(localmean)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                    :: localmean(self%allBlockStat%vsize)
      localmean = self%totalStat%localmean()
   end function bvstat_localmean


   function bvstat_locallastblockmean(self) result(lbmean)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                    :: lbmean(self%allBlockStat%vsize)
      lbmean = self%lastBlockMeanValue
   end function bvstat_locallastblockmean


   function bvstat_mean(self) result(mean)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                    :: mean(self%allBlockStat%vsize)
      mean = self%totalStat%mean()
   end function bvstat_mean


   function bvstat_lastblockmean(self) result(lbmean)
#ifdef MPI
      use MPI_F08
#endif
      class(blockvectorstat), intent(in)  :: self
      real(r8)                    :: lbmean(self%allBlockStat%vsize)
      real(r8) sendbuf(self%allBlockStat%vsize),recvbuf(self%allBlockStat%vsize)
      integer nproc
#ifdef MPI
      integer ierr, n
#endif
      sendbuf = self%lastBlockMeanValue
      recvbuf = sendbuf  ! default for serial run
      nproc = 1
#ifdef MPI
      n = self%allBlockStat%vsize
      call MPI_REDUCE(sendbuf,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
#endif
      lbmean = recvbuf / nproc
   end function bvstat_lastblockmean


   function bvstat_localvar(self) result(var)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                              :: var(self%allBlockStat%vsize)
      var = self%totalStat%localvar()
   end function bvstat_localvar


   function bvstat_var(self) result(var)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                              :: var(self%allBlockStat%vsize)
      var = self%totalStat%var()
   end function bvstat_var


   function bvstat_localstddev(self) result(stddev)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                              :: stddev(self%allBlockStat%vsize)
      stddev = self%allBlockStat%localstddev()
   end function bvstat_localstddev


   function bvstat_stddev(self) result(stddev)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                              :: stddev(self%allBlockStat%vsize)
      stddev = self%allBlockStat%stddev()
   end function bvstat_stddev


   integer(i8) function bvstat_localcount(self)
      class(blockvectorstat), intent(in)  :: self
      bvstat_localcount = self%totalStat%localcount()
   end function bvstat_localcount


   integer(i8) function bvstat_count(self)
      class(blockvectorstat), intent(in)  :: self
      bvstat_count = self%totalStat%count()
   end function bvstat_count


   integer(i8) function bvstat_localblockcount(self)
      class(blockvectorstat), intent(in)  :: self
      bvstat_localblockcount = self%allBlockStat%localcount()
   end function bvstat_localblockcount


   integer(i8) function bvstat_blockcount(self)
      class(blockvectorstat), intent(in)  :: self
      bvstat_blockcount = self%allBlockStat%count()
   end function bvstat_blockcount


   function bvstat_localweight(self) result(w)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                              :: w
      w = self%totalStat%localweight()
   end function bvstat_localweight


   function bvstat_weight(self) result(w)
      class(blockvectorstat), intent(in)  :: self
      real(r8)                              :: w
      w = self%totalStat%weight()
   end function bvstat_weight



!!!!!!!! mstat methods

   subroutine mstat_create(self,m,n)
      class(matrixstat), intent(inout) :: self
      integer, intent(in)         :: m,n
      allocate(self%averageSum(m,n),self%varianceSum(m,n))
      self%vsize = m*n
      self%m = m; self%n = n
      call mstat_reset(self)
   end subroutine

   logical function mstat_iscreated(self)
      class(matrixstat), intent(in) :: self
      mstat_iscreated = self%vsize > 0
   end function

   subroutine mstat_destroy(self)
      class(matrixstat), intent(inout) :: self
      deallocate(self%averageSum,self%varianceSum)
   end subroutine

   subroutine mstat_reset(self)
      class(matrixstat), intent(inout) :: self
      self%averageSum = 0
      self%varianceSum = 0
      self%weightSum = 0
      self%dataCounter = 0
   end subroutine

   subroutine mstat_addData(self,vdata,vwgt)
      class(matrixstat), intent(inout) :: self
      real(r8), intent(in)         :: vdata(:,:)
      real(r8), intent(in), optional :: vwgt

      if (.not.present(vwgt)) then
         self%averageSum = self%averageSum + vdata
         self%varianceSum = self%varianceSum + vdata**2
         self%weightSum = self%weightSum + 1.d0
         self%dataCounter = self%dataCounter + 1
      else
         self%averageSum = self%averageSum + vwgt*vdata
         self%varianceSum = self%varianceSum + vwgt*vdata**2
         self%weightSum = self%weightSum + vwgt
         self%dataCounter = self%dataCounter + 1
      endif
   end subroutine mstat_addData


! access functions


   function mstat_localmean(self) result(localmean)
      class(matrixstat), intent(in)  :: self
      real(r8)                    :: localmean(self%m,self%n)
      call assert(self%weightSum > 0d0,"mstat:mean no data or no weight")
      localmean = self%averageSum / self%weightSum
   end function mstat_localmean


   function mstat_mean(self) result(mean)
#ifdef MPI
      use MPI_F08
#endif
      class(matrixstat), intent(in)  :: self
      real(r8)                    :: mean(self%m,self%n)
      real(r8) sendbuf(self%vsize+1),recvbuf(self%vsize+1)
      integer taskid
#ifdef MPI
      integer ierr
#endif
      sendbuf(1:self%vsize) = reshape( self%averageSum, (/ self%vsize /) )
      sendbuf(self%vsize+1) = self%weightSum
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,self%vsize+1,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0) then
         call assert(recvbuf(self%vsize+1) > 0,"mstat:mean no data")
         mean = reshape( recvbuf(1:self%vsize) / recvbuf(self%vsize+1), (/ self%m,self%n/) )
      else
        mean = 0.d0
      end if
   end function mstat_mean


   function mstat_localvar(self) result(var)
      class(matrixstat), intent(in)  :: self
      real(r8)                    :: var(self%m,self%n)
      call assert(self%weightSum > 0,"mstat:variance no data or no weight")
      var = max(0.d0, self%varianceSum/self%weightSum - (self%averageSum / self%weightSum)**2)
   end function mstat_localvar


   function mstat_var(self) result(var)
#ifdef MPI
      use MPI_F08
#endif
      class(matrixstat), intent(in)  :: self
      real(r8)                    :: var(self%m,self%n)
      real(r8) sendbuf(2*self%vsize+1),recvbuf(2*self%vsize+1)
      integer taskid
#ifdef MPI
      integer ierr
#endif

      sendbuf(1:self%vsize) = reshape( self%averageSum, (/ self%vsize /) )
      sendbuf(self%vsize+1:2*self%vsize) = reshape( self%varianceSum, (/ self%vsize /) )
      sendbuf(2*self%vsize+1) = self%weightSum
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,2*self%vsize+1,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0) then
         call assert(recvbuf(2*self%vsize+1)> 0,"mstat:variance no data or no weight")
         var = reshape( max( 0.d0, recvbuf(self%vsize+1:2*self%vsize) / recvbuf(2*self%vsize+1) &
               - ( recvbuf(1:self%vsize) / recvbuf(2*self%vsize+1) )**2 ) , (/ self%m, self%n /) )
      else
         var = 0.d0
      end if
   end function mstat_var


   function mstat_localstddev(self) result(stddev)
      class(matrixstat), intent(in)  :: self
      real(r8)                    :: stddev(self%m,self%n)
      if (self%dataCounter > 1) then
         stddev = sqrt(mstat_localvar(self) / (self%dataCounter - 1))
      else
         stddev = 0.d0
      end if
   end function mstat_localstddev


   function mstat_stddev(self) result(stddev)
#ifdef MPI
      use MPI_F08
#endif
      class(matrixstat), intent(in)  :: self
      real(r8)                    :: stddev(self%m,self%n)
      real(r8) sendbuf(2*self%vsize+2),recvbuf(2*self%vsize+2)
      integer taskid
#ifdef MPI
      integer ierr
#endif
      sendbuf(1:self%vsize) = reshape( self%averageSum, (/ self%vsize /) )
      sendbuf(self%vsize+1:2*self%vsize) = reshape( self%varianceSum, (/ self%vsize /) )
      sendbuf(2*self%vsize+1) = self%weightSum
      sendbuf(2*self%vsize+2) = self%dataCounter
      recvbuf = sendbuf  ! default for serial run
      taskid = 0
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,2*self%vsize+2,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
#endif
      if (taskid == 0 .and. recvbuf(2*self%vsize+1) > 0 .and. recvbuf(2*self%vsize+2) > 1) then
         stddev = reshape( &
               sqrt ( max(0.d0, ( recvbuf(self%vsize+1:2*self%vsize) / recvbuf(2*self%vsize+1) &
                  - ( recvbuf(1:self%vsize) / recvbuf(2*self%vsize+1) )**2 ) / (recvbuf(2*self%vsize+2) - 1.d0) ) ),  &
               (/ self%m, self%n /) )
      else
         stddev = 0.d0
      end if
   end function mstat_stddev


   integer(i8) function mstat_localcount(self)
      class(matrixstat), intent(in)  :: self
      mstat_localcount = self%dataCounter
   end function mstat_localcount


   integer(i8) function mstat_count(self)
#ifdef MPI
      use MPI_F08
#endif
      class(matrixstat), intent(in)  :: self
      integer(i8) send,recv
#ifdef MPI
      integer ierr
#endif
      send = self%dataCounter
      recv = send  ! default for serial run
#ifdef MPI
      call MPI_REDUCE(send,recv,1,MPI_INTEGER8,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(recv,1,MPI_INTEGER8,MASTER,MPI_COMM_WORLD,ierr)
#endif
      mstat_count = recv
   end function mstat_count


   function mstat_localweight(self) result(wgt)
      class(matrixstat), intent(in)  :: self
      real(r8)                    :: wgt
      wgt = self%weightSum
   end function mstat_localweight


   function mstat_weight(self) result(wgt)
#ifdef MPI
      use MPI_F08
#endif
      class(matrixstat), intent(in)  :: self
      real(r8)                    :: wgt
      real(r8) sendbuf,recvbuf
#ifdef MPI
      integer ierr
#endif
      sendbuf = self%weightSum
      recvbuf = sendbuf  ! default for serial run
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
#endif
      wgt = recvbuf
   end function mstat_weight


 !!!!!!!!! block matrix stat methods

   subroutine bmstat_create(self,m,n,blockLen)
      class(blockmatrixstat), intent(inout) :: self
      integer, intent(in) :: m,n
      integer, intent(in) :: blockLen
      call self%allBlockStat%create(m,n)
      call self%thisBlockStat%create(m,n)
      call self%totalStat%create(m,n)
      allocate(self%lastBlockMeanValue(m,n))
      self%blockLen = blockLen
   end subroutine

   logical function bmstat_iscreated(self)
      class(blockmatrixstat), intent(in) :: self
      bmstat_iscreated = self%blockLen > 0
   end function

   subroutine bmstat_destroy(self)
      class(blockmatrixstat), intent(inout) :: self
      call self%allBlockStat%destroy()
      call self%thisBlockStat%destroy()
      call self%totalStat%destroy()
      deallocate(self%lastBlockMeanValue)
      self%blockLen = 0
   end subroutine

   subroutine bmstat_reset(self)
      class(blockmatrixstat), intent(inout) :: self
      call self%allBlockStat%reset()
      call self%thisBlockStat%reset()
      call self%totalStat%reset()
      self%lastBlockMeanValue = 0
   end subroutine

   integer function bmstat_blocklen(self)
      class(blockmatrixstat), intent(in) :: self
      bmstat_blocklen = self%blockLen
   end function

   subroutine bmstat_addData(self,mdata,mwgt)
      class(blockmatrixstat), intent(inout) :: self
      real(r8), intent(in)           :: mdata(:,:)
      real(r8), intent(in), optional :: mwgt

      real(r8) mean(size(mdata,1),size(mdata,2))
      if (.not.present(mwgt)) then
         call self%thisBlockStat%add(mdata)
         call self%totalStat%add(mdata)
         if (self%thisBlockStat%localcount()==self%blockLen) then
            mean = self%thisBlockStat%localmean()
            call self%allBlockStat%add(mean)
            call self%thisBlockStat%reset()
            self%lastBlockMeanValue = mean
         end if
      else
         call self%thisBlockStat%add(mdata,mwgt)
         call self%totalStat%add(mdata,mwgt)
         if (self%thisBlockStat%localcount()==self%blockLen) then
            mean = self%thisBlockStat%localmean()
            call self%allBlockStat%add(mean)
            call self%thisBlockStat%reset()
            self%lastBlockMeanValue = mean
         end if
      endif
   end subroutine bmstat_addData


! access functions

   function bmstat_localmean(self) result(localmean)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                    :: localmean(self%allBlockStat%m,self%allBlockStat%n)
      localmean = self%totalStat%localmean()
   end function bmstat_localmean


   function bmstat_locallastblockmean(self) result(lbmean)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                    :: lbmean(self%allBlockStat%m,self%allBlockStat%n)
      lbmean = self%lastBlockMeanValue
   end function bmstat_locallastblockmean


   function bmstat_mean(self) result(mean)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                    :: mean(self%allBlockStat%m,self%allBlockStat%n)
      mean = self%totalStat%mean()
   end function bmstat_mean


   function bmstat_lastblockmean(self) result(lbmean)
#ifdef MPI
      use MPI_F08
#endif
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                    :: lbmean(self%allBlockStat%m,self%allBlockStat%n)
      real(r8) sendbuf(self%allBlockStat%vsize),recvbuf(self%allBlockStat%vsize)
      integer nproc,n
#ifdef MPI
      integer ierr
#endif
      n = self%allBlockStat%vsize
      sendbuf = reshape( self%lastBlockMeanValue , (/ n /) )
      recvbuf = sendbuf  ! default for serial run
      nproc = 1
#ifdef MPI
      call MPI_REDUCE(sendbuf,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
#endif
      lbmean = reshape( recvbuf / nproc , (/ self%allBlockStat%m,self%allBlockStat%n /) )
   end function bmstat_lastblockmean


   function bmstat_localvar(self) result(var)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                              :: var(self%allBlockStat%m,self%allBlockStat%n)
      var = self%totalStat%localvar()
   end function bmstat_localvar


   function bmstat_var(self) result(var)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                              :: var(self%allBlockStat%m,self%allBlockStat%n)
      var = self%totalStat%var()
   end function bmstat_var


   function bmstat_localstddev(self) result(stddev)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                              :: stddev(self%allBlockStat%m,self%allBlockStat%n)
      stddev = self%allBlockStat%localstddev()
   end function bmstat_localstddev


   function bmstat_stddev(self) result(stddev)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                              :: stddev(self%allBlockStat%m,self%allBlockStat%n)
      stddev = self%allBlockStat%stddev()
   end function bmstat_stddev


   integer(i8) function bmstat_localcount(self)
      class(blockmatrixstat), intent(in)  :: self
      bmstat_localcount = self%totalStat%localcount()
   end function bmstat_localcount


   integer(i8) function bmstat_count(self)
      class(blockmatrixstat), intent(in)  :: self
      bmstat_count = self%totalStat%count()
   end function bmstat_count


   integer(i8) function bmstat_localblockcount(self)
      class(blockmatrixstat), intent(in)  :: self
      bmstat_localblockcount = self%allBlockStat%localcount()
   end function bmstat_localblockcount


   integer(i8) function bmstat_blockcount(self)
      class(blockmatrixstat), intent(in)  :: self
      bmstat_blockcount = self%allBlockStat%count()
   end function bmstat_blockcount


   function bmstat_localweight(self) result(w)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                              :: w
      w = self%totalStat%localweight()
   end function bmstat_localweight


   function bmstat_weight(self) result(w)
      class(blockmatrixstat), intent(in)  :: self
      real(r8)                              :: w
      w = self%totalStat%weight()
   end function bmstat_weight



end module newStatistics_m



