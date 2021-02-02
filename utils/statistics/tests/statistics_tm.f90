! TODO UNKNOWN ORIGIN AND AUTHOR (OLDER THAN 2010?) (see amolqc2010-0917/utils/statisticst.f90)
! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module statistics_tm
  use kinds_m, only: r8
  use statistics
#ifdef MPI
  use MPI_F08
#endif

contains
  subroutine statistic_test
    ! this tests generic interface to the many modules in statistics_m.f90
    type(simpleStat)      :: stat
    type(weightStat)      :: wstat
    type(IntHistogram)    :: ihist
    type(DoubleHistogram) :: dhist
    real(r8)                :: xMin=0.d0,xMax=100.d0
    integer               :: iMin=0,iMax=100,nBins=10
    integer               :: i
    integer, pointer      :: ih(:,:)

    call reset(stat)
    call reset(wstat)

    do i=1,5
       call addData(stat,dble(i))
       call addData(wstat,dble(i),1.d0)
    enddo

    call assertEqualAbsolute(mean(stat), 3._r8, 1e-12_r8, "stattest-1")
    call assertEqualAbsolute(mean(wstat), 3._r8, 1e-12_r8, "stattest-2")

    call assert(dataCount(stat) == 5, "stattest-3")
    call assert(dataCount(wstat) == 5, "stattest-4")

    call createHistogram(ihist,iMin,iMax,nBins)
    call createHistogram(dhist,xMin,xMax,nBins)

    do i=0,100
       call addData(ihist,i)
       call addData(dhist,dble(i))
    enddo

    call assertEqualAbsolute(mean(ihist), 50._r8, 1e-12_r8, "stattest-5")
    call assertEqualAbsolute(mean(dhist), 50._r8, 1e-12_r8, "stattest-6")

    call assert(dataCount(ihist) == 101, "stattest-7")
    call assert(dataCount(dhist) == 101, "stattest-8")

  end subroutine statistic_test

  subroutine statistic_par_test()
    integer, parameter :: MASTER=0
    integer ierr,mytid,nproc,i, dataCount
    real(r8)  meanAll
    type(simpleStat)  :: stat

    mytid = 0
    nproc = 1
#ifdef MPI
    ! MPI initialization
    call mpi_comm_rank(MPI_COMM_WORLD, mytid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
#endif

    call reset(stat)

    call assert(nproc==1 .or. nproc==2, "statistic_par_test only with 1 or 2 processes")

    do i=1,100
      call addData(stat,mytid+1.0d-3*dble(i))
    enddo

   meanAll = meanAllNodes(stat)
   dataCount = dataCountAllNodes(stat)
   if (nproc == 1) then
      call assertEqualAbsolute(meanAll, 0.0505_r8, 1.e-12_r8, "stattest-par-1")
      call assert(dataCount == 100, "stattest-par-2")
   else
      if (mytid == MASTER) then
         call assertEqualAbsolute(meanAll, 0.5505_r8, 1.e-12_r8, "stattest-par-3")
         call assert(dataCount == 200, "stattest-par-4")
      end if
   end if
  end subroutine statistic_par_test
end module statistics_tm
