! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

program statisticsTests_p
#ifdef MPI
    use MPI_F08
#endif
    use intList_tm,       only: IntList_test
    use statistics_tm,     only: statistic_test, statistic_par_test
    use newStatistics_tm, only: newstatistics_test

    implicit none
#ifdef MPI
    call MPI_INIT()
#endif

    call IntList_test()
    call statistic_test()
    call statistic_par_test()
    call newstatistics_test()

#ifdef MPI
      call MPI_FINALIZE()
#endif
end program statisticsTests_p
