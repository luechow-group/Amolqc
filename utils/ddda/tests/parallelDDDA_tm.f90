! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module parallelDDDA_tm
#ifdef MPI
    use MPI_F08
    use parallelDDDA_m, only: myMPIStatisticStruct, myMPIDecorrelationStruct, myMPIReduceDecorr
#endif
    use statistic_m
    use decorrelation_m
    use error_m,  only: assert

    implicit none
contains
subroutine ParallelDdda_Test(passed)
    logical :: passed
#ifdef MPI
    integer n_procs      ! number of MPI processes, determined at start up time
    integer my_id        ! id in [0, .. nproc-1] of current MPI process
    integer p
    type(mpi_status) status

    type(Statistic) :: stat
    type(StatisticData_t) :: statD
    type(MPI_Datatype) :: mpi_statistic

    type(Decorrelation_t) :: decorr, decorr2
    type(DecorrelationData_t) :: decorrD, decorrD2, decorrD3
    type(MPI_Datatype) :: mpi_decorrelation

    type(mpi_op) :: mpi_decorr_sum

    call MPI_Init                                  ! initialize MPI
    call MPI_Comm_size(MPI_COMM_WORLD, n_procs)    ! set up MPI communicator return nprocs
    call MPI_Comm_rank(MPI_COMM_WORLD, my_id)      ! return task id of this process

    ! build and send statistic objects
    mpi_statistic = myMPIStatisticStruct()
    call stat%Create(statD)

    do p = 0, my_id
        call stat%Add_data(p * 1d0)
    end do

    if (my_id == 0) then
        write(*,*) "DDDA MPI test:"
        write(*,'(A40)', advance='no') [character(36) :: "TestSendStat:"]
    end if
    call MPI_BARRIER(MPI_COMM_WORLD)

    if (my_id > 0) then
        call MPI_SEND(statD, 1, mpi_statistic, 0, 0, MPI_COMM_WORLD)
    else
        do p = 1, n_procs - 1
            call MPI_Recv(statD, 1, mpi_statistic, p, 0, MPI_COMM_WORLD, status)
            !write(*,*) "   stat received from my_id=", p
            !call stat%write
        end do
    end if

    ! build and send decorrelation objects
    call myMPIDecorrelationStruct(mpi_decorrelation, mpi_decorr_sum)

    call decorr%Create(decorrD)
    do p = 0, my_id * 10
        call decorr%Add_data(p * 1d0)
    end do

    if (my_id == 0) then
        write(*,*) ' '//achar(27)//'[32m PASSED'//achar(27)//'[0m'
        write(*,'(A40)', advance='no') [character(36) :: "TestSendDecorr:"]
    end if
    call MPI_BARRIER(MPI_COMM_WORLD)

    if (my_id > 0) then
        call MPI_SEND(decorrD, 1, mpi_decorrelation, 0, 1, MPI_COMM_WORLD)
    else
        decorrD2 = decorrD

        do p = 1, n_procs - 1
            call MPI_Recv(decorrD3, 1, mpi_decorrelation, p, 1, MPI_COMM_WORLD, status)
            decorrD2 = decorrD2 + decorrD3
            !write(*,*) "decorr received from my_id=", p
        end do
        !call decorr%Create(decorrD2)
        !write(*,*) "summed decorr received:"
        !call decorr%write(6)
    end if

    if (my_id == 0) then
        write(*,*) ' '//achar(27)//'[32m PASSED'//achar(27)//'[0m'
        write(*,'(A40)', advance='no') [character(36) :: "TestReduceDecorr:"]
    end if
    call MPI_BARRIER(MPI_COMM_WORLD)

    call myMPIReduceDecorr(decorrD, decorrD3)

    if (my_id == 0) then
        call decorr%Create(decorrD3)
        ! write(6,*) "reduced decorr:"
        ! call decorr%write(6)
        call decorr2%Create(decorrD2)
        call assert(decorr ** decorr2)
    end if

    if (my_id == 0) then
        write(*,*) ' '//achar(27)//'[32m PASSED'//achar(27)//'[0m'
    end if

    call MPI_Finalize()
    call stat%destroy()
    call decorr%destroy()
    call decorr2%destroy()
#endif
end subroutine ParallelDdda_Test
end module parallelDDDA_tm