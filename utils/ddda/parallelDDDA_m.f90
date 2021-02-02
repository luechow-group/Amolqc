! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module parallelDDDA_m
    use statistic_m, only: StatisticData_t
    use decorrelation_m, only: DecorrelationData_t, operator (+)
    use globalUtils_m, only: MASTER, nproc
#ifdef MPI
    use, intrinsic :: iso_c_binding,   only : c_ptr, c_f_pointer, c_associated
    use kinds_m, only: i4
    use MPI_F08

    implicit none

contains
        subroutine myMPIReduceDecorr(decorrD, decorrDtot)
            type(DecorrelationData_t), intent(in)   :: decorrD
            type(DecorrelationData_t), intent(out)  :: decorrDtot
            type(mpi_datatype)                    :: mpi_decorrelation
            type(mpi_op)                          :: mpi_decorr_sum
            type(mpi_status)                      :: status
            integer(i4)                        :: k

            call myMPIDecorrelationStruct(mpi_decorrelation, mpi_decorr_sum)
            !TODO real reduce not working with intel compiler
            ! call MPI_Reduce(decorrD, decorrDtot, 1, mpi_decorrelation, mpi_decorr_sum, 0, MPI_COMM_WORLD)

            if (.not. MASTER) then
                call MPI_SEND(decorrD, 1, mpi_decorrelation, 0, 1, MPI_COMM_WORLD)
            else
                decorrDtot = decorrD
                do k = 1, nproc - 1
                    call MPI_RECV(decorrD, 1, mpi_decorrelation, k, 1, MPI_COMM_WORLD, status)
                    decorrDtot = decorrDtot + decorrD
                end do
            end if
        end subroutine myMPIReduceDecorr

        function myMPIStatisticStruct() result(mpi_statistic)
            type(mpi_datatype)             :: mpi_statistic
            type(StatisticData_t)            :: stat
            integer(kind=MPI_Address_kind) :: disp(3)

            call MPI_Get_address(stat%nSamp, disp(1))
            call MPI_Get_address(stat%sum,   disp(2))
            call MPI_Get_address(stat%sqSum, disp(3))

            disp(3) = disp(3) - disp(1)
            disp(2) = disp(2) - disp(1)
            disp(1) = 0

            call MPI_TYPE_CREATE_STRUCT(3, [1,1,1], disp, [mpi_integer8, mpi_real8, mpi_real8], mpi_statistic)
            call MPI_TYPE_COMMIT(mpi_statistic)
        end function myMPIStatisticStruct

        subroutine myMPIDecorrelationStruct(mpi_decorrelation, mpi_decorr_sum)
            type(mpi_datatype), intent(inout) :: mpi_decorrelation
            type(mpi_op),       intent(inout) :: mpi_decorr_sum ! mpi_decorr_sum is a mpi operation
            type(mpi_datatype)                :: mpi_statistic, mpi_statarry, mpi_singledecorr
            type(DecorrelationData_t)           :: decorr, decorrArray(2)
            integer(kind=MPI_Address_kind)    :: disp(5), lb, extent
            integer                           :: i

            ! get mpi_statistic as new datatype
            mpi_statistic = myMPIStatisticStruct()

            call MPI_Get_address(decorr%stats(1), disp(1))
            call MPI_Get_address(decorr%stats(2), disp(2))
            lb = 0
            extent = disp(2) - disp(1)
            call MPI_TYPE_CREATE_RESIZED(mpi_statistic, lb, extent, mpi_statarry)
            call MPI_TYPE_COMMIT(mpi_statarry)

            !get mpi_singledecorr
            call MPI_Get_address(decorr%size,     disp(1))
            call MPI_Get_address(decorr%nSamp,    disp(2))
            call MPI_Get_address(decorr%wSamps,   disp(3))
            call MPI_Get_address(decorr%wSampsEx, disp(4))
            call MPI_Get_address(decorr%stats,    disp(5))

            do i = 2, 5
                disp(i) = disp(i) - disp(1)
            end do
            disp(1) = 0

            call MPI_TYPE_CREATE_STRUCT(5, [1,1,64,64,64], disp,&
                    [mpi_integer8, mpi_integer8, mpi_real8, mpi_logical, mpi_statarry],&
                    mpi_singledecorr)
            call MPI_TYPE_COMMIT(mpi_singledecorr)

            ! get mpi_decorrelation
            call MPI_Get_address(decorrArray(1), disp(1))
            call MPI_Get_address(decorrArray(2), disp(2))
            lb = 0
            extent = disp(2) - disp(1)
            call MPI_TYPE_CREATE_RESIZED(mpi_singledecorr, lb, extent, mpi_decorrelation)
            call MPI_TYPE_COMMIT(mpi_decorrelation)

            ! not compiling with mpiifort
            ! call MPI_OP_CREATE(decorrSum, .true., mpi_decorr_sum)
        end subroutine myMPIDecorrelationStruct

        subroutine decorrSum(invec, inoutvec, len, datatype)
            type(c_ptr), value                             :: invec, inoutvec
            integer                                        :: len, i
            type(mpi_datatype)                             :: datatype
            type(DecorrelationData_t), dimension(:), pointer :: invec_d, inoutvec_d

            call c_f_pointer(invec,    invec_d,    [len] )
            call c_f_pointer(inoutvec, inoutvec_d, [len] )


            do i = 1, len
                inoutvec_d(i) = invec_d(i) + inoutvec_d(i)
            end do

            nullify(inoutvec_d)
            nullify(invec_d)
        end subroutine decorrSum
#endif
end module parallelDDDA_m
