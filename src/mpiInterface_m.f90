! Copyright (C) 1996, 2006, 2012, 2018 Arne Luechow
! Copyright (C) 2013 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later


! collection of routines that encapsulate MPI
! for speed let compiler inline these routines

module mpiInterface_m
  use kinds_m, only: r8, i8
  use global_m
#ifdef MPI
  use MPI_F08
#endif
  use globalUtils_m, only: mytid, nproc, MASTER, PARALLEL_RUN

  implicit none

  interface myMPIBcastInteger
    module procedure myMPIBcastIntegerArray
    module procedure myMPIBcastIntegerScalar
  end interface myMPIBcastInteger

  interface myMPIBcastDouble
    module procedure myMPIBcastDoubleArray
    module procedure myMPIBcastDoubleScalar
  end interface myMPIBcastDouble

  interface myMPIReduceSumDouble
    module procedure myMPIReduceSumDoubleArray
    module procedure myMPIReduceSumDoubleScalar
  end interface myMPIReduceSumDouble

  interface myMPIReduceSumInteger
    module procedure myMPIReduceSumIntegerArray
    module procedure myMPIReduceSumIntegerScalar
  end interface myMPIReduceSumInteger

  interface myMPIAllReduceSumDouble
    module procedure myMPIAllReduceSumDoubleArray
    module procedure myMPIAllReduceSumDoubleScalar
  end interface myMPIAllReduceSumDouble

  interface myMPIAllReduceSumInteger
    module procedure myMPIAllReduceSumIntegerArray
    module procedure myMPIAllReduceSumIntegerScalar
  end interface myMPIAllReduceSumInteger

  interface myMPIGatherInteger
    module procedure myMPIGatherIntegerArray
    module procedure myMPIGatherIntegerScalar
  end interface myMPIGatherInteger

contains

subroutine myMPIInitialize(ierr)
  integer :: ierr, nn

#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, mytid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
#else
  ierr = 0
  mytid = 0
  nproc = 1
#endif
  if (mytid > 0) MASTER = .false.
  if (nproc > 0) PARALLEL_RUN = .true.
end subroutine myMPIInitialize


subroutine myMPIFinalize(ierr)
  integer :: ierr
  
#ifdef MPI
  call mpi_finalize(ierr)
#else
  ierr = 0
#endif
end subroutine myMPIFinalize


subroutine myMPIBarrier(ierr)
  integer :: ierr
  
#ifdef MPI
  call mpi_barrier(MPI_COMM_WORLD,ierr)
#else
  ierr = 0
#endif
end subroutine myMPIBarrier


subroutine myMPIBcastString(string,n)
  integer :: n, ierr
  character*(*) :: string

#ifdef MPI
  call mpi_bcast(string,n,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif
end subroutine myMPIBcastString


subroutine myMPIBcastIntegerArray(array,n)
  integer n,ierr
  integer array(n)

#ifdef MPI
  call mpi_bcast(array,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
end subroutine myMPIBcastIntegerArray

subroutine myMPIBcastIntegerScalar(array,n)
  integer n,ierr
  integer array

#ifdef MPI
  call mpi_bcast(array,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
end subroutine myMPIBcastIntegerScalar


subroutine myMPIBcastDoubleArray(array,n)
  integer n,ierr
  double precision array(n)

#ifdef MPI
  call mpi_bcast(array,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
end subroutine myMPIBcastDoubleArray

subroutine myMPIBcastDoubleScalar(array,n)
  integer n,ierr
  double precision array

#ifdef MPI
  call mpi_bcast(array,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
end subroutine myMPIBcastDoubleScalar


subroutine myMPIReduceSumDoubleArray(sendbuf,recvbuf,n)
  integer n,ierr
  real(r8) sendbuf(n),recvbuf(n)

#ifdef MPI
  call mpi_reduce(sendbuf,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
#endif
end subroutine myMPIReduceSumDoubleArray

subroutine myMPIReduceSumDoubleScalar(sendbuf,recvbuf,n)
  integer n,ierr
  real(r8) sendbuf,recvbuf

#ifdef MPI
  call mpi_reduce(sendbuf,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
#endif
end subroutine myMPIReduceSumDoubleScalar


subroutine myMPIReduceSumIntegerArray(sendbuf,recvbuf,n)
  integer n,ierr
  integer sendbuf(n),recvbuf(n)

#ifdef MPI
  call mpi_reduce(sendbuf,recvbuf,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
#endif
end subroutine myMPIReduceSumIntegerArray

subroutine myMPIReduceSumIntegerScalar(sendbuf,recvbuf,n)
  integer n,ierr
  integer sendbuf,recvbuf

#ifdef MPI
  call mpi_reduce(sendbuf,recvbuf,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
#endif
end subroutine myMPIReduceSumIntegerScalar


subroutine myMPIAllReduceSumDoubleArray(sendbuf,recvbuf,n)
  integer n,ierr
  real(r8) sendbuf(n),recvbuf(n)

#ifdef MPI
  call mpi_allreduce(sendbuf,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
#endif
end subroutine myMPIAllReduceSumDoubleArray

subroutine myMPIAllReduceSumDoubleScalar(sendbuf,recvbuf,n)
  integer n,ierr
  real(r8) sendbuf,recvbuf

#ifdef MPI
  call mpi_allreduce(sendbuf,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
#endif
end subroutine myMPIAllReduceSumDoubleScalar


subroutine myMPIAllReduceSumIntegerArray(sendbuf,recvbuf,n)
  integer n,ierr
  integer sendbuf(n),recvbuf(n)

#ifdef MPI
  call mpi_allreduce(sendbuf,recvbuf,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
#endif
end subroutine myMPIAllReduceSumIntegerArray

subroutine myMPIAllReduceSumIntegerScalar(sendbuf,recvbuf,n)
  integer n,ierr
  integer sendbuf,recvbuf

#ifdef MPI
  call mpi_allreduce(sendbuf,recvbuf,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
#endif
end subroutine myMPIAllReduceSumIntegerScalar


subroutine myMPIGatherIntegerArray(sendbuf,n,recvbuf,ierr)
  integer n,ierr
  integer sendbuf(n),recvbuf(nproc*n)

#ifdef MPI
  call MPI_GATHER(sendbuf,n,MPI_INTEGER,recvbuf,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
  ierr = 0
#endif
end subroutine myMPIGatherIntegerArray

subroutine myMPIGatherIntegerScalar(sendbuf,n,recvbuf,ierr)
  integer n,ierr
  integer sendbuf,recvbuf(nproc*n)

#ifdef MPI
  call MPI_GATHER(sendbuf,n,MPI_INTEGER,recvbuf,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
  ierr = 0
#endif
end subroutine myMPIGatherIntegerScalar

subroutine myMPIScatterInteger(sendbuf,n,recvbuf,ierr)
  integer n,ierr
  integer sendbuf(nproc*n),recvbuf(n)

#ifdef MPI
  call MPI_SCATTER(sendbuf,n,MPI_INTEGER,recvbuf,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
  ierr = 0
#endif
end subroutine myMPIScatterInteger


subroutine myMPIGatherDouble(sendbuf,n,recvbuf,ierr)
  integer n,ierr
  real(r8) sendbuf(n),recvbuf(nproc*n)

#ifdef MPI
  call MPI_GATHER(sendbuf,n,MPI_DOUBLE_PRECISION,recvbuf,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
  ierr = 0
#endif
end subroutine myMPIGatherDouble

! simple wrapper for MPI_GATHERV, which gathers a different amount of data from
! each process. this wrapper assumes that the data should be received in order,
! i.e. the "displacement" (array offset) in the receive buffer for a thread is
! the sum of all element counts of the previous threads.
!TODO: use integer(kind=MPI_Address_kind), however seems to be working atm
subroutine myMPIGatherDoubleV(sendbuf,n,recvbuf,recvcnt,ierr)
  integer n,ierr
  real(r8) sendbuf(n),recvbuf(nproc*n)
  integer displacements(nproc), recvcnt(nproc), displ, i
#ifdef MPI
  displ = 0
  do i = 1, nproc
    displacements(i) = displ
    displ = displ + recvcnt(i)
  enddo
  call MPI_GATHERV(sendbuf,n,MPI_DOUBLE_PRECISION,recvbuf,recvcnt,displacements,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
  ierr = 0
#endif
end subroutine myMPIGatherDoubleV

subroutine myMPIAllGatherDouble(sendbuf,n,recvbuf,ierr)
  integer n,ierr
  real(r8) sendbuf(n),recvbuf(nproc*n)

#ifdef MPI
  call MPI_ALLGATHER(sendbuf,n,MPI_DOUBLE_PRECISION,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
  ierr = 0
#endif
end subroutine myMPIAllGatherDouble


subroutine myMPIScatterDouble(sendbuf,n,recvbuf,ierr)
  integer n,ierr
  real(r8) sendbuf(nproc*n),recvbuf(n)

#ifdef MPI
  call MPI_SCATTER(sendbuf,n,MPI_DOUBLE_PRECISION,recvbuf,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
  recvbuf = sendbuf
  ierr = 0
#endif
end subroutine myMPIScatterDouble


subroutine myMPISendDouble(mpiV,n,id,tag)
  integer                :: n        ! size of array
  real(r8)                 :: mpiV(n)
  integer                :: id       ! id of receiving node
  integer                :: tag      ! message tag
  integer ierr

#ifdef MPI
  call MPI_SEND(mpiV,n,MPI_DOUBLE_PRECISION,id,tag,MPI_COMM_WORLD,ierr)
#endif
end subroutine myMPISendDouble

subroutine myMPIReceiveDouble(mpiV,n,id,tag)
  integer                :: n        ! size of array
  real(r8)                 :: mpiV(n)
  integer                :: id       ! id of sending node
  integer                :: tag      ! message tag

#ifdef MPI
  integer                :: ierr

  type(MPI_STATUS)       :: status
  call MPI_RECV(mpiV,n,MPI_DOUBLE_PRECISION,id,tag,MPI_COMM_WORLD,status,ierr)
#endif
end subroutine myMPIReceiveDouble


subroutine myMPISendInteger(mpiV,n,id,tag)
  integer                :: n        ! size of array
  integer                :: mpiV(n)
  integer                :: id       ! id of receiving node
  integer                :: tag      ! message tag

#ifdef MPI
  integer ierr

  call MPI_SEND(mpiV,n,MPI_INTEGER,id,tag,MPI_COMM_WORLD,ierr)
#endif
end subroutine myMPISendInteger

subroutine myMPIReceiveInteger(mpiV,n,id,tag)
  integer                :: n        ! size of array
  integer                :: mpiV(n)
  integer                :: id       ! id of sending node
  integer                :: tag      ! message tag

#ifdef MPI
  integer                :: ierr

  type(MPI_STATUS)       :: status
  call MPI_RECV(mpiV,n,MPI_INTEGER,id,tag,MPI_COMM_WORLD,status,ierr)
#endif
end subroutine myMPIReceiveInteger

subroutine myMPIabort(ierr)
  integer, intent(inout) :: ierr

#ifdef MPI
  integer errorcode
  call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
#else
  ierr = 1
#endif
end subroutine myMPIabort

real(r8) function myMPIWallTime()
#ifdef MPI
  myMPIWallTime = MPI_Wtime()
#else
  integer(i8) count,count_rate

  call system_clock(count,count_rate)
  myMPIWallTime = dble(count)/dble(count_rate)
#endif
end function myMPIWallTime

subroutine myMPI_Recv_DynamicReal8(buf, source, tag)
  real(r8), allocatable, intent(out) :: buf(:)
  integer, intent(in) :: source, tag
#ifdef MPI
  type(MPI_STATUS) :: status
  integer :: srcID, count

  call MPI_PROBE(source, tag, MPI_COMM_WORLD, status)
  call MPI_Get_count(status, MPI_REAL8, count)

  if (allocated(buf)) then
    deallocate(buf)
  end if

  allocate(buf(count))
  srcID = status%MPI_SOURCE

  call MPI_RECV(buf, count, MPI_REAL8, srcID, tag, MPI_COMM_WORLD, status)
#endif
end subroutine myMPI_Recv_DynamicReal8

subroutine myMPI_Recv_DynamicInteger4(buf, source, tag)
    integer, allocatable, intent(out) :: buf(:)
    integer, intent(in) :: source, tag
#ifdef MPI
    type(MPI_STATUS) :: status
    integer :: srcID, count

    call MPI_PROBE(source, tag, MPI_COMM_WORLD, status)
    call MPI_Get_count(status, MPI_INTEGER4, count)

    if (allocated(buf)) then
      deallocate(buf)
    end if

    allocate(buf(count))
    srcID = status%MPI_SOURCE

    call MPI_RECV(buf, count, MPI_INTEGER4, srcID, tag, MPI_COMM_WORLD, status)
#endif
end subroutine myMPI_Recv_DynamicInteger4

end module mpiInterface_m
