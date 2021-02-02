! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module globalUtils_m
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   ! output/logging unit
   ! MPI n_proc/my_task_id
   ! abort subroutine (needs to know MPI and output unit)

#ifdef MPI
   use MPI_F08
#endif

   implicit none

   integer      :: nproc = 1      ! number of parallel processors (1: serial)
   integer      :: mytid = 0      ! task id (0 master or serial, 1..nproc-1)
   integer      :: logmode = 2    ! output verbosity (0=none,1,2,..6)
   integer      :: iul = output_unit       ! stdout
   integer      :: iue = error_unit        ! stderr
   integer      :: iull = 100     ! logging unit for node specific output
   logical      :: MASTER = .true.! .true. if MPI Master or serial
   logical      :: PARALLEL_RUN = .false. ! .t. if parallel run (nproc > 1)

   logical      :: useLAlib = .true.   ! use external BLAS LAPACK

contains

   subroutine abortp(s)
      ! general abnormal thread termination
      character(*), intent(in) :: s
#ifdef MPI
      integer ierr, iexit
#endif

      write(iul,*) ' *** ABNORMAL TERMINATION *** '
      write(iul,*) ' from rank ',mytid
      write(iul,*) s
#ifdef MPI
      iexit = 1
      call MPI_ABORT(MPI_COMM_WORLD, iexit, ierr)
#endif
      error stop
   end subroutine abortp

   subroutine stopp(s)
      ! general normal thread termination, abnormal only in parallel runs
      character(*), intent(in) :: s
#ifdef MPI
      integer ierr, iexit

      write(iul,*) ' *** ABNORMAL TERMINATION *** '
      write(iul,*) ' from rank ',mytid
      write(iul,*) s
      iexit = 1
      call MPI_ABORT(MPI_COMM_WORLD, iexit, ierr)
      error stop
#else
      write(iul,*) s
      stop 0
#endif
   end subroutine stopp

   integer pure function getNprocs()
      getNprocs = nproc
   end function

   subroutine setNprocs(n)
      integer, intent(in) :: n
      nproc = n
   end subroutine

   integer pure function getMyTaskId()
      getMyTaskId = mytid
   end function

   subroutine getMyTaskIdChar(mytaskid)
      character(len=*), intent(inout) :: mytaskid
      character(len=8) fs
      integer l
      l = len(mytaskid)
      fs = "(i"//char(48+l)//"."//char(48+l)//")"
      write(mytaskid,fs) mytid
   end subroutine

   subroutine setMyTaskId(n)
      integer, intent(in) :: n
      mytid = n
   end subroutine

   integer pure function getIul()
      getIul = iul
   end function getIul

   subroutine setIul(log_unit)
      integer, intent(in) :: log_unit
      iul = log_unit
   end subroutine setIul

   integer pure function getIull()
      getIull = iull
   end function getIull

   subroutine setIull(log_unit)
      integer, intent(in) :: log_unit
      iull = log_unit
   end subroutine setIull

   integer pure function getLogmode()
      getLogmode = logmode
   end function getLogmode

   subroutine setLogMode(verb)
      integer, intent(in) :: verb
      logmode = verb
   end subroutine setLogMode

   integer pure function getAltLogmode()
      getAltLogmode = logmode
   end function getAltLogmode

   subroutine setAltLogMode(verb)
      integer, intent(in) :: verb
      logmode = verb
   end subroutine setALtLogMode


end module globalUtils_m



