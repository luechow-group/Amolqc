! Copyright (C) 2006-2008, 2012-2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

MODULE init_m

  use kinds_m, only: r8
  use compilerStrings_m, only: getCompilerOptions, getCompilerVersion
  use global_m
  use error_m
  use mpiInterface_m, only: myMPIInitialize, myMPIfinalize
  use machine_m, only: myGetArgLocal,myGetDate,myGetHost,mygetEnv
  use random_m, only: init_ran
  use utils_m, only: readFileParallel,printHeader,getFormattedHeader
#ifdef MPI
  use MPI_F08, only: MPI_GET_LIBRARY_VERSION, MPI_MAX_LIBRARY_VERSION_STRING
#endif

  implicit none

CONTAINS

!-----------------------!
subroutine initAmolqc()
!-----------------------!

  character(len=6) :: name
  character(len=180) :: path
  integer ierr,nn

  call myMPIInitialize(ierr)
  call getEnvironmentVariableName(name)
  call myGetEnv(name,path)
  call setAmolqcPath(path)
  call initSaveResults()

end subroutine initAmolqc


!-------------------------!
subroutine finalizeAmolqc()
!-------------------------!


  integer ierr
  character(len=26) :: date

  if (MASTER) then
     call myGetDate(date)
     write(iul,'(/2A//A//)') 'Amolqc run finished on ',date,'Bye!'
     close(iul)
  endif

  ! Properly finish MPI
  call myMPIfinalize(ierr)
end subroutine finalizeAmolqc



!--------------------------------!
subroutine initFiles(lines,nl)
!--------------------------------!


  character(len=*), intent(inout) :: lines(:)
  integer, intent(out)            :: nl
  character(len=40)               :: date,host
  character(len=4)                :: outtype = ".out"
  character(len=4)                :: intype = ".in"
  character(len=180)              :: path
  integer counter, i
  logical fileExists
  character(len=2) cntstr

  integer id,io,iflag,ierr,nn,length

#ifdef MPI
  character (len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpiVersion
  integer :: nameLength

  call MPI_GET_LIBRARY_VERSION(mpiVersion, nameLength)
#endif

  nl = 0

  call myGetArgLocal(1,baseName,ierr)
  if (MASTER) then
      if (ierr /= 0) call abortp(" amolqc requires a command line argument (name of .in file)")
      ! replacing 'string.in' with 'string'
      length = scan(baseName,' ') - 1
      if (length > 3) then
          if (baseName(length-2:length) == '.in') then
              baseName = baseName(1:length-3)
          else if (baseName(length-3:length) == '.ami') then
              baseName = baseName(1:length-4)
              intype = '.ami'
              outtype = '.amo'
          end if
      end if
      inquire(file=trim(baseName)//trim(intype),exist=fileExists)
      if (.not. fileExists) then
         intype = '.ami'
         outtype = '.amo'
         inquire(file=trim(baseName)//trim(intype),exist=fileExists)
         if (.not. fileExists) call abortp('input file '//trim(baseName)//trim(intype)//&
                                           ' or '//trim(baseName)//'.in not found')
      end if

      iul=9
      iull=100
      inquire(file=trim(baseName)//outtype,exist=fileExists)
      if (fileExists) then
         ! do not overwrite but find new out file name
         do counter=1,99
            write(cntstr,'(I2)') counter
            cntstr = adjustl(cntstr)
            inquire(file=trim(baseName)//outtype//'-'//trim(cntstr),exist=fileExists)
            if (.not.fileExists) exit
         end do
         open(iul,file=trim(baseName)//outtype//'-'//trim(cntstr),status='new',iostat=io)
         call assert(io==0,'(initFiles): opening output file '//trim(baseName)//outtype//'-'//trim(cntstr)//' failed')
      else
         open(iul,file=trim(baseName)//outtype,status='new',iostat=io)
         call assert(io==0,'(initFiles): opening output file '//trim(baseName)//outtype//' failed')
      end if

      call myGetDate(date)
      call myGetHost(host)

      write(iul,'(/)')
      write(iul,*) '             __  __    ____    _         ____     _____   '
      write(iul,*) '     /\     |  \/  |  / __ \  | |       / __ \   / ____|  '
      write(iul,*) '    /  \    | \  / | | |  | | | |      | |  | | | |       '
      write(iul,*) '   / /\ \   | |\/| | | |  | | | |      | |  | | | |       '
      write(iul,*) '  / ____ \  | |  | | | |__| | | |____  | |__| | | |____   '
      write(iul,*) ' /_/    \_\ |_|  |_|  \____/  |______|  \___\_\  \_____|  '

      write(iul,'(//A/)') ' Atoms and Molecules with Quantum Monte Carlo -- electron structure code'
      write(iul,'(A)')    ' initial version:'
      write(iul,'(A/)')   '  Arne Luechow, Penn State University, 2/1996'
      write(iul,'(A)')    ' main author:'
      write(iul,'(A/)')   '  Arne Luechow, RWTH Aachen University, 52056 Aachen, Germany'
      write(iul,'(A)')    ' with contributions from:'
      write(iul,'(A)')    '  Sebastian Manten, Christian Diedrich, Annika Bande, Tony C. Scott,'
      write(iul,'(A)')    '  Annett Schwarz, Rene Petz, Raphael Berner, Alexander Sturm,'
      write(iul,'(A)')    '  Marko Hermsen, Kaveh Haghighi Mood, Christoph Schulte,'
      write(iul,'(A)')    '  Leonard Reuter, Michael A. Heuer, Jil Ludovicy'
 
#ifdef VERSION
      write(iul,'(//2A)')  ' version:          ', VERSION
#else
      write(iul,'(//2A)')  ' version:          ', 'UNKNOWN'
#endif
      write(iul,'(2A)')    ' compiler version: ',getCompilerVersion()
      write(iul,'(2A)')    ' compiler options: ',getCompilerOptions()
#ifdef MPI
      write(iul,'(2A)')    ' mpi version:      ',trim(mpiVersion)
#endif
      write(iul,'(//5A,I4,A)') ' run started on ',trim(host),' at ',trim(date), &
          ' on ',nproc,' processor(s)'
      call getAmolqcPath(path)
      write(iul,'(2A)') ' using path: ',trim(path)

  else
      iull = 100+getMyTaskId()
      iul = iull
  endif

  ! Attention: slave processes may have the wrong baseName!
  ! only MASTER reads input file
  call readFileParallel(mytid,trim(baseName)//trim(intype),lines,nl)

  ! write input in .out
  if (MASTER) then
     call printHeader(iul, 'reading input')
     do i = 1, len(lines)
        if (lines(i) /= '') then
           write(iul,'(A)') trim(lines(i))
        end if
     end do
  end if
end subroutine initFiles



!--------------------------!
subroutine initGen(lines,nl)
!--------------------------!

  ! read and initialize $gen section

  integer                     :: nl
  character(len=120)          :: lines(nl)
  integer                     :: seed,iflag

  real(r8)                      :: dummy

  ! get seed and initialize random number generator

  call getinta(lines,nl,'seed=',seed,iflag)
  if (iflag /= 0 .or. seed <= 0) call abortp('$gen: seed > 0 required')
  seed = abs(seed) + mytid
  dummy = init_ran(seed)

  call getinta(lines,nl,'verbose=',logmode,iflag)
  if (.not. MASTER .and. logmode<3) logmode=0
  if (logmode>=2) then
     call printHeader(iul, getFormattedHeader('gen', 'initializing RNG and setting general parameters'))
     write(iul,'(A,I7,A,I2)') ' seed =',seed,'     verbose level =',logmode
  endif

end subroutine initGen

end module init_m

