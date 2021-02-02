! Copyright (C) 2013, 2015, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module subloop_m

   use kinds_m, only: r8
#ifdef NAG
   use, intrinsic :: f90_unix_io, only: flush
#endif
   use global_m
   use error_m
   use parsing_m, only: getNextBlock, getinta, finda
   use utils_m, only: getToken,getTimeString,readFileParallel,expandMacro,printHeader,getFormattedHeader
   use init_m, only: initAmolqc,initFiles,initGen,finalizeAmolqc
   use rWSample_m, only: RWSample, writeSampleCommand, recalculateSample
   use randomWalker_m, only: setNumberOfElectrons,setNumberOfCenters,setEpart
   use qmcSample_m, only: sample, initInitialWalker
   use qmc_m, only: qmc_run, qmc_init, initWalkerStat, initTrajectory
   use properties_m
   use waveFunction_m
   use elocTest_m
   use optDerivsTest_m
   use jastrow_m, only: jasChangeType
   use wfParameters_m, only: wfparams_change
   ! use ce_new_m
   use references_m, only: references_check, references_check1, references_optimize
   use psiMax_m
   use maxRawData_m, only: maxraw_init
   use maxAnalysis_m, only: maxana_init
   use maxBasins_m, only: maxbas_init
   use mpiInterface_m, only: myMPIBcastInteger, myMPIWallTime, myMPIAllReduceSumInteger

   implicit none

contains


subroutine subloop_callSubroutine(lines, nl, smpl, psimax_obj)
   character(len=120), intent(in) :: lines(:)
   integer, intent(in), optional  :: nl
   type(RWSample), intent(inout)  :: smpl
   type(psimax), intent(inout)    :: psimax_obj
   character(len=120) :: subName
   integer iflag
   logical :: tmp = .false.

   call getstra(lines,nl,'name=',subName,iflag)
   call subloop(subName, smpl, tmp, psimax_obj)
end subroutine subloop_callSubroutine



subroutine subloop(subname, smpl, exitSubLoop, psimax_obj)
!-------------------------------------------!
   ! subloop enters the event script loop under the name "subname" which is either
   ! a subroutine (i.e. globally stored series of blocks) or a macro command
   ! (i.e. blocks stored in a file)
   character(len=*), intent(in)          :: subname
   type(RWSample), intent(inout)         :: smpl
   logical, intent(out)                   :: exitSubloop
   type(psimax), intent(inout), optional :: psimax_obj

   integer, parameter :: MAXLINES=100

   character(len=120)          :: token
   character(len=MAXLEN)       :: inLines(MAXLINES)=''
   character(len=MAXLEN)       :: blockLines(MAXLINES)=''
   character(len=MAXLEN)       :: macroLines(MAXLINES)=''
   character(len=15)           :: si = ' =======>      '
   character(len=14)           :: sf = '      <======='
   character(len=180)          :: macropath,macrofile
   integer                     :: idx,nbl,nil,i,iuf,io,nnew,iflag,fileExistsInt,exitLoopInt,exitLoopSum
   integer                     :: loopIdx,loopIter,currentLoopIter,subIdx
   integer                     :: mLines
   real(r8)                      :: start,startCPU,endCPU
   real(r8)                      :: tstart,tstartCPU,tendCPU,sendbuf(1),recvbuf(1),t
   logical                     :: wfRead,found,exitLoop,fileExists,wout

   exitSubloop = .false.
   call getAmolqcPath(macropath)
   call assert(len(trim(macropath)) < 116,"amolqc: amolqc path exceeds definition")
   macropath = trim(macropath)//"/cmds"

   do subIdx=1,MAXSUBS
      if (subNames(subIdx)==subname) then
         if (logmode >= 2) write(iul,'(/2a)') ' ================>     calling subroutine ',trim(subname)
         exit
      end if
   end do
   if (subIdx>MAXSUBS) then ! not found in subroutine list, look for .cmd file
      macrofile = trim(macropath)//'/'//trim(subname)//'.cmd'
      fileExistsInt = 1
      if (MASTER) then
         inquire(file=macrofile,exist=fileExists)
         if (.not.fileExists) then
            fileExistsInt = 0
            write(iul,'(//2a/)') "ERROR in .in file: unknown subroutine ",trim(subname)
         else
            if (logmode >= 2) write(iul,'(/2a)') ' * * *  calling macro ',trim(subname)
         end if
      end if
      call myMPIBcastInteger(fileExistsInt,1)
      if (fileExistsInt==0) call abortp("unkown command in infile")
      call readFileParallel(mytid,macrofile,macroLines,mLines)
      nil = 0
      idx = 1
      call expandMacro(inLines,nil,macroLines,mLines,idx)
      if (logmode >= 2) then
         write(iul,'(//3a)') '  ---  Expanding macro cmd ',trim(subname),' to:'
         do i=1,mLines
            write(iul,'(A)') trim(inlines(idx+i-1))
            if ( index(inlines(idx+i-1), 'REQUIRED') > 0) then
               call abortp('Macro used without giving all REQUIRED arguments.')
            end if
         enddo
      endif
   else ! found in subroutine list
      do i=1,subLen(subIdx)
         inLines(i) = subLines(i,subIdx)
      end do
      nil = subLen(subIdx)
   end if


   idx = 1
   wout = .false.   ! write output
   if (MASTER .and. logmode > 1) wout = .true.
   do
      call getNextBlock(inLines,nil,idx,'$',')','!',MAXLINES,blockLines,nbl)

      if (nbl == 0) exit

      token = getToken(blockLines(1),'$','(')

      start = myMPIWallTime()
      call cpu_time(startCPU)

      if (token=='change_jastrow') then
         if (wout) write(iul,'(/a/)') si//'$change_jastrow - changing Jastrow terms'//sf
         call jasChangeType(blocklines,nbl)
         call recalculateSample(smpl)
      else if (token=='change_parameters') then
         if (wout) write(iul,'(/a/)') si//'$change_parameters'//sf
         call wfparams_change(blocklines,nbl)
      else if (token=='eloctest') then
         if (wout) write(iul,'(/a/)') si//'$eloctest - testing local energies'//sf
         call runEloctest(blockLines,nbl,smpl)
      else if (token=='init_rawdata_generation') then
         if (wout) write(iul,'(/a/)') si//'$init_rawdata_generation - initializing raw data generation'//sf
         call maxraw_init(blocklines,nbl)
      else if (token=='init_max_analysis') then
         if (wout) write(iul,'(/a/)') si//'$init_max_analysis - initializing maximum analysis'//sf
         call maxana_init(blocklines,nbl)
      else if (token=='init_basin_analysis') then
         if (wout) write(iul,'(/a/)') si//'$init_basin_analysis - initializing basin analysis'//sf
         call maxbas_init(blocklines,nbl)
      else if (token=='init_max_search') then
         if (wout) write(iul,'(/a/)') si//'$init_max_search - initializing maxima search'//sf
         if (present(psimax_obj)) then
            call psimax_obj%init(blocklines, nbl)
         else
            call error("subloop: psimax not available")
         end if
      else if (token=='init_walker') then
         if (wout) write(iul,'(/a/)') si//'$init_walker - setting an initial walker'//sf
         call initInitialWalker(blockLines,nbl)
      else if (token=='optimize_refs') then
         if (wout) write(iul,'(/a/)') si//'$optimize_refs - optimizing references'//sf
         call references_optimize(blockLines,nbl)
      else if (token=='optimtest') then
         if (wout) write(iul,'(/a/)') si//'$optimtest - testing the optimizers'//sf
         call optimizeTest(blockLines,nbl,smpl)
      else if (token=='print_results') then
         if (wout) write(iul,'(/a/)') si//'$print_results - printing stored results'//sf
         call global_printSavedResults(blockLines,nbl)
      else if (token=='props') then
         if (wout) write(iul,'(/a/)') si//'$props - calculting properties'//sf
         call propInit(blocklines,nbl)
      else if (token=='qmc') then
         if (wout) write(iul,'(/a/)') si//'$qmc - running a qmc calculation'//sf
         if (present(psimax_obj)) then
            call qmc_init(blockLines,nbl, smpl, psimax_obj)
            call qmc_run(smpl, psimax_obj)
         else
            call qmc_init(blockLines,nbl, smpl)
            call qmc_run(smpl)
         endif
      else if (token=='sample') then
         if (wout) write(iul,'(/a/)') si//'$sample - creating or modifying the walker sample'//sf
         call sample(blockLines,nbl,smpl)
      else if (token=='save_result') then
         if (wout) write(iul,'(/a/)') si//'$save_results - storing current results'//sf
         call global_saveResult(blockLines,nbl)
      else if (token=='sed') then
         if (wout) write(iul,'(/a/)') si//'$sed - running a sed calculation'//sf
         if (present(psimax_obj)) then
            call qmc_init(blockLines,nbl, smpl, psimax_obj)
            call qmc_run(smpl, psimax_obj)
         else
            call qmc_init(blockLines,nbl, smpl)
            call qmc_run(smpl)
         endif
      else if (token=='stop_if') then
         exitLoop = global_exitIf(blocklines,nbl)
         if (exitLoop) then
            call abortp("stop on condition in subloop_m.f90")
         end if
      else if (token=='test_balance') then
         call testBalance(smpl)
      else if (token=='write_sample') then
         if (wout) write(iul,'(/a/)') si//'$write_sample - writing the walker sample'//sf
         call writeSampleCommand(blockLines,nbl,smpl)
      ! loop commands
      else if (token=='begin_loop') then
         call assert(nbl==1,"begin_loop calls may consist of only one line")
         loopIdx = idx
         call getinta(blocklines,nbl,'count=',loopIter,iflag)
         currentLoopIter = 1
         call setCurrentLoopIdx(currentLoopIter)
         call assert(iflag==0,'begin_loop requires count value')
      else if (token=='end_loop') then
         if (loopIter > 1) then
            idx = loopIdx
            loopIter = loopIter - 1
            currentLoopIter = currentLoopIter + 1
            call setCurrentLoopIdx(currentLoopIter)
         else
            call setCurrentLoopIdx(1)
         endif
      else if (token=='exit_if') then
         exitLoop = global_exitIf(blocklines,nbl)
         exitLoopInt = 0
         exitLoopSum = 0
         if (exitLoop) then
            exitLoopInt = 1
         endif
         call myMPIAllReduceSumInteger(exitLoopInt,exitLoopSum,1)
         if (exitLoopSum >= 1) then
            if (finda(blocklines,nbl,'stop')) then
               call abortp("exit_if on condition and stop in subloop_m.f90")
            else
               !found = .false.
               if (MASTER .and. logmode>=2) then
                  call printHeader(iul, 'exit from subloop: continue', '-')
               end if
               exitSubLoop = .true.
            end if
         end if
      ! macro expansion
      else
         macrofile = trim(macropath)//'/'//trim(token)//'.cmd'
         fileExistsInt = 1
         if (MASTER) then
            inquire(file=macrofile,exist=fileExists)
            if (.not.fileExists) then
               fileExistsInt = 0
               write(iul,'(//2a/)') "ERROR in .in file: unknown command ",trim(token)
            end if
         end if
         call myMPIBcastInteger(fileExistsInt,1)
         if (fileExistsInt==0) call abortp("unkown command in infile")
         call readFileParallel(mytid,macrofile,macroLines,mLines)
         call assert(nbl==1,' currently only one-line macros allowed')
         call expandMacro(inLines,nil,macroLines,mLines,idx)
         if (logmode >= 2) then
            write(iul,'(//3a)') '  ---  Expanding macro cmd ',trim(token),' to:'
            do i=1,mLines
               write(iul,'(A)') trim(inlines(idx+i-1))
               if ( index(inlines(idx+i-1), 'REQUIRED') > 0) then
                  call abortp('Macro used without giving all REQUIRED arguments.')
               end if
            enddo
         endif
      endif

      if (logmode >= 2) then
         if (token=='qmc' .or. token=='sample') then
            write(iul,'(/3A,F18.2,A)') ' wall clock time for   ',trim(token),' : ', &
              myMPIWallTime()-start,' s'
            call cpu_time(endCPU)
            write(iul,'(3A,F18.2,A//)') ' cpu time (master) for ',trim(token),' : ', &
            endCPU-startCPU,' s'
         end if
         call flush(iul)
      end if

   end do

   if (logmode >= 2) then
      if (subIdx>MAXSUBS) then
         write(iul,'(/2a/)') ' =============>    end macro ',trim(subname)
      else
         write(iul,'(/2a/)') ' =============>    end subroutine ',trim(subname)
      end if
   end if


end subroutine subloop

end module subloop_m
