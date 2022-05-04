! Copyright (C) 2013-2015, 2018 Arne Luechow
! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module mainLoop_m

private
public :: mainloop

contains

subroutine mainloop

   use kinds_m, only: r8
#ifdef NAG
   use, intrinsic :: f90_unix_io, only: flush
#endif
   use global_m
   use error_m
   use parsing_m, only: getNextBlock, getinta, finda
   use utils_m, only: getToken,getTimeString,readFileParallel,expandMacro,printHeader,getFormattedHeader
   use init_m, only: initAmolqc,initFiles,initGen,finalizeAmolqc
   use rWSample_m, only: RWSample, writeSampleCommand, compareSample, recalculateSample, &
       getFirst, isNext, getNext, getSampleSize
   use randomWalker_m, only: RandomWalker, setNumberOfElectrons, setNumberOfCenters, setEpart
   use qmcSample_m, only: sample, initInitialWalker, displaceInitialWalker
   use qmc_m, only: qmc_run, qmc_init, initWalkerStat, initTrajectory
   use optimizeParams_m
   use properties_m
   use eloc_m, only: wf_init, wf_ecp_init
   use elocTest_m
   use optDerivsTest_m
   use jastrow_m, only: jasChangeType
   use wfParameters_m, only: wfparams_change
   ! use ce_new_m
   use references_m, only: references_check, references_check1, references_optimize
   use psiMax_m, only: psimax
   use rhoMax_m, only: rhoMax_t
   use moMax_m, only: moMax_run, moMax_plot, moMax_plot_plane
   use maximizeSample_m, only: maximizeSample, maximizeWalker
   use maximizeSampleRho_m, only: maximizeSampleRho, maximizeWalkerRho
   use maxRawData_m, only: maxraw_init
   use maxAnalysis_m, only: maxana_init
   use maxBasins_m, only: maxbas_init
   use subloop_m
   use mpiInterface_m, only: myMPIReduceSumInteger, myMPIBcastInteger,&
           myMPIReduceSumDouble, myMPIWallTime
   use electronDensity_m, only: calculateDensity, scanBondDensity
   use vxc_m, only: CalcVxc
   use rhoGrid_m, only: RhoGrid_t, InitRhoGrid
   use coulombDensity_m, only: coulomb_density

   implicit none

   integer, parameter :: MAXLINES=1000

   character(len=120)          :: token
   character(len=MAXLEN)       :: inLines(MAXLINES)=''
   character(len=MAXLEN)       :: blockLines(MAXLINES)=''
   character(len=MAXLEN)       :: macroLines(MAXLINES)=''
   character(len=180)          :: macropath,macrofile,subName
   type(RWSample)              :: smpl
   type(psimax)                :: psimax_obj
   type(rhoMax_t)              :: rhoMax
   type(RhoGrid_t)             :: rhoGrid
   integer                     :: bs,idx,nbl,nil,i,iuf,io,nnew,iflag,fileExistsInt
   integer                     :: loopIdx,loopIter,currentLoopIter,subIdx,subLine
   integer                     :: mLines
   real(r8)                    :: start,startCPU,endCPU
   real(r8)                    :: tstart,tstartCPU,tendCPU,sendbuf(1),recvbuf(1),t
   logical                     :: wfRead,found,exitLoop,fileExists,subMode,converged,wout, is_converged

   call initFiles(inLines,nil)
   call getAmolqcPath(macropath)
   call assert(len(trim(macropath)) < 116,"amolqc: amolqc path exceeds definition")
   macropath = trim(macropath)//"/cmds"

   tstart = myMPIWallTime()
   call cpu_time(tstartCPU)

   idx = 1
   wfRead = .false.
   converged = .true.
   subIdx = 0
   wout = .false.   ! write output
   do

      if (.not.converged) then
         if (MASTER) call printHeader(iul, 'not converged: terminating run', '-')
         exit
      end if

      call getNextBlock(inLines,nil,idx,'$',')','!',MAXLINES,blockLines,nbl)

      if (nbl == 0) exit

      token = getToken(blockLines(1),'$','(')

      start = myMPIWallTime()
      call cpu_time(startCPU)

      select case(token)
      case('gen')
         if (wfRead) call abortp('$gen block must precede $wf in .in file')
         call initGen(blockLines,nbl)
         if (MASTER .and. logmode>1) wout = .true.
      case('wf')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'wave function'))
         call wf_init(blockLines,nbl)
         if (.not.wfRead) then
            call setNumberOfElectrons(ne)
            call setNumberOfCenters(ncenter)
            call setEpart(do_epart)
            wfRead = .true.
         endif
      case('ecp')
         if (.not.wfRead) call abortp('$ecp block must follow $wf in .in file')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'effective core potential settings'))
         call wf_ecp_init(blockLines,nbl)
      case('analyze_refs')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'analyze references'))
         if (finda(blockLines,nbl,'max_mode=')) then
           call references_check1(blockLines,nbl)
         else
           call references_check(blockLines,nbl)
         end if
      case('calculate_density')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'get the determinant density at position'))
         call calculateDensity(blockLines, nbl)
      case('call_subroutine')
         call subloop_callSubroutine(blocklines, nbl, smpl, psimax_obj)
      case('change_jastrow')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'changing Jastrow terms'))
         call jasChangeType(blocklines,nbl)
         call recalculateSample(smpl)
      case('change_parameters')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'changing parameters'))
         call wfparams_change(blocklines,nbl)
      case('compare_sample')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'comparing samples'))
         call compareSample(blocklines, nbl, smpl)
      case('coulomb_density')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'coulomb_density'))
         call coulomb_density(blocklines, nbl)
      case('displace_walker')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'displace_walker'))
         call displaceInitialWalker(blocklines, nbl)
      case('eloctest')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'testing local energies'))
         call runEloctest(blockLines,nbl,smpl)
      case('init_rawdata_generation')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'initializing raw data generation'))
         call maxraw_init(blocklines,nbl)
      case('init_max_analysis')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'initializing maximum analysis'))
         call maxana_init(blocklines,nbl)
      case('init_basin_analysis')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'initializing basin analysis'))
         call maxbas_init(blocklines,nbl)
      case('init_max_search')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'initializing maxima search'))
         call psimax_obj%init(blocklines,nbl)
      case('init_rho_analysis')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'initializing electronic density analysis'))
         call rhoMax%init(blocklines,nbl)
      case('init_rho_grid')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'initializing electronic density grid'))
         call InitRhoGrid(blocklines,nbl,smpl,rhoGrid)
      case('init_walker')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'setting an initial walker'))
         call initInitialWalker(blockLines,nbl)
      case('maximize_sample')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'get psi^2 maxima of sample'))
         call maximizeSample(smpl, psimax_obj, wout)
      case('maximize_sample_rho')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'get rho maxima of sample'))
         call maximizeSampleRho(smpl, rhoMax)
      case('maximize_walker')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'get psi^2 maxima of walker'))
         call maximizeWalker(blockLines, nbl, smpl, psimax_obj, wout)
      case('maximize_walker_rho')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'get rho maxima of walker'))
         call maximizeWalkerRho(blockLines, nbl, smpl, rhoMax)
      case('maximize_mos')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'get maxima/minima of mos'))
         call moMax_run(blockLines, nbl)
      case('optimize_parameters')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'optimizing wave function parameters'))
         call optimizeParameters(blockLines,nbl,smpl,converged)
      case('optimize_refs')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'optimizing references'))
         call references_optimize(blockLines,nbl)
      case('plot_mos')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'write MO values on a grid'))
         call moMax_plot(blockLines, nbl)
      case('plot_mo_in_plane')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'write MO values on a grid with plane vectors'))
         call moMax_plot_plane(blockLines, nbl)
      case('print_results')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'printing stored results'))
         call global_printSavedResults(blockLines,nbl)
      case('props')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'calculting properties'))
         call propInit(blocklines,nbl)
      case('qmc')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'running a qmc calculation'))
         call qmc_init(blockLines,nbl,smpl, psimax_obj)
         call qmc_run(smpl, psimax_obj, rhoMax)
      case('sample')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'creating or modifying the walker sample'))
         call sample(blockLines,nbl,smpl)
      case('save_result')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'storing current results'))
         call global_saveResult(blockLines,nbl)
      case('scan_bond_density')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'scans the electronic density on a bond'))
         call scanBondDensity(blockLines, nbl)
      case('sed')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'running a sed calculation'))
         call qmc_init(blockLines, nbl, smpl, psimax_obj)
         call qmc_run(smpl, psimax_obj)
      case('stop_if')
         exitLoop = global_exitIf(blocklines,nbl)
         if (exitLoop) then
            if (MASTER) then
               call printHeader(iul, 'stop on condition', '-')
            end if
            exit
         end if
      case('test_balance')
         call testBalance(smpl)
      case('trajectory')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'initializing a trajectory run'))
         call initTrajectory(blockLines,nbl)
      case('vxc')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'calculating Vxc'))
         call CalcVxc(blockLines,nbl,smpl,rhoMax,rhoGrid)
      case('walker_stat')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'initializing walker statistics'))
         call initWalkerStat(blockLines,nbl)
      case('wf_param_deriv_test')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'testing the parameter derivatives'))
         call optimizeTest(blockLines,nbl,smpl)
      case('write_sample')
         if (wout) call printHeader(iul, getFormattedHeader(token, 'writing the walker sample'))
         call writeSampleCommand(blockLines,nbl,smpl)

      ! subroutine definition   
      case('begin_subroutine')
         call getstra(blockLines,nbl,'name=',subName,iflag)
         if (iflag /= 0) call abortp('begin_subroutine requires name argument')
         subIdx = subIdx + 1
         subNames(subIdx) = subName
         subLine = 0
         do 
            if (subLine > MAXSUBLINES) call abortp('subroutine has too many lines')
            token = getToken(inLines(idx+subLine),'$','(')
            if (token == 'end_subroutine') then
               idx = idx + subLine + 1
               subLen(subIdx) = subLine
               exit
            end if
            subLine = subLine + 1
            subLines(subLine,subIdx) = inLines(idx+subLine-1)
         end do
         if (wout) call printHeader(iul, getFormattedHeader('subroutine', 'storing subroutine' // trim(subName)))
      
      ! loop commands
      case('begin_loop')
         call assert(nbl==1,"begin_loop calls may consist of only one line")
         loopIdx = idx
         call getinta(blocklines,nbl,'count=',loopIter,iflag)
         currentLoopIter = 1
         call setCurrentLoopIdx(currentLoopIter)
         call assert(iflag==0,'begin_loop requires count value')
      case('end_loop')
         if (loopIter > 1) then
            idx = loopIdx
            loopIter = loopIter - 1
            currentLoopIter = currentLoopIter + 1
            call setCurrentLoopIdx(currentLoopIter)
         else
            call setCurrentLoopIdx(1)
         endif
      case('exit_if')
         exitLoop = global_exitIf(blocklines,nbl)
         if (exitLoop) then
            if (finda(blocklines,nbl,'stop')) then
               if (MASTER) then
                  call printHeader(iul, 'exit from loop: stop', '-')
               end if
               exit
            else 
               found = .false.
               if (MASTER .and. logmode>=2) then
                  call printHeader(iul, 'exit from loop: continue', '-')
               end if
               do i=idx,nil
                  if (index(inLines(i),'$end_loop') > 0) then
                     idx = i+1
                     found = .true.
                     exit
                  end if
               end do
               if (.not.found) call abortp("$exit_if only allowed between $begin_loop and $end_loop")
            end if
         end if

      ! macro expansion
      case default
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
            call printHeader(iul, getFormattedHeader(token, 'expanding the macro'))
            do i=1,mLines
               write(iul,'(A)') trim(inlines(idx+i-1))
               if ( index(inlines(idx+i-1), 'REQUIRED') > 0) then
                  call abortp('Macro used without giving all REQUIRED arguments.')
               end if
            enddo
         endif
      end select

      if (logmode >= 2) then
         if (token=='qmc' .or. token=='sample' .or. token=='optimize_parameters') then
            write(iul,'(/3A,F18.2,A)') ' wall clock time for   ',trim(token),' : ', &
              myMPIWallTime()-start,' s'
            call cpu_time(endCPU)
            write(iul,'(3A,F18.2,A//)') ' cpu time (master) for ',trim(token),' : ', &
            endCPU-startCPU,' s'
         end if
         call flush(iul)
      end if

   end do

   call cpu_time(tendCPU)
   sendbuf(1) = tendCPU - tstartCPU
   call myMPIReduceSumDouble(sendbuf,recvbuf,1)
   t = myMPIWallTime()-tstart

   if (MASTER .and. logmode >= 2) then
      write(iul,'(//a,a17)') ' wall clock time for run         : ',getTimeString(t)
      write(iul,'(a,f17.4)') ' total cpu time for run (core-h) : ',recvbuf(1)/3600.d0
      write(iul,'(a,f17.4)') ' cpu time per mpi process (h)    : ',recvbuf(1)/3600.d0/nproc
   end if

end subroutine mainloop

end module mainLoop_m

