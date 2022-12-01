! Copyright (C) 2016, 2018 Arne Luechow
! Copyright (C) 2016-2018 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optParamsPOpt_m

! optimize the energy with a modified Genear method (Toulouse/Umrigar) for a fixed sample w.r.t.

   use kinds_m, only: r8, i8
   use globalUtils_m, abortprog => abortp
   use global_m
   use subloop_m, only: subloop
   use sorting_m, only: quickSortIndex
   use rwSample_m
   use elocAndPsiTermsGen_m
   use wfParameters_m
   use waveFunction_m, only: writeWF
   use ecp_m, only: EcpType
   use utils_m, only: intToStr, polyfit
   use random_m
   use mpiInterface_m, only: myMPIBcastDouble, myMPIBcastString, myMPIWallTime
   implicit none

   private
   public :: eminpopt_optimizeSample



contains


   subroutine eminpopt_optimizeSample(lines,nl,WFP,sample,converged)
   !---------------------------------------------------------------!

   integer, intent(in)                  :: nl
   character(len=120), intent(in)       :: lines(nl)
   type(WFParamDef), pointer            :: WFP
   type(RWSample), intent(inout)        :: sample  ! (fixed) sample for optimization
   logical, intent(out)                 :: converged
   integer                              :: np,npCI,npMO,npJ
   integer(i8)                            :: ii, deltaESampleSize,stpSampleSize
   real(r8), allocatable                  :: p(:),p0(:)       ! parameter vector
   real(r8), allocatable                  :: delta_p(:)       ! change of parameter vector
   real(r8), allocatable                  :: delta_e(:)       ! perturbation theory denominator
   real(r8), allocatable                  :: g(:) ,S(:,:)     ! gradient and overlap matrix
   real(r8) e0, var, lambdaOpt
   real(r8), allocatable                  :: fi(:),ELi(:),fiEL(:),fifj(:,:),fifjEL(:,:),fiELj(:,:)
   real(r8), allocatable                  :: bb(:),Imat(:,:)
   integer lwork,i,j,info,iflag, eqIter, eqStep, nSize,io,optIter
   integer, allocatable                 :: ipiv(:)
   real(r8), allocatable                  :: work(:)
   type(ElocAndPsiTermsGen)             :: EPsiTGen
   real(r8)                               :: targetE, targetVar,cffac,dmax,maxVar,lambdaMax
   character(len=80)                    :: subName, fname, deltaEFileName
   logical                              :: doWriteWF, doWriteDE, doReadDE, subSample, deltaEOnce, doDeltaE, found
   logical                              :: stbl,stpl
   real(r8)                               :: eRef, pe0, pvar, allDeltaE, nu
   real(r8)                               :: x(ne),y(ne),z(ne)
   real(r8)                               :: tstart, t
   type(RandomWalker), pointer          :: rwp
   type(WFParamDerivTerms)              :: wfpDT
   type(eConfigArray)                   :: ec
   type(EcpType)                        :: ecp
   logical                              :: exitSubloop = .false.


   converged = .true.

   call internal_readInput()         ! internal subroutine after 'contains'

   eRef = 0
   call ElocAndPsiTermsGen_create(EPsiTGen,eRef,WFP)

   if (logmode>=2) then
      write(iul,'(/A/)') '   - -  energy minimization using perturbative method: initialization  - -'
   endif

   np    =  ElocAndPsiTermsGen_nParams(EPsiTGen)
   npCI  =  ElocAndPsiTermsGen_nPCI(EPsiTGen)
   npJ   =  ElocAndPsiTermsGen_npJ(EPsiTGen)
   npMO  =  ElocAndPsiTermsGen_npmo(EPsiTGen)
   call assert(np>0,'eminpopt_optimizeSample: no parameters')
   allocate(p(np), p0(np), delta_p(np), delta_e(np), S(np,np), g(np), bb(np))
   allocate(ipiv(np),work(np*np))
   if(stbl) allocate(Imat(np,np))
   p = 0; p0 = 0; delta_p = 0; S = 0; g = 0
   allocate(fi(np),ELi(np),fiEL(np),fifj(np,np),fifjEL(np,np),fiELj(np,np))
   lwork = np*np

   if (logmode >= 2) then
      write(iul,*) ' starting wf parameter optimization with optType=',WFP%optType
      write(iul,'(a,g12.3,a,g12.3)') ' max_var =', maxVar
   end if

   call eConfigArray_new(ec,ne,1)
   allocate(wfpDT%fi(np), wfpDT%Eli(np), wfpDT%fij(np,np))

   !read delta E from deltaE.dat
   if(doReadDE) call internal_readdE()
   ! use constant Delta E instead of calculated value
   if (allDeltaE > 1.d-9) then
      delta_e = allDeltaE
      if (logmode >= 2) write(iul,'(/a,g10.3/)') ' using constant delta E throughout optimisation:', allDeltaE
   end if

   eqStep = 1
   do
      if(doWriteWF) call getPlusOptIter(optIter)
      tstart = myMPIWallTime()
      call ElocAndPsiTermsGen_reset(EPsiTGen)

      doDeltaE = (abs(allDeltaE) < 1.d-9) .and. ( ( deltaEOnce .and. (eqStep==1) ) .or. .not. deltaEOnce )
      if(doReadDE) doDeltaE = .false.
      if (logmode >= 3) write(iul,*) ' calculate delta E this iteration:',doDeltaE

      wfpDT%fiCalc = .true.
      wfpDT%ELiCalc = doDeltaE
      wfpDT%fijCalc = .false.

      rwp => getFirst(sample)
      ii = 1
      do
         call pos(rwp,x,y,z)
         call eConfigArray_set(ec,1,x,y,z)
         call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
         !!!!call resetTo(rwp,x,y,z,WFP%optType)  ! recalculate and reset walker
         call ElocAndPsiTermsGen_add(EPsiTGen,wfpDT)
         if (doDeltaE .and. ii == deltaESampleSize) then
            ! calculate delta E with smaller sample
            e0 = ElocAndPsiTermsGen_EmeanALL(EPsiTGen)
            call ElocAndPsiTermsGen_resultALL(EPsiTGen,fi=fi,ELi=ELi,fiEL=fiEL,fifj=fifj,   &
                                              fifjEL=fifjEL,fiELj=fiELj)
            if (MASTER) then
               do j=1,np
                  S(j,j) = fifj(j,j) - fi(j)*fi(j)
                  delta_e(j) = ( fifjEL(j,j) - fi(j)*fiEL(j) - fi(j)*fiEL(j) &
                             + fi(j)*fi(j)*e0 + fiELj(j,j) - fi(j)*ELi(j) ) / S(j,j) - e0
               end do
               if (logmode >= 3) then
                  write(iul,*) ' delta E calculated with sample size:',ii
                  write(iul,'(10G11.3)') delta_e(:)
               endif
            end if
            wfpDT%ELiCalc = .false.
         end if
      if (.not.isNext(sample)) exit
         rwp => getNext(sample)
         ii = ii + 1
      end do

      e0 = ElocAndPsiTermsGen_EmeanALL(EPsiTGen)
      var = ElocAndPsiTermsGen_varALL(EPsiTGen)
      nSize = getSampleSizeAllNodes(sample)
      if (logmode >= 2) then
         write(iul,'(a,f15.5,a,f12.5,a,f12.3,a,i10)') &
            ' with Emean=',e0,' +/- ',sqrt(var/nSize),' var=',var,' size=',nSize
         if (eqStep > 0) then
            write(iul,'(a,f15.5,a,f12.3)') ' Difference to projection: Delta E=',e0-pe0,' Delta var =',var-pvar
         end if
      end if

      if (doDeltaE .and. deltaESampleSize == 0) then
         call ElocAndPsiTermsGen_resultALL(EPsiTGen,fi=fi,ELi=ELi,fiEL=fiEL,fifj=fifj,   &
                                           fifjEL=fifjEL,fiELj=fiELj)

         if (MASTER) then
            do j=1,np
               S(j,j) = fifj(j,j) - fi(j)*fi(j)
               delta_e(j) = ( fifjEL(j,j) - fi(j)*fiEL(j) - fi(j)*fiEL(j) &
                          + fi(j)*fi(j)*e0 + fiELj(j,j) - fi(j)*ELi(j) ) / S(j,j) - e0
            end do
            if (logmode >= 3) then
               write(iul,*) ' delta_e (full sample):'
               write(iul,'(10G11.3)') delta_e(:)
            end if
         end if
      else
         call ElocAndPsiTermsGen_resultALL(EPsiTGen,fi=fi,fiEL=fiEL,fifj=fifj)
      end if

      if (var > maxVar) then
         converged = .false.
         exit
      !! else traceback to p=p0 and reduce e.g. trust radius
      end if


      if (MASTER) then
         p0 = wfparams_get(WFP)
         do j=1,np
            S(1:np,j) = fifj(1:np,j) - fi(1:np)*fi(j)
         enddo
         g(1:np) = 2*(fiEL(1:np)-fi(1:np)*e0)

         bb = -g
         if(stbl) then
            call dpotrf('L',np,S,np,info)
            if (info == 0) then
               write(iul,*) ' -> Aprx Hessian positive definite!'
            else if (info > 0) then
               write(iul,*) ' -> Aprx Hessian is not positive definite!'
               write(iul,*) ' -> adding nu to the diagonal elements'

               S=S+nu*Imat
               do i=1,10
                  call dpotrf('L',np,S,np,info)
                  if (info == 0) then
                     write(iul,*) ' -> Aprx Hessian positive definite!'
                  else if (info > 0) then
                     write(iul,'(i3,a,f15.6,a)') i,':  nu = ',nu, ' aprx Hessian not positive definite'
                     nu=nu*4
                     S=S+nu*Imat
                  else
                     write(iul,*) ' !!! WARNING !!! Cholesky decomp error'
                  end if
               end do
            else
               write(iul,*) ' !!! WARNING !!! Cholesky decomp error'
            end if
         end if

         ! calculate bb = S^-1 * (-g) by solving S*bb = -g
         call DSYSV('L',np,1,S,np,ipiv,bb,np,work,lwork,info)

         if (logmode >= 3) write(iul,*) ' info = ',info

         delta_p = bb / (2*delta_e)
      end if

      ! select best parameter set and update sample once again
      if (stpl .eqv. .false.) then
         lambdaOpt = 1.d0
      else
         if (stpSampleSize ==0 ) then
           call internal_chooseStepLength(lambdaOpt)
         else
           call internal_chooseStepLength(lambdaOpt,stpSampleSize)
        endif
      endif

      if (MASTER) then
         p = p0 + lambdaOpt*delta_p(1:np)
      end if
      call myMPIBcastDouble(p,np)
      call wfparams_set(WFP,p,.true.) ! with normalized CI coeefs
      call ElocAndPsiTermsGen_reset(EPsiTGen)

      wfpDT%fiCalc = .false.
      wfpDT%ELiCalc = .false.
      wfpDT%fijCalc = .false.

      rwp => getFirst(sample)
      do
         call pos(rwp,x,y,z)
         call eConfigArray_set(ec,1,x,y,z)
         call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
         call resetTo_without_Calc(rwp,x,y,z) ! reset rw
         !!!!call resetTo(rwp,x,y,z,WFP%optType)  ! recalculate and reset walker
         call ElocAndPsiTermsGen_add(EPsiTGen,wfpDT)
      if (.not.isNext(sample)) exit
         rwp => getNext(sample)
      enddo

      pe0 = ElocAndPsiTermsGen_EmeanALL(EPsiTGen)
      pvar = ElocAndPsiTermsGen_varALL(EPsiTGen)
      nSize = getSampleSizeAllNodes(sample)

      if (logmode >= 2) then
         write(iul,*) ' new parameter vector:'
         write(iul,'(10g12.4)') p
         write(iul,'(a,f15.5,a,f12.5,a,f12.3,a,i10)') &
            ' with projected Emean=',pe0,' +/- ',sqrt(var/nSize),' var=',pvar,' size=',nSize
      end if
      if (doWritewf) then
         fname = trim(baseName)//'-'//trim(intToStr(optIter))//'.wf'

         call writeWF(fname,.false.,ecp)
      end if

      t = myMPIWallTime()-tstart
      if (MASTER .and. logmode >= 2) then
         write(iul,'(/a,f18.2,a)') ' wall clock time for optimisation step : ',t,' s'
      end if

     if (doWriteDE .and. MASTER) then
      if(.not.doWriteWF) call getPlusOptIter(optIter)
      fname = trim(baseName)//'-'//trim(intToStr(optIter))//'.dat'
        open(21,file=fname,iostat=io)
        if (io /= 0) call abortprog('could not open deltaE.dat')
        do i=1,size(delta_e)
          write(21,*) delta_e(i)
         enddo
         close(21)
     endif

   if (eqStep >= eqIter) exit

      eqStep = eqStep + 1
      ! call "subroutine" subName in .in or macro subName.cmd that should contain
      ! code for equilibrating the sample with the new wave function
      call subloop(subName,sample,exitSubloop)
      if (exitSubloop) exit
   end do

   call setCurrentResult(pe0,0.d0,pvar)

   call ElocAndPsiTermsGen_destroy(EPsiTGen)

   contains


      subroutine internal_readInput()
      !-----------------------------!
         stbl=.false.
         stpl=.false.
         lambdaMax=1.0d0
         call getdbla(lines,nl,'lambdamax=',lambdaMax,iflag)
         cffac = 0.001
         call getdbla(lines,nl,'cffac=',cffac,iflag)
         found = finda(lines,nl,'quad')
         if (found) stpl=.true.
         call getdbla(lines,nl,'target_E=',targetE,iflag)
         call getdbla(lines,nl,'target_var=',targetVar,iflag)
         if (stpl .eqv. .true. .and. iflag /= 0) call abortprog('pOpt: target_E and target_var are required if quad is specified')
         nu = 0.01d0
         call getdbla(lines,nl,'nu=',nu,iflag)
         if (iflag ==0) stbl=.true.
         doWriteWF = finda(lines,nl,'write_wf')
         doWriteDE = finda(lines,nl,'write_de')
         doReadDE=.false.
         call getstra(lines,nl,'delta_e_filename=',deltaEFileName,iflag)
         if(iflag==0) doReadDE=.true.

         maxVar = 1.d9
         call getdbla(lines,nl,'max_var=',maxVar,iflag)
         dmax = 1.d9
         allDeltaE = 0.d0
         call getdbla(lines,nl,'delta_e=',allDeltaE,iflag)
         subName = 'equilibrate'
         call getstra(lines,nl,'eq_call=',subName,iflag)
         eqIter = 0
         call getinta(lines,nl,'eq_iter=',eqIter,iflag)
         stpSampleSize=0
         call getint8a(lines,nl,'stp_sample_size=',stpSampleSize,iflag)
         deltaESampleSize = 0
         call getint8a(lines,nl,'delta_e_sample_size=',deltaESampleSize,iflag)

         deltaEOnce = .false.
         found = finda(lines,nl,'delta_e_once')
         if (found) deltaEOnce = .true.
         found = finda(lines,nl,'delta_e_always')
         if (found) deltaEOnce = .false.
         if (iflag==0) subSample=.true.
      end subroutine internal_readInput


      subroutine internal_calcEPsiTerms(noDrv)
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for Genear method
         logical,intent(in),optional :: noDrv
         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(WFParamDerivTerms) :: wfpDT
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)
       if(present(noDrv)) then
          wfpDT%fiCalc = .false.
          wfpDT%ELiCalc = .false.
          wfpDT%fijCalc = .false.
          wfpDT%noCalc=.true.
          rwp => getFirst(sample)
           do
            call pos(rwp,x,y,z)
            call eConfigArray_set(ec,1,x,y,z)
            call eloc(0, ec, 'none', wfpDef=WFP, wfpDT=wfpDT)
            !!!!call resetTo(rwp,x,y,z,WFP%optType)  ! recalculate and reset walker
            call ElocAndPsiTermsGen_add(EPsiTGen,wfpDT)
           if (.not.isNext(sample)) exit
              rwp => getNext(sample)
           enddo
           wfpDT%noCalc=.false.
       else
           rwp => getFirst(sample)
           do
            call pos(rwp,x,y,z)
            call eConfigArray_set(ec,1,x,y,z)
            call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
            !!!!call resetTo(rwp,x,y,z,WFP%optType)  ! recalculate and reset walker
            call ElocAndPsiTermsGen_add(EPsiTGen,wfpDT)
           if (.not.isNext(sample)) exit
              rwp => getNext(sample)
           enddo
      endif
      end subroutine internal_calcEPsiTerms

      subroutine internal_SubcalcEPsiTerms(n,noDrv)
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for linear method
         ! only for first n samples
         integer(i8),intent(in) ::  n
         logical,intent(in),optional :: noDrv
         integer(i8)           ::  i
         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(WFParamDerivTerms) :: wfpDT
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)

         if (present(noDrv)) then
               i=0
         wfpDT%fiCalc = .false.
         wfpDT%ELiCalc = .false.
         wfpDT%fijCalc = .false.
               wfpDT%noCalc=.true.
               rwp => getFirst(sample)
               do
                  i=i+1
                  call pos(rwp,x,y,z)
                  call eConfigArray_set(ec,1,x,y,z)
                  call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
                  call ElocAndPsiTermsGen_add(EPsiTGen,wfpDT)
               if (.not.isNext(sample) .or. i>n) exit
                  rwp => getNext(sample)
               enddo

            wfpDT%noCalc=.false.
         else
               i=0
               rwp => getFirst(sample)
               do
                  i=i+1
                  call pos(rwp,x,y,z)
                  call eConfigArray_set(ec,1,x,y,z)
                  call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
                  call ElocAndPsiTermsGen_add(EPsiTGen,wfpDT)
               if (.not.isNext(sample) .or. i>n) exit
                  rwp => getNext(sample)
               enddo
         endif
      end subroutine internal_SubcalcEPsiTerms
      subroutine internal_readdE()
      !---------------------------------!
            if(MASTER) then
               open(21,file=deltaEFileName,iostat=io)
                 if (io /= 0) then
                      call abortprog(' could not open delta E file')
                 else
                        write(iul,*)""
                        write(iul,*)"  reading delta E from ",trim(deltaEFileName)," file ... "
                        write(iul,*)""

                       do i=1,size(delta_e)
                         read(21,*) delta_e(i)
                       enddo
                 endif
               close(21)
            endif
         call myMPIBcastDouble(delta_e,size(delta_e))

      end subroutine internal_readdE

      subroutine internal_chooseStepLength(lambdaOpt,stpSampleSize)
      !------------------------------------------------------------------------
         ! choose the optimal step length based on a trust radius and a cost function

         real(r8),intent(out)           :: lambdaOpt ! calculated optimal step length
         integer(i8),intent(in), optional :: stpSampleSize! calculated optimal step length

         integer n
         real(r8) cf(4),var,lambda(4),a(3)

         if (NPCI >0 .and. npJ==0 .and. npMO==0) then
            lambdaOpt=1.d0
            if (logmode >= 2) write(iul,'(a)') ' lambda=1.0 for ci coeffs'
         else
            lambda = (/ 0.02d0, 0.1d0, 0.6d0, 1.d0 /)


               do n=1,4

                  if (MASTER) p(1:npJ+npMO) = p0(1:npJ+npMO) + lambda(n)*lambdaMax*delta_p(1:npJ+npMO)
                  if (MASTER .and. npCI>0) p(npJ+npMO+1:)= p0(npJ+npMO+1:)+ delta_p(npJ+npMO+1:)

                  call myMPIBcastDouble(p,np)
                  call wfparams_set(WFP,p)

                  call ElocAndPsiTermsGen_reset(EPsiTGen)
                  if (present(stpSampleSize)) then
                     call internal_SubcalcEPsiTerms(stpSampleSize,.true.)
                  else
                     call internal_calcEPsiTerms(.true.)
                  endif

                  e0 = ElocAndPsiTermsGen_EmeanALL(EPsiTGen)
                  var = ElocAndPsiTermsGen_varALL(EPsiTGen)

                  cf(n) = abs(e0-targetE) + cffac* abs(var-targetVar)
                  if (logmode >= 2) write(iul,'(I5,A,F10.2,A,F15.5,A,F15.5,A,F15.5)') &
                     n,': lambda=',lambda(n)*lambdamax,' Emean =',e0,' var = ',var,' cf = ',cf(n)

               end do
               a=polyfit(lambda*lambdaMax, cf, 2)
               if (a(3)<0) then
                  write(iul,*) '  a=', a(3)
                  write(iul,*) ' WARNING: there is no min for quadratic model!!!'
                  lambdaOpt=lambdaMax
                  write(iul,*) ' lambda=', lambdaOpt
               else
                lambdaOpt=-a(2)/(2.0*a(3))
                write(iul,*) ' lambda=', lambdaOpt
               endif
        endif


      end subroutine internal_chooseStepLength

   end subroutine eminpopt_optimizeSample





end module optParamsPOpt_m
