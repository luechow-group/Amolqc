! Copyright (C) 2011-2013, 2015, 2018 Arne Luechow
! Copyright (C) 2014-2018 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optParamsELin_m

! optimize the energy with a modified linear method (Toulouse/Umrigar) for a fixed sample w.r.t.

   use kinds_m, only: r8
   use global_m
   use subloop_m, only: subloop
   use sorting_m, only: quickSortIndex
   use rwSample_m
   use elocAndPsiTermsLin_m
   use wfParameters_m
   use waveFunction_m, only: writeWF
   use ecp_m, only: EcpType
   use utils_m, only: intToStr,polyfit
   use mpiInterface_m, only: myMPIBcastDouble, myMPIBcastString
   implicit none

   private
   public :: eminlin_optimizeSample



contains


   subroutine eminlin_optimizeSample(lines,nl,WFP,sample,converged)
   !--------------------------------------------------------------!

   integer, intent(in)                  :: nl
   character(len=120), intent(in)       :: lines(nl)
   type(WFParamDef), pointer            :: WFP
   type(RWSample), intent(inout)        :: sample  ! (fixed) sample for optimization
   logical, intent(out)                 :: converged
   integer                              :: np, subSampleSize, npCI,npJ,npMO
   real(r8), allocatable                  :: p(:),p0(:)       ! parameter vector
   real(r8), allocatable                  :: delta_p(:,:) ! change of parameter vector (allow several)
   real(r8), allocatable                  :: H(:,:),S(:,:)  ! gradient and Hessian
   real(r8)                               :: e0, var, minVar, minE, lambda(6), lambdaOpt
   real(r8), allocatable                  :: fi(:),ELi(:),fiEL(:),fifj(:,:),fifjEL(:,:),fiELj(:,:)
   integer                              :: i,j,idxMin,maxDP, iflag, eqIter, eqStep, optIter
   integer                              :: nSize,root
   type(ElocAndPsiTermsLin)             :: EPsiTLin
   real(r8)                               :: targetE, targetVar,cffac,dmax,maxVar,lambdaMax
   character(len=80)                    :: subName, fname
   logical                              :: doWriteWF,subSample,manualEv,proJ,manualLambda,max_prj
   logical                              :: quad
   type(WFParamDerivTerms)              :: wfpDT
   real(r8)                               :: eRef,pe0,pvar
   type(EcpType)                        :: ecp
   logical                              :: exitSubloop = .false.

   converged = .true.
   manualEv = .false.
   call internal_readInput()         ! internal subroutine after 'contains'

   eRef = 0
   call ElocAndPsiTermsLin_create(EPsiTLin,eRef,WFP)

   if (logmode>=2) then
      write(iul,'(/A/)') '   - -  energy minimization using linear method: initialization  - -'
   endif

   np   = ElocAndPsiTermsLin_nParams(EPsiTLin)
   npCI = ElocAndPsiTermsLin_nPCI(EPsiTLin)
   npJ  = ElocAndPsiTermsLin_nPJ(EPsiTLin)
   npMO = ElocAndPsiTermsLin_nPMO(EPsiTLin)
   allocate(wfpDT%fi(np),wfpDT%fij(np,np),wfpDT%ELi(np))
   call assert(np>0,'eminlin_optimizeSample: no parameters')
   allocate(p(np),p0(np),delta_p(0:np,min(maxDP,np)),H(0:np,0:np),S(0:np,0:np))
   p = 0; p0 = 0; delta_p = 0; H = 0; S = 0
   allocate(fi(np),ELi(np),fiEL(np),fifj(np,np),fifjEL(np,np),fiELj(np,np))

   if (logmode >= 2) write(iul,*) ' starting wf parameter optimization with optType=',WFP%optType

   eqStep = 1
   do
      if(doWriteWF) call getPlusOptIter(optIter)
      call ElocAndPsiTermsLin_reset(EPsiTLin)
      call internal_calcEPsiTerms()
      e0 = ElocAndPsiTermsLin_EmeanALL(EPsiTLin)
      var = ElocAndPsiTermsLin_varALL(EPsiTLin)
      nSize = getSampleSizeAllNodes(sample)
      call ElocAndPsiTermsLin_resultALL(EPsiTLin,fi,ELi,fiEL,fifj,fifjEL,fiELj)
      if (logmode >= 2) then
         write(iul,'(a,f15.5,a,f12.5,a,f12.3,a,i10)') &
            ' with Emean=',e0,' +/- ',sqrt(var/nSize),' var=',var,' size=',nSize
         if (eqStep > 0) then
            write(iul,'(a,f15.5,a,f12.3)') ' Difference to projection: Delta E=',e0-pe0,' Delta var =',var-pvar
         end if
      end if

      if (var > maxVar) then
         converged = .false.
         exit
      !! else traceback to p=p0 and reduce e.g. trust radius
      end if


      if (MASTER) then
         p0 = wfparams_get(WFP)
         p = p0
         call internal_calcHandS()
         delta_p = rightEigenvector(H,S,np+1,1,maxDP)
      end if

       if (subSample) then
         if(manualEv) then
              call internal_setEV(root,minE,minVar,subSampleSize)
              if(.not. manualLambda) call internal_chooseStepLength(root,minE,minVar,lambdaOpt,subSampleSize)
              idxMin=root
              manualEv=.false.
         else
              call internal_chooseEV(maxDP,idxMin,minE,minVar,subSampleSize)
              if(.not. manualLambda) call internal_chooseStepLength(idxMin,minE,minVar,lambdaOpt,subSampleSize)
         endif
       else
         if(manualEv) then
              call internal_setEV(root,minE,minVar)
              if(.not. manualLambda) call internal_chooseStepLength(root,minE,minVar,lambdaOpt)
              idxMin=root
              manualEv=.false.
           else
              call internal_chooseEV(maxDP,idxMin,minE,minVar)
              if(.not. manualLambda) call internal_chooseStepLength(idxMin,minE,minVar,lambdaOpt)
         endif
       endif


      ! select best parameter set and update sample once again
      if (MASTER) then

         if (logmode >= 3) then
            write(iul,'(/a)') ' delta_p before symmetrisation:'
            write(iul,'(10F10.5)') delta_p(1:,idxMin)
         endif

         call wfparams_symmetriseMOs(WFP,delta_p(1:,idxMin))

         if (logmode >= 3) then
            write(iul,'(/a)') ' delta_p after symmetrisation:'
            write(iul,'(10F10.5)') delta_p(1:,idxMin)
         endif

         if(npCI>0) then
           p(1:npJ+npMO) = p0(1:npJ+npMO) + lambdaOpt*delta_p(1:npJ+npMO,idxMin)
           p(npJ+npMO+1:) = p0(npJ+npMO+1:) + delta_p(npJ+npMO+1:,idxMin)
        else
          p = p0 + lambdaOpt*delta_p(1:np,idxMin)
       endif

      end if
      call myMPIBcastDouble(p,np)
      call wfparams_set(WFP,p,.true.) ! with normalized CI coeefs
      call ElocAndPsiTermsLin_reset(EPsiTLin)
      call internal_reCalcSample()
      pe0 = ElocAndPsiTermsLin_EmeanALL(EPsiTLin)
      pvar = ElocAndPsiTermsLin_varALL(EPsiTLin)
      nSize = getSampleSizeAllNodes(sample)

      if (logmode >= 2) then
         write(iul,'(/a,i3,a,f10.2)')  ' final parameter set ev',idxMin,' lambda=',lambdaOpt
         write(iul,*) ' new parameter vector:'
         write(iul,'(10g12.4)') p
         write(iul,'(a,f15.5,a,f12.5,a,f12.3,a,i10)') &
            ' with projected Emean=',pe0,' +/- ',sqrt(var/nSize),' var=',pvar,' size=',nSize
             if (npCI>0) then
                        write(iul,*) ''
                        write(iul,*) 'ci coefficients are normalized'
                        write(iul,*) ''
              endif
      end if
      if (doWritewf) then
         fname = trim(baseName)//'-'//trim(intToStr(optIter))//'.wf'
         call writeWF(fname,.false.,ecp)
      end if
   if (eqStep >= eqIter) exit

      eqStep = eqStep + 1
      ! call "subroutine" subName in .in or macro subName.cmd that should contain
      ! code for equilibrating the sample with the new wave function
      call subloop(subName, sample,exitSubloop)
      if (exitSubloop) exit
   end do

   call setCurrentResult(pe0,0.d0,pvar)

   deallocate(p,delta_p,H,S)
   deallocate(fi,ELi,fiEL,fifj,fifjEL,fiELj)

   call ElocAndPsiTermsLin_destroy(EPsiTLin)

   contains


      subroutine internal_readInput()
      !-----------------------------!
         maxDP = 3
         call getinta(lines,nl,'max_ev=',maxDP,iflag)
         call getdbla(lines,nl,'target_E=',targetE,iflag)
         if (iflag /= 0) call abortp('eminlin: target_E option required')
         call getdbla(lines,nl,'target_var=',targetVar,iflag)
         if (iflag /= 0) call abortp('eminlin: target_var option required')
         maxVar = 1.d9
         call getdbla(lines,nl,'max_var=',maxVar,iflag)
         dmax = 1.d9
         call getdbla(lines,nl,'dmax=',dmax,iflag)
         cffac = 0.001d0
         call getdbla(lines,nl,'cffac=',cffac,iflag)
         doWriteWF = finda(lines,nl,'write_wf')
         subName = 'equilibrate'
         call getstra(lines,nl,'eq_call=',subName,iflag)
         eqIter = 0
         call getinta(lines,nl,'eq_iter=',eqIter,iflag)
         subSample=.false.
         call getinta(lines,nl,'ev_sample_size=',subSampleSize,iflag)
         if (iflag==0) subSample=.true.
         call getinta(lines,nl,'root=',root,iflag)
         if(iflag==0) manualEv=.true.
         proJ =.false.
         proJ = finda(lines,nl,'prt_prj')
         max_prj=.false.
         max_prj = finda(lines,nl,'max_prj')
         manualLambda=.false.
         call getdbla(lines,nl,'lambda=',lambdaOpt,iflag)
         if (iflag == 0) manualLambda=.true.
         lambdaMax=1.0d0
         call getdbla(lines,nl,'lambdamax=',lambdaMax,iflag)
         quad=.false.
         quad = finda(lines,nl,'quad')
      end subroutine internal_readInput


      subroutine internal_calcEPsiTerms(noDrv)
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for linear method
         logical,intent(in),optional :: noDrv
         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)

         if(present(noDrv)) then
          wfpDT%noCalc=.true.
            rwp => getFirst(sample)
               do
                  call pos(rwp,x,y,z)
                  call eConfigArray_set(ec,1,x,y,z)
                  call eloc(0, ec, 'none', wfpDef=WFP, wfpDT=wfpDT)
                  call ElocAndPsiTermsLin_add(EPsiTLin,wfpDT)
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
                  call ElocAndPsiTermsLin_add(EPsiTLin,wfpDT)
               if (.not.isNext(sample)) exit
                  rwp => getNext(sample)
               enddo
         endif
      end subroutine internal_calcEPsiTerms

      subroutine internal_SubcalcEPsiTerms(n,noDrv)
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for linear method
         ! only for first n samples
         integer,intent(in) ::  n
         logical,intent(in),optional :: noDrv
         integer           ::  i
         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)

         if (present(noDrv)) then
               i=0
               wfpDT%noCalc=.true.
               rwp => getFirst(sample)
               do
                  i=i+1
                  call pos(rwp,x,y,z)
                  call eConfigArray_set(ec,1,x,y,z)
                  call eloc(0, ec, 'none', wfpDef=WFP, wfpDT=wfpDT)
                  call ElocAndPsiTermsLin_add(EPsiTLin,wfpDT)
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
                  call ElocAndPsiTermsLin_add(EPsiTLin,wfpDT)
               if (.not.isNext(sample) .or. i>n) exit
                  rwp => getNext(sample)
               enddo
         endif
      end subroutine internal_SubcalcEPsiTerms

      subroutine internal_reCalcSample()
      !---------------------------------!

         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)

          wfpDT%noCalc=.true.
            rwp => getFirst(sample)
               do
                  call pos(rwp,x,y,z)
                  call eConfigArray_set(ec,1,x,y,z)
                  call eloc(0, ec, 'none', wfpDef=WFP, wfpDT=wfpDT)
                  call resetTo_without_Calc(rwp,x,y,z) ! reset rw
                  call ElocAndPsiTermsLin_add(EPsiTLin,wfpDT)! for ElocAndPsiTermsLin_EmeanALL
               if (.not.isNext(sample)) exit
                  rwp => getNext(sample)
               enddo
             wfpDT%noCalc=.false.
      end subroutine internal_reCalcSample


      subroutine internal_calcHandS()
      !-----------------------------!
         ! calculate H and S matrices of linear method
         H(0,0) = e0
         H(1:np,0) = fiEL(1:np)-fi(1:np)*e0
         H(0,1:np) = H(1:np,0) + ELi(1:np)
         do j=1,np
            H(1:np,j) = fifjEL(1:np,j) - fi(1:np)*fiEL(j) - fi(j)*fiEL(1:np) &
                      + fi(1:np)*fi(j)*e0 + fiELj(1:np,j) - fi(1:np)*ELi(j)
         enddo

         S = 0
         S(0,0) =1
         do j=1,np
            S(1:np,j) = fifj(1:np,j) - fi(1:np)*fi(j)
         enddo

         if (logmode >= 3) then
            write(iul,*) ' H:'
            do i=0,np
               write(iul,'(15G11.3)') H(i,:)
            enddo
            write(iul,*) ' S:'
            do i=0,np
               write(iul,'(15G11.3)') S(i,:)
            enddo
            write(iul,*) ' fi:'
            write(iul,'(15G11.3)') fi(:)
            write(iul,*) ' fiEL:'
            write(iul,'(15G11.3)') fiEL(:)
            write(iul,*) ' ELi:'
            write(iul,'(15G11.3)') ELi(:)
            write(iul,*) ' fifj:'
            do i=1,np
               write(iul,'(15G11.3)') fifj(i,:)
            enddo
            write(iul,*) ' fifjEL:'
            do i=1,np
               write(iul,'(15G11.3)') fifjEL(i,:)
            enddo
         endif
      end subroutine internal_calcHandS

      subroutine internal_chooseEV(n,idxMin,minE,minVar,subSampleSize)
      !------------------------------------------------!
         ! choose the correct eigenvalue (not necessary the lowest)
         integer :: n        ! test the first n eigenvalues
         integer :: idxMin
         real(r8)  :: minE, minVar
         integer,intent(in),optional :: subSampleSize
         integer i,idxPrj
         real(r8) e0,var,prj,maxPrj
         idxPrj=0
         maxPrj=0
         minVar = 1.d9
         idxMin = 0
         if (present(subSampleSize) .and. logmode >= 2) then
             write(iul,'(/a)') ' A subset of samples are used to idetify correct ev.'
             write(iul,'(/a,i5)') ' sub sample size =',subSampleSize
         endif
         if(.not. max_prj) then

               if (logmode >= 2) then
                 if (proJ ) then
                    write(iul,'(/a)') ' sample Emean and var and projection for ev:'
                 else
                    write(iul,'(/a)') ' sample Emean and var for lowest ev:'
                 endif
               endif
               do i=1,min(np,n)
                  if (MASTER) then
                     p = p0 + delta_p(1:np,i)
                  end if
                  call myMPIBcastDouble(p,np)
                  call wfparams_set(WFP,p)

                  call ElocAndPsiTermsLin_reset(EPsiTLin)
                  if (present(subSampleSize)) then
                     call internal_SubcalcEPsiTerms(subSampleSize,.true.)
                  else
                     call internal_calcEPsiTerms(.true.)
                  endif
                  e0 = ElocAndPsiTermsLin_EmeanALL(EPsiTLin)
                  var = ElocAndPsiTermsLin_varALL(EPsiTLin)
                  if (proJ) prj = internal_Projection(p0,p)
                  if (logmode >= 2) then
                     if(proJ) then
                        write(iul,'(i5,F15.5,F15.5,F15.5)') i,e0,var,prj
                     else
                        write(iul,'(i5,F15.5,F15.5)') i,e0,var
                     endif
                  endif
                  if (proJ) then
                     if (abs(prj)>abs(maxPrj)) then
                       maxPrj=prj
                       idxPrj=i
                     endif
                   endif
                  if (var < minVar) then
                     minVar = var
                     minE = e0
                     idxMin = i
                  end if
               end do
               if (proJ .and. MASTER) then
                   write(iul,*)
                   write(iul,*) "max projection= ",maxPrj," for ev= ",idxPrj
               endif
         else
               if (logmode >= 2) then
                 if (proJ ) then
                    write(iul,'(/a)') '  projection for ev:'
                 endif
               endif
               do i=1,min(np,n)
                  if (MASTER) then
                     p = p0 + delta_p(1:np,i)
                  end if
                  call myMPIBcastDouble(p,np)
                  call wfparams_set(WFP,p)
                  prj = internal_Projection(p0,p)
                     if (abs(prj)>abs(maxPrj)) then
                       maxPrj=prj
                       idxPrj=i
                     endif
                  if (logmode >= 2 .and. MASTER) then
                     if(proJ) then
                        write(iul,'(10x,i5,F15.5)') i,prj
                     endif
                  endif
               end do

               if (MASTER) then
                  p = p0 + delta_p(1:np,idxPrj)
               end if
               call myMPIBcastDouble(p,np)
               call wfparams_set(WFP,p)
               if (present(subSampleSize)) then
                  call internal_SubcalcEPsiTerms(subSampleSize,.true.)
               else
                  call internal_calcEPsiTerms(.true.)
               endif
               e0 = ElocAndPsiTermsLin_EmeanALL(EPsiTLin)
               var = ElocAndPsiTermsLin_varALL(EPsiTLin)
               minVar = var
               minE = e0
               idxMin = idxPrj
               if (logmode >= 2 .and. MASTER) then
                  write(iul,'(/a)') ' sample Emean and var for ev with laget projection:'
                  write(iul,'(F15.5,F15.5)') e0,var
                  write(iul,*)
                  write(iul,*) "max projection= ",maxPrj," for ev= ",idxPrj
               endif
         endif

      end subroutine internal_chooseEV

      subroutine internal_setEV(root,E,Var,subSampleSize)
      !------------------------------------------------!
         ! choose the correct eigenvalue (not necessary the lowest)
         integer,intent(in) :: root        ! desired  eigen value for the for the first iteration
         real(r8),intent(out)  :: E, Var
         integer,intent(in),optional :: subSampleSize



         if (logmode >= 2) write(iul,'(/a)') ' sample Emean and var for chosen ev:'


            if (MASTER) then
               p = p0 + delta_p(1:np,root)
            end if
            call myMPIBcastDouble(p,np)
            call wfparams_set(WFP,p)

            call ElocAndPsiTermsLin_reset(EPsiTLin)
            if (present(subSampleSize)) then
               call internal_SubcalcEPsiTerms(subSampleSize,.true.)
            else
               call internal_calcEPsiTerms(.true.)
            endif
            E = ElocAndPsiTermsLin_EmeanALL(EPsiTLin)
            Var = ElocAndPsiTermsLin_varALL(EPsiTLin)

            if (logmode >= 2) write(iul,'(i5,F15.5,F15.5)') root,E,var

      end subroutine internal_setEV


      subroutine internal_chooseStepLength(idx,e0,var0,lambdaOpt,subSampleSize)
      !------------------------------------------------------------------------
         ! choose the optimal step length based on a trust radius and a cost function
         integer :: idx       ! selected eigenvector
         real(r8)  :: e0,var0   ! E and Var for unit step length
         real(r8)  :: lambdaOpt ! calculated optimal step length
         integer, intent(in), optional :: subSampleSize
         integer n,cfIdx
         real(r8) d,cfmin,cf,var,cfa(4),a(4),lambdaa(4),cfs(6),x1,x2,x3,y1,y2,y3,r1,r2,r3


         if (NPCI >0 .and. npJ==0 .and. npMO==0) then
            lambdaOpt=1.d0
            if (logmode >= 2) write(iul,'(a)') ' lambda=1.0 for ci coeffs'
         elseif(quad .eqv. .false.) then
            lambda = (/ 0.02d0, 0.1d0, 0.3d0, 0.5d0, 0.7d0, 1.d0 /)

            d = sum(abs(delta_p(1:npJ+npMO,idx)))/(npJ+npMO)
            if (logmode >= 2) write(iul,'(a,i3,a,f15.5,a,f10.2,a,f15.5)')  &
               'best ev ',idx,' with E=',e0,' var=',var0,' d=',d

            if (logmode >= 2) write(iul,'(A,F15.5)') 'cffac=',cffac

            cfmin = abs(e0-targetE) + cffac* abs(var0-targetVar)
            cfIdx = 6
            cfs(cfIdx) = cfmin

            if (logmode >= 2) write(iul,'(A,F15.5)') &
                     ' lambda=1.0 with cf=',cfmin
               do n=1,size(lambda)-1

                  if (MASTER) p(1:npJ+npMO) = p0(1:npJ+npMO) + lambda(n)*delta_p(1:npJ+npMO,idx)
                  if (MASTER .and. npCI>0) p(npJ+npMO+1:)= p0(npJ+npMO+1:)+ delta_p(npJ+npMO+1:,idx)

                  call myMPIBcastDouble(p,np)
                  call wfparams_set(WFP,p)

                  call ElocAndPsiTermsLin_reset(EPsiTLin)
                  if (present(subSampleSize)) then
                     call internal_SubcalcEPsiTerms(subSampleSize,.true.)
                  else
                     call internal_calcEPsiTerms(.true.)
                  endif

                  e0 = ElocAndPsiTermsLin_EmeanALL(EPsiTLin)
                  var = ElocAndPsiTermsLin_varALL(EPsiTLin)

                  cf = abs(e0-targetE) + cffac* abs(var-targetVar)
                  cfs(n) = cf
                  if (logmode >= 2) write(iul,'(I5,A,F10.2,A,F15.5,A,F15.5,A,F15.5)') &
                     n,': lambda=',lambda(n),' Emean =',e0,' var = ',var,' cf = ',cf
                  if (cf < cfmin) then
                     cfmin = cf
                     cfIdx = n
                  end if

               end do

               if (lambda(cfIdx) < dmax/d) then
                  if (cfIdx > 1 .and. cfIdx < 6) then
                     x1 = lambda(cfIdx - 1)
                     x2 = lambda(cfIdx)
                     x3 = lambda(cfIdx + 1)
                     y1 = cfs(cfIdx - 1)
                     y2 = cfs(cfIdx)
                     y3 = cfs(cfIdx + 1)
                     r1 = x1 * (y3 - y2)
                     r2 = x2 * (y1 - y3)
                     r3 = x3 * (y2 - y1)
                     lambdaOpt = (x3 * r3 + x2 * r2 + x1 * r1) / (2 * (r1 + r2 + r3))
                  else
                     lambdaOpt = lambda(cfIdx)
                  end if
                  if (logmode >= 2) write(iul,'(a,f10.2)') ' choosing min cost function: lambda=',lambdaOpt
               else
                  lambdaOpt = dmax/d
                  if (logmode >= 2) write(iul,'(a,f10.2)') ' choosing trust radius: lambda=',lambdaOpt
               end if
             else
               lambdaa = (/ 0.02d0, 0.1d0, 0.6d0, 1.d0 /)
               if (logmode >= 2) write(iul,'(A,F15.5)') 'cffac=',cffac

               do n=1,4

                  if (MASTER) p(1:npJ+npMO) = p0(1:npJ+npMO) + lambdaa(n)*lambdaMax*delta_p(1:npJ+npMO,idx)
                  if (MASTER .and. npCI>0) p(npJ+npMO+1:)= p0(npJ+npMO+1:)+ delta_p(npJ+npMO+1:,idx)

                  call myMPIBcastDouble(p,np)
                  call wfparams_set(WFP,p)

                  call ElocAndPsiTermsLin_reset(EPsiTLin)
                  if (present(subSampleSize)) then
                     call internal_SubcalcEPsiTerms(subSampleSize,.true.)
                  else
                     call internal_calcEPsiTerms(.true.)
                  endif

                  e0 = ElocAndPsiTermsLin_EmeanALL(EPsiTLin)
                  var = ElocAndPsiTermsLin_varALL(EPsiTLin)

                  cfa(n) = abs(e0-targetE) + cffac* abs(var-targetVar)
                  if (logmode >= 2) write(iul,'(I5,A,F10.2,A,F15.5,A,F15.5,A,F15.5)') &
                     n,': lambda=',lambdaMax*lambdaa(n),' Emean =',e0,' var = ',var,' cf = ',cfa(n)


               end do
               a=polyfit(lambdaa*lambdaMax, cfa, 2)
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


      real(r8) function internal_Projection(p0,pNew) !projects new CI vector to initial one
      real(r8), intent(in) :: p0(npCI),pNew(npCI)
      real(r8)             :: ci(npCI+1),cif(npCI+1)

      if (npCI<1 .or. npJ>0 .or. npMO>0 ) call abortp('eminlin: internal_printProjection is only valid for ci params')
       ci=cicoeffs_get(WFP)
       ci(2:)=p0
       cif(1)=ci(1)
       cif(2:)=pNew
       ci  = ci/sqrt(dot_product(ci,ci))
       cif = cif/sqrt(dot_product(cif,cif))
       internal_Projection=dot_product(ci,cif)

       end function internal_Projection




   end subroutine eminlin_optimizeSample




   function rightEigenvector(H,S,n,debugLevel,maxDP)
   !  ----------------------------------------------

   ! get the right eigenvector to the lowest (almost) real right eigenvalue
   real(r8),intent(inout):: H(:,:)
   real(r8),intent(inout):: S(:,:)
   integer,intent(in)  :: n        ! dimension of H and S
   integer,intent(in)  :: debugLevel
   integer,intent(in)  :: maxDP    ! how many eigenvalues (sorted) to return
   real(r8)              :: rightEigenvector(size(H,1),min(n-1,maxDP))

   integer            :: lwork
   real(r8)             :: alphar(n),alphai(n),beta(n),work(1,8*n),vr(n,n),vl(n,n)
   real(r8)             :: ev(n),evi(n),mav
   integer            :: info
   integer            :: i,idx(n)

   call assert(size(H,1) == n,'rightEigenvector: wrong dimension 1')

   lwork = 8*n

   alphar = 0d0; alphai = 0d0; beta =0d0; ev = 0d0; evi=0d0
   work = 0d0; vr = 0d0; vl = 0d0

   ! Lapack general nonsymmetric eigenvalue problem
   call DGGEV('N','V',n,H,n,S,n,alphar,alphai,beta,vl,n,vr,n,work,lwork,info)

   if (info < 0) then
      write(iul,*) 'the', info, '-th argument had an illegal value'
      call abortp('DGGEV failed')
   else if (info > 0) then
      write(iul,*) 'QZ iteration failed. Try it with one more iteration or a better sample'
      call abortp('DGGEV failed')
   end if

   do i=1,n
      ev(i) = alphar(i)/beta(i)
      evi(i) = alphai(i)/beta(i)
      vr(:,i) = vr(:,i) / vr(1,i)    ! normalize to 1st element==1
      idx(i) = i
   end do

   call quickSortIndex(isGreaterEqual,isSmallerEqual,idx)

   if (debugLevel > 0 .and. mytid==0) then
      write(901,'(/A)') 'DGGEV eigenvalues: Re, Im, beta'
      do i=1,min(n,10)
         mav = sum(abs(vr(2:n,idx(i)))) / (n-1)
         write(901,'(i5,4g18.10)') i,ev(idx(i)),evi(idx(i)),beta(idx(i)),mav
      enddo
   endif

   do i=1,min(n-1,maxDP)
      rightEigenvector(:,i) = vr(:,idx(i))
   end do

   contains

      logical function isGreaterEqual(i,j)
         integer, intent(in) :: i,j
         isGreaterEqual = ev(i) >= ev(j)
      end function isGreaterEqual

      logical function isSmallerEqual(i,j)
         integer, intent(in) :: i,j
         isSmallerEqual = ev(i) <= ev(j)
      end function isSmallerEqual


   end function rightEigenvector

end module optParamsELin_m
