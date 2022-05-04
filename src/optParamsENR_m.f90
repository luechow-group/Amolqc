! Copyright (C) 2011-2016, 2018 Arne Luechow
! Copyright (C) 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module optParamsENR_m

! optimize the energy for a fixed sample using the Newton-Raphson alg.

   use kinds_m, only: r8
   use global_m
   use subloop_m, only: subloop
   use rWSample_m
   use elocAndPsiTermsENR_m
   use wfParameters_m
   use waveFunction_m, only: writeWF
   use ecp_m, only: EcpType
   use utils_m, only: intToStr
   use mpiInterface_m, only: myMPIBcastDouble, myMPIBcastString
   implicit none

   private
   public :: eminNR_optimizeSample


contains



   subroutine eminNR_optimizeSample(lines,nl,WFP,sample,converged)
   !-------------------------------------------------------------!

   integer, intent(in)                  :: nl
   character(len=120), intent(in)       :: lines(nl)
   type(WFParamDef), pointer            :: WFP              !
   type(RWSample), intent(inout)        :: sample           ! (fixed) sample for optimization
   logical, intent(out)                 :: converged
   integer                              :: nParams, np,npCI
   real(r8), allocatable                  :: p(:),p0(:)       ! parameter vector
   real(r8), allocatable                  :: delta_p(:)       ! change of parameter vector
   real(r8), allocatable                  :: g(:),bb(:),H(:,:),H0(:,:)  ! gradient and Hessian
   real(r8)                               :: e0,var,pe0,pvar
   real(r8), allocatable                  :: fi(:),ELi(:),fiEL(:)
   real(r8), allocatable                  :: fifj(:,:),fifjEL(:,:),fiELj(:,:),fij(:,:),fijEL(:,:)
   real(r8), allocatable                  :: A(:,:),B(:,:),D(:,:),Imat(:,:)
   !!!real(r8), allocatable                  :: eval(:),evec(:,:)
   integer lwork,i,j,info,n,ierr,eqIter, eqStep, nSize, NRMode, iflag, iter,optIter
   integer, allocatable                 :: ipiv(:)
   real(r8), allocatable                  :: work(:)
   real(r8)                               :: targetE, targetVar,cffac,dmax,maxVar, lambda(6), lambdaOpt, eRef, gfac
   real(r8)                               :: nu,r,delta_q,delta_f,normdp,mabsdp,nuStart,deltaFmin
   type(ElocAndPsiTermsENR)             :: EPsiTENR
   character(len=80)                    :: subName,fname
   logical                              :: doWriteWF, keepnu
   type(WFParamDerivTerms)              :: wfpDT
   type(EcpType)                        :: ecp
   logical                              :: exitSubloop = .false.

   converged = .true.
   keepnu = .false.



   call internal_readInput()        ! internal subroutine after contains



   eRef = 0.d0
   call ElocAndPsiTermsENR_create(EPsiTENR,eRef,WFP)

   if (logmode>=2) then
      write(iul,'(/A/)') '   - -  energy minimization using Newton-Raphson: initialization  - -'
      write(iul,'(a,i3,a,f12.3)') ' parameters:  nrmethod = ',NRMode,'   gradient factor = ',gfac
   endif

   np = ElocAndPsiTermsENR_nParams(EPsiTENR)
   npCI = ElocAndPsiTermsENR_nPCI(EPsiTENR)
   WFP => ElocAndPsiTermsENR_getWFP(EPsiTENR)


   call assert(np>0,'eminNR_optimizeSample: no parameters')
   allocate(p(np),p0(np),delta_p(np),bb(np),g(np),H(np,np),H0(np,np),ipiv(np),work(np*np))
   allocate(fi(np),ELi(np),fiEL(np))
   allocate(fifj(np,np),fifjEL(np,np),fiELj(np,np),fij(np,np),fijEL(np,np))
   allocate(A(np,np),B(np,np),D(np,np),Imat(np,np))
   !!!allocate(eval(np),evec(np,np))

   allocate(wfpDT%fi(np),wfpDT%fij(np,np),wfpDT%ELi(np))

   Imat = 0;
   do i=1,np
      Imat(i,i) = 1.d0
   end do

   nu = nuStart

   p = 0; delta_p = 0

   if (logmode >= 2) write(iul,*) ' starting wf parameter optimization with optType=',WFP%optType

   eqStep = 1
   do
      if(doWriteWF) call getPlusOptIter(optIter)
      call ElocAndPsiTermsENR_reset(EPsiTENR)
      call internal_calcEPsiTerms()

      e0 = ElocAndPsiTermsENR_EmeanALL(EPsiTENR)
      var = ElocAndPsiTermsENR_varALL(EPsiTENR)
      nSize = getSampleSizeAllNodes(sample)
      call ElocAndPsiTermsENR_resultALL(EPsiTENR,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)
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

         call internal_calcGradAndHessian(g,H)
         H0 = H

         if (logmode >=2) then
            call dpotrf('L',np,H,np,info)
            if (info == 0) then
               write(iul,*) ' -> Hessian positive definite!'
            else if (info > 0) then
               write(iul,*) ' -> Hessian not positive definite!'
            else
               write(iul,*) ' !!! WARNING !!! Cholesky decomp error'
            end if
         end if

         if (NRMode == 1 .or. NRMode == 2) then

            H = H0

         else if (NRMode == 3) then

            H = H0 + nu*Imat

            if (logmode >= 2) write(iul,'(/a)') ' find Newton step:'
            do iter=1,10
               call dpotrf('L',np,H,np,info)  ! Cholesky decomposition
               if (info == 0) then
                  write(iul,'(i3,a,f15.6,a)') iter,':  nu = ',nu, ' Hessian positive definite'
               else if (info > 0) then
                  write(iul,'(i3,a,f15.6,a)') iter,':  nu = ',nu, ' Hessian not positive definite'
               else
                  write(iul,*) ' !!! WARNING !!! Cholesky decomp error'
               end if
               if (info > 0) then !! not positive definite
                  nu = 4*nu
                  H = H0 + nu*Imat
               else
                  exit
               end if
            end do
            H = H0 + nu*Imat

         else if (NRMode == 4) then

            H = H0
            do i=1,np
               H(i,i) = H(i,i)*(1+nu)
            enddo

            if (logmode >= 2) write(iul,'(/a)') ' find Newton step:'
            do iter=1,10
               call dpotrf('L',np,H,np,info)  ! Cholesky decomposition
               if (info == 0) then
                  write(iul,'(i3,a,f15.6,a)') iter,':  nu = ',nu, ' Hessian positive definite'
               else if (info > 0) then
                  write(iul,'(i3,a,f15.6,a)') iter,':  nu = ',nu, ' Hessian not positive definite'
               else
                  write(iul,*) ' !!! WARNING !!! Cholesky decomp error'
               end if
               if (info > 0) then !! not positive definite
                  nu = 4*nu
                  H = H0
                  do i=1,np
                     H(i,i) = H(i,i)*(1+nu)
                  enddo
               else
                  exit
               end if
            end do
            H = H0
            do i=1,np
               H(i,i) = H(i,i)*(1+nu)
            enddo
         end if

         p = wfparams_get(WFP)
         p0 = p

         lwork = np*np
         bb = -g

         call internal_writeDataToLogFile()

         ! calculate bb = -H^-1 * g by solving H*bb = -g
         call DSYSV('L',np,1,H,np,ipiv,bb,np,work,lwork,info)

         if (info/=0 .and. logmode>=2) write(iul,*) ' !!! DSYSV failed: INFO=',info
         if (logmode >= 3) then
            write(iul,'(/a)') ' DSYSV result b:'
            write(iul,'(10F10.5)') bb
         endif

         delta_p = bb                                    ! change of parameter vector
         call wfparams_symmetriseMOs(WFP,delta_p)

         if (logmode >= 3) then
            write(iul,'(/a)') ' delta_p after symmetrisation:'
            write(iul,'(10F10.5)') delta_p
         endif

         normdp = sqrt(dot_product(delta_p,delta_p))     ! 2-norm
         mabsdp = sum(abs(delta_p))/np                   ! mean abs of vector components

         ! calculate energy change delta_q for the quadratic model
         delta_q = 0.5d0*dot_product(delta_p,matmul(H0,delta_p)) + dot_product(g,delta_p)

         if (logmode >= 2) then
            write(iul,*)
            write(iul,*) ' norm(delta_p)=',normdp
            write(iul,*) ' mean abs(delta_p_i)=',mabsdp
            write(iul,*) ' delta_q =',delta_q
         end if

         call internal_writeNewVectorToLogFile()
      endif

      lambdaOpt = 1.d0
      if (NRMode == 2) then
         ! use line search in Newton direction
         call internal_chooseStepLength(lambdaOpt)
      end if

      if (MASTER) p = p0 + lambdaOpt*delta_p

      call myMPIBcastDouble(p,np)
      call wfparams_set(WFP,p,.true.)!ci coeefs are going to be normilized

      call ElocAndPsiTermsENR_reset(EPsiTENR)
      call internal_reCalcSample()
      pe0 = ElocAndPsiTermsENR_EmeanALL(EPsiTENR)
      pvar = ElocAndPsiTermsENR_varALL(EPsiTENR)
      if (MASTER) then
         ! true energy change for fixed sample
         delta_f = pe0 - e0
         r = delta_f/delta_q

         if (logmode >= 2) then
            write(iul,'(/A,F15.5,A,F15.5)') ' projected Emean =',pe0,' var = ',pvar
            write(iul,'(/3(a,f15.5))') ' r = ',r,'   delta_f = ',delta_f,'  e0 = ',e0
         end if

         ! adaptation of nu
         if (keepnu .eqv. .false.) then
            if (r < 0.25 .and. abs(delta_f) > deltaFmin) then
                 nu = 4*nu
                 if (logmode >= 2) write(iul,'(a,f15.6/)') ' -> increasing nu to ',nu
            else if (r > 0.75) then
                 nu = nu/2
                 if (logmode >= 2) write(iul,'(a,f15.6/)') ' -> decreasing nu to ',nu
            end if
              if (r < 0 .and. abs(delta_f) > deltaFmin) then
                 if (logmode >= 2) write(iul,'(/A/)') ' !!! going back to previous parameter vector !!!'
                 p = p0
            end if
         else
            nu=nuStart
         end if
      end if

      call myMPIBcastDouble(r,1)
      if (r < 0) then
         call myMPIBcastDouble(p,np)
         call wfparams_set(WFP,p,.true.)  !ci coeffs are going to be normalized
         call ElocAndPsiTermsENR_reset(EPsiTENR)
         call internal_reCalcSample()
         pe0 = ElocAndPsiTermsENR_EmeanALL(EPsiTENR)
         pvar = ElocAndPsiTermsENR_varALL(EPsiTENR)
         if (logmode >= 2) then
            write(iul,'(A,F15.5,A,F15.5)') ' going back: projected Emean =',pe0,' var = ',pvar
         end if
      end if
      if (logmode >= 2 .and. npCI>0) then
         write(iul,*) ''
         write(iul,*) 'ci coefficients are normalized'
         write(iul,*) ''
      endif

      if (doWritewf) then
         fname = trim(baseName)//'-'//trim(intToStr(optIter))//'.wf'
         call writeWF(fname,.false.,ecp)
      end if

   if (eqStep >= eqIter) exit

      eqStep = eqStep + 1
      ! call "subroutine" subName in .in or macro subName.cmd that should contain
      ! code for equilibrating the sample with the new wave function
      call subloop(subName,sample,exitSubloop)

      if (exitSubloop) exit

   end do

   call setCurrentResult(e0,0.d0,var)

!    deallocate(p,p0,delta_p,bb,g,H,H0,Imat,ipiv,work)
!    deallocate(fi,ELi,fiEL)
!    deallocate(fifj,fifjEL,fiELj,fij,fijEL)
!    deallocate(A,B,D)

   call ElocAndPsiTermsENR_destroy(EPsiTENR)

   contains

      subroutine internal_readInput()
      !-----------------------------!
         character(len=40) optMethod
         NRMode = 3
         call getstra(lines,nl,'method=',optMethod,iflag)
         if (iflag==0) then
            if (optMethod=='nr' .or. optMethod=='newton') then  ! Newton-Raphson
               NRMode = 1
            else if (optMethod=='scaled_nr' .or. optMethod=='scaled_newton') then
               NRMode = 2
            else if (optMethod=='snr' .or. optMethod=='lm_newton') then  ! stabilized Newton-Raphson
               NRMode = 3
            else if (optMethod=='lm') then  ! Levenberg-Marquardt
               NRMode = 4
            else
               call abortp("$optimize_parameters: illegal newton method name")
            end if
         end if
         gfac = 0
         call getdbla(lines,nl,'gfac=',gfac,iflag)
         call getdbla(lines,nl,'target_E=',targetE,iflag)
         if (iflag /= 0 .and. NRmode==2) call abortp('scaled_newton: target_E option required')
         call getdbla(lines,nl,'target_var=',targetVar,iflag)
         if (iflag /= 0 .and. NRmode==2) call abortp('scaled_newton: target_var option required')
         maxVar = 1.d9
         call getdbla(lines,nl,'max_var=',maxVar,iflag)
         dmax = 1.d9
         call getdbla(lines,nl,'dmax=',dmax,iflag)
         cffac = 1.0d0
         call getdbla(lines,nl,'cffac=',cffac,iflag)
         nuStart = 0.01d0
         call getdbla(lines,nl,'nu=',nuStart,iflag)
         deltaFmin = 0.0d0
         call getdbla(lines,nl,'delta_f_min=',deltaFmin,iflag)
         doWriteWF = finda(lines,nl,'write_wf')
         keepnu    = finda(lines,nl,'keep_nu')
         subName = 'equilibrate'
         call getstra(lines,nl,'eq_call=',subName,iflag)
         eqIter = 0
         call getinta(lines,nl,'eq_iter=',eqIter,iflag)
      end subroutine internal_readInput


      subroutine internal_calcEPsiTerms()
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for linear method
         real(r8) :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)

         rwp => getFirst(sample)
         do
            call pos(rwp,x,y,z)
            call eConfigArray_set(ec,1,x,y,z)
            call eloc(0, ec, WFP%optType, wfpDef=WFP, wfpDT=wfpDT)
            call ElocAndPsiTermsENR_add(EPsiTENR,wfpDT)
         if (.not.isNext(sample)) exit
            rwp => getNext(sample)
         enddo

      end subroutine internal_calcEPsiTerms

      subroutine internal_reCalcSample()
      !---------------------------------!
         ! calculate sample average for E_loc  without params drivatives and updates rw
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
            call ElocAndPsiTermsENR_add(EPsiTENR,wfpDT) ! for ElocAndPsiTermsENR_EmeanALL
         if (.not.isNext(sample)) exit
            rwp => getNext(sample)
         enddo
         wfpDT%noCalc=.false.
      end subroutine internal_reCalcSample


      subroutine internal_calcGradAndHessian(g,H)
         real(r8) :: g(:)
         real(r8) :: H(:,:)

         g = 2 * (fiEL - fi * e0)

         A = 2 * (fijEL - fij * e0 - fifjEL + fifj * e0)
         do j = 1, np
            B(:, j) = -2 * (fi(:) * g(j) + fi(j) * g(:))
            D(:, j) = - fi(:) * ELi(j) - fi(j) * ELi(:)
         enddo
         B = B + 4 * (fifjEL - fifj * e0)
         D = D + fiELj + TRANSPOSE(fiELj)

         H = A + B + D

         g = g + gfac * ELi
      end subroutine internal_calcGradAndHessian

      subroutine internal_chooseStepLength(lambdaOpt)
      !---------------------------------------------------------!
         ! choose the optimal step length based on a trust radius and a cost function
         integer :: idx       ! selected eigenvector
         real(r8)  :: e0,var0   ! E and Var for unit step length
         real(r8)  :: lambdaOpt ! calculated optimal step length
         integer n,cfIdx
         real(r8) d,cfmin,cf,var,cf1,cfs(6),x1,x2,x3,y1,y2,y3,r1,r2,r3

         lambda = (/ 0.02d0, 0.1d0, 0.3d0, 0.5d0, 0.7d0, 1.d0 /)

         d = sum(abs(delta_p(1:np)))/np
         if (logmode >= 2) write(iul,'(3(a,f15.5))')  &
            ' lambda=dmax/d=',dmax/d,' d=',d,' dmax=',dmax

         cfmin = 1.d99
         cf1 = 1.d99
         do n=1,size(lambda)

            if (MASTER) p = p0 + lambda(n)*delta_p(1:np)

            call myMPIBcastDouble(p,np)
            call wfparams_set(WFP,p)

            call ElocAndPsiTermsENR_reset(EPsiTENR)
            call internal_calcEPsiTerms()
            e0 = ElocAndPsiTermsENR_EmeanALL(EPsiTENR)
            var = ElocAndPsiTermsENR_varALL(EPsiTENR)

            cf = abs(e0-targetE) + cffac* abs(var-targetVar)
            cfs(n) = cf
            if (logmode >= 2) write(iul,'(I5,A,F10.2,A,F15.5,A,F15.5,A,F15.5)') &
               n,': lambda=',lambda(n),' Emean =',e0,' var = ',var,' cf = ',cf
            if (cf < cfmin) then
               cfmin = cf
               cfIdx = n
            end if
            if (lambda(n) == 1.d0) cf1=cf

         end do

         if (lambda(cfIdx) < dmax/d) then
            if ((cf1-cfmin)/cfmin > 0.05) then
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
            else
               lambdaOpt = 1.d0
            end if
            if (logmode >= 2) write(iul,'(a,f10.2)') ' choosing min cost function (5% threshold): lambda=',lambdaOpt
         else
            lambdaOpt = dmax/d
            if (logmode >= 2) write(iul,'(a,f10.2)') ' choosing trust radius: lambda=',lambdaOpt
         end if
      end subroutine internal_chooseStepLength

      subroutine internal_writeDataToLogFile()
         real(r8) maxgrad,meangrad
         integer i

         maxgrad = maxval(abs(g))
         meangrad = sum(abs(g))/np
         if (logmode >= 2) then
            write(iul,'(/2(A,G12.4))') ' gradient with abs mean = ',meangrad,' and abs max =',maxgrad
            write(iul,'(10g12.4)') (g(i),i=1,np)
         endif
         if (logmode >= 3) then
            write(iul, *)
            write(iul,'(A,G13.5,7G12.4)') 'ONE:',e0,2*(fiEL(1)-e0*fi(1)),2*(fiEL(1)-e0*fi(1))+ELi(1), &
              2*(fiEL(1)-e0*fi(1))+2*ELi(1),H(1,1),A(1,1),B(1,1),D(1,1)
            write(iul,*) 'Hessian:'
            do i=1,np
               write(iul,'(15G10.3)') H(i,:)
            enddo
            if (logmode >= 4) then
               write(iul,*) 'fi:'
               write(iul,'(10G10.3)') fi(:)
               write(iul,*) 'ELi:'
               write(iul,'(10G10.3)') ELi(:)
               write(iul,*) 'fiEL:'
               write(iul,'(10G10.3)') fiEL(:)
               write(iul,*) 'fij:'
               do i=1,np
                  write(iul,'(15G10.3)') fij(i,:)
               enddo
               write(iul,*) 'fijEL:'
               do i=1,np
                  write(iul,'(15G10.3)') fijEL(i,:)
               enddo
               write(iul,*) 'fifj:'
               do i=1,np
                  write(iul,'(15G10.3)') fifj(i,:)
               enddo
               write(iul,*) 'fifjEL:'
               do i=1,np
                  write(iul,'(15G10.3)') fifjEL(i,:)
               enddo
               write(iul,*) 'fiELj:'
               do i=1,np
                  write(iul,'(15G10.3)') fiELj(i,:)
               enddo
               write(iul,*) 'A:'
               do i=1,np
                  write(iul,'(15G10.3)') A(i,:)
               enddo
               write(iul,*) 'B:'
               do i=1,np
                  write(iul,'(15G10.3)') B(i,:)
               enddo
               write(iul,*) 'D:'
               do i=1,np
                  write(iul,'(15G10.3)') D(i,:)
               enddo
            endif
         endif
      end subroutine internal_writeDataToLogFile

      subroutine internal_writeNewVectorToLogFile()
         if (logmode >= 2) then
            write(iul,*)
            if (logmode >= 3) then
               write(iul,*) ' initial parameter vector:'
               write(iul,'(10g12.4)') p
               write(iul,*) ' delta_p:'
               write(iul,'(10g12.4)') delta_p
            endif
            write(iul,*) ' new parameter vector:'
            write(iul,'(10g12.4)') p + delta_p
         endif
      end subroutine internal_writeNewVectorToLogFile

   end subroutine eminNR_optimizeSample


end module optParamsENR_m
