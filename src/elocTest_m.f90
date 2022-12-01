! Copyright (C) 1996, 2015, 2018 Arne Luechow
! Copyright (C) 2011, 2013 Alexander Sturm
!
! SPDX-License-Identifier: GPL-3.0-or-later

MODULE elocTest_m

   use kinds_m, only : r8
   use globalUtils_m, only: iul, logmode
   use rwSample_m
   use wfdata_m
   use parsing_m
   use elocData_m
   use eConfigs_m, only: eConfigArray, eConfigArray_new, eConfigArray_set, eConfigArray_destroy, eConfigArray_getPtr
   use eloc_m, only : eloc
   use numderivs_m

   implicit none

   private
   public :: runEloctest


CONTAINS

   subroutine runEloctest(lines, nl, sample)

      integer, intent(in) :: nl
      character(len = 120), intent(in) :: lines(nl)
      type(rwSample), intent(inout) :: sample

      integer :: points = 7          ! points for differentiation rule
      integer :: sampleSize = 0      ! number of points from sample
      real(r8) :: h = 1.d-5           ! denominator for numerical differentiation
      logical :: numDerivs = .false. ! do numerical derivatives?
      logical :: extrapDerivs = .false.
      logical :: extrap1Derivs = .false.
      logical :: update = .false.
      logical passed, allPassed
      real(r8) x(ne), y(ne), z(ne), tolEloc, tolGrad
      integer iflag, i, savelogmode, verb
      type(randomWalker), pointer :: rwp

      call assert(ne>0, &
         ' runEloctest: wave function must be initialized')
      call assert(getSampleSize(sample) > 0, &
         ' runEloctest: no sample available')

      allPassed = .true.
      savelogmode = logmode
      h = 1.d-5
      call getdbla(lines, nl, 'h=', h, iflag)
      points = 7
      call getinta(lines, nl, 'rule=', points, iflag)
      sampleSize = 1
      call getinta(lines, nl, 'points=', sampleSize, iflag)
      numDerivs = .false.
      numDerivs = finda(lines, nl, 'derivatives')
      extrapDerivs = finda(lines, nl, 'extrap_derivs')
      extrap1Derivs = finda(lines, nl, 'extrap1_derivs')
      update = .false.
      update = finda(lines, nl, 'updatetest')
      tolGrad = 1.e-11_r8
      tolEloc = 1.e-8_r8
      call getdbla(lines, nl, 'tol_grad=', tolGrad, iflag)
      call getdbla(lines, nl, 'tol_eloc=', tolEloc, iflag)
      call getinta(lines, nl, 'verbose=', verb, iflag)
      if (iflag == 0) logmode = verb

      if (logmode >= 2) then
         write(iul, '(/A/)') '  * * *  Eloctest  * * *'
         write(iul, *) ' calculating psi and E_local '
         write(iul, *) '  with contributions for sample points'
      end if

      rwp => getFirst(sample)
      do i = 1, sampleSize
         call pos(rwp, x, y, z)
         call elocContribs(x, y, z)
         if (numDerivs) then
            call elocDerivTest(x, y, z, h, points, tolGrad, tolEloc, passed)
            if (.not.passed) allPassed = .false.
         end if
         if (extrapDerivs) then
            call elocDerivExtrapTest(x, y, z, h, tolGrad, tolEloc, passed)
            if (.not.passed) allPassed = .false.
         end if
         if (extrap1Derivs) then
            call elocDerivExtrap1Test(x, y, z, h, tolGrad, tolEloc, passed)
            if (.not.passed) allPassed = .false.
         end if
         if (update) then
            call updatetest(x, y, z)
         end if
         if (.not.isNext(sample)) exit
         rwp => getNext(sample)
      end do

      write(iul,*)
      if (allPassed) then
         write(iul,*) ' all eloctests passed!'
      else
         write(iul,*) ' eloctests failed!'
      end if

      logmode = savelogmode

   end subroutine


   subroutine elocContribs(x, y, z)

      real(r8), intent(in) :: x(:), y(:), z(:)

      integer i
      real(r8) psiv, vpot, psilapl
      real(r8) jas, phi, jlapl, flapl, jgrad(3 * ne), fgrad(3 * ne)
      real(r8) expU
      type(eConfigArray) :: ec

      if (logmode >= 3) then
         write(iul, '(//)')
         write(iul, *) ' electron coordinates:'
         do i = 1, ne
            write(iul, '(i5,3F8.3)') i, x(i), y(i), z(i)
         enddo
      endif

      call eConfigArray_new(ec, ne, 1)
      call eConfigArray_set(ec, 1, x, y, z)
      call eloc(0, ec, 'none')

      if (logmode >= 3) then
         write(iul, '(3(A8,G14.8))') ' Eloc = ', elEloc(1), &
            ' Phi = ', elPhi(1), ' U = ', elU(1)
         write(iul, '(3(A8,G14.8))') ' Psi = ', elPhi(1) * exp(elU(1)), &
            ' Vpot = ', elVen(1) + elVee(1) + vpot0, ' Vpp = ', elECPPot
         write(iul, '(A8,G14.8)') ' Vnuc = ', vpot0
      end if

      expU = exp(elU(1))
      vpot = elVen(1) + elVee(1)
      jas = elU(1)
      phi = elPhi(1)
      psiv = elPhi(1) * expU
      jlapl = elUlapl(1)
      flapl = elFlapl(1)
      do i = 1, ne
         jgrad(3 * i - 2) = elUgrad(3 * i - 2, 1)
         jgrad(3 * i - 1) = elUgrad(3 * i - 1, 1)
         jgrad(3 * i) = elUgrad(3 * i, 1)
         fgrad(3 * i - 2) = elFgrad(3 * i - 2, 1)
         fgrad(3 * i - 1) = elFgrad(3 * i - 1, 1)
         fgrad(3 * i) = elFgrad(3 * i, 1)
      enddo
      psilapl = flapl * jas + phi * jlapl + 2 * dot_product(fgrad, jgrad)

      if (logmode >= 3) then
         write(iul, '(2(A8,G17.7))') ' jas =   ', elU(1), ' jlapl = ', jlapl
         write(iul, '(3(A8,G17.7))') ' Ekin1 = ', elEloc(1) - elVen(1) - elVee(1) - vpot0, &
            ' Ekin2 = ', -0.5d0 * psilapl / psiv, ' Ekin3 = ', sum(elEkini(1:ne, 1))

         write(iul, *) ' grad Psi / Psi:'
         do i = 1, ne
            write(iul, *) i, elxDrift(i, 1), elyDrift(i, 1), elzDrift(i, 1)
         enddo
         write(iul, *) ' grad Phi / Phi:'
         do i = 1, ne
            write(iul, *) i, fgrad(3 * i - 2) / phi, fgrad(3 * i - 1) / phi, fgrad(3 * i) / phi
         enddo
         write(iul, *) ' grad U:'
         do i = 1, ne
            write(iul, *) i, jgrad(3 * i - 2), jgrad(3 * i - 1), jgrad(3 * i)
         enddo
      endif

      call eConfigArray_destroy(ec)

   end subroutine


   subroutine elocDerivtest(xx, yy, zz, h, points, tolGrad, tolEloc, passed)

      real(r8), intent(inout) :: xx(:), yy(:), zz(:)
      real(r8), intent(in)    :: h                    ! denominator in numerical derivatives
      integer, intent(in)     :: points               ! 3|5|7|9 point rule
      real(r8), intent(in)    :: tolGrad, tolEloc     ! tolerance for max(abs_grad_error) and max(abs_eloc_error)
      logical, intent(out)    :: passed               ! .true. if abs_errors smaller than tolGrad and tolEloc
      integer i, k, pnts
      real(r8), pointer :: x(:), y(:), z(:)
      real(r8) psiv, expU
      real(r8) jas, phi, jlapl, flapl, jgrad(3 * ne), fgrad(3 * ne)
      real(r8) psi0(8, 3 * ne), fx(3 * ne), fxx(3 * ne), f1
      real(r8) psif(8, 3 * ne), ff1
      real(r8) psij(8, 3 * ne), fj1
      real(r8) vpots0, vpot, grad(3 * ne)
      real(r8) gr0(3 * ne), grj0(3 * ne), grf0(3 * ne)
      real(r8) numEloc, eloc0, lapl0, flapl0, jlapl0
      real(r8) maxError, maxErrorPhi, maxErrorU
      type(eConfigArray) :: ec

      passed = .true.

      call eConfigArray_new(ec, ne, 1)
      call eConfigArray_set(ec, 1, xx, yy, zz)
      call eConfigArray_getPtr(ec, 1, x, y, z)   ! x,y,z pointer to ec positions
      call eloc(0, ec, 'none')

      maxError = 0
      maxErrorPhi = 0
      maxErrorU = 0
      expU = exp(elU(1))
      vpot = elVen(1) + elVee(1)
      jas = elU(1)
      phi = elPhi(1)
      psiv = elPhi(1) * expU
      jlapl = elUlapl(1)
      flapl = elFlapl(1)
      do i = 1, ne
         grad(3 * i - 2) = elxDrift(i, 1)
         grad(3 * i - 1) = elyDrift(i, 1)
         grad(3 * i) = elzDrift(i, 1)
         jgrad(3 * i - 2) = elUgrad(3 * i - 2, 1)
         jgrad(3 * i - 1) = elUgrad(3 * i - 1, 1)
         jgrad(3 * i) = elUgrad(3 * i, 1)
         fgrad(3 * i - 2) = elFgrad(3 * i - 2, 1)
         fgrad(3 * i - 1) = elFgrad(3 * i - 1, 1)
         fgrad(3 * i) = elFgrad(3 * i, 1)
      end do


      !     // original values
      f1 = psiv
      fj1 = jas
      ff1 = phi
      eloc0 = elEloc(1)
      jlapl0 = jlapl
      flapl0 = flapl
      vpots0 = elVen(1) + elVee(1) + vpot0
      lapl0 = -2d0 * (elEloc(1) - vpots0) * psiv
      do i = 1, 3 * ne
         gr0(i) = grad(i) * psiv
         grj0(i) = jgrad(i)
         grf0(i) = fgrad(i)
      end do

      !     // Calculate points for numerical derivaties

      do k = 1, ne
         do pnts = 1, points / 2
            x(k) = x(k) + pnts * h
            call eloc(0, ec, 'none')
            psi0(2 * pnts - 1, 3 * k - 2) = elPhi(1) * exp(elU(1))
            psij(2 * pnts - 1, 3 * k - 2) = elU(1)
            psif(2 * pnts - 1, 3 * k - 2) = elPhi(1)
            x(k) = x(k) - 2 * pnts * h
            call eloc(0, ec, 'none')
            psi0(2 * pnts, 3 * k - 2) = elPhi(1) * exp(elU(1))
            psij(2 * pnts, 3 * k - 2) = elU(1)
            psif(2 * pnts, 3 * k - 2) = elPhi(1)
            x(k) = x(k) + pnts * h

            y(k) = y(k) + pnts * h
            call eloc(0, ec, 'none')
            psi0(2 * pnts - 1, 3 * k - 1) = elPhi(1) * exp(elU(1))
            psij(2 * pnts - 1, 3 * k - 1) = elU(1)
            psif(2 * pnts - 1, 3 * k - 1) = elPhi(1)
            y(k) = y(k) - 2 * pnts * h
            call eloc(0, ec, 'none')
            psi0(2 * pnts, 3 * k - 1) = elPhi(1) * exp(elU(1))
            psij(2 * pnts, 3 * k - 1) = elU(1)
            psif(2 * pnts, 3 * k - 1) = elPhi(1)
            y(k) = y(k) + pnts * h

            z(k) = z(k) + pnts * h
            call eloc(0, ec, 'none')
            psi0(2 * pnts - 1, 3 * k) = elPhi(1) * exp(elU(1))
            psij(2 * pnts - 1, 3 * k) = elU(1)
            psif(2 * pnts - 1, 3 * k) = elPhi(1)
            z(k) = z(k) - 2 * pnts * h
            call eloc(0, ec, 'none')
            psi0(2 * pnts, 3 * k) = elPhi(1) * exp(elU(1))
            psij(2 * pnts, 3 * k) = elU(1)
            psif(2 * pnts, 3 * k) = elPhi(1)
            z(k) = z(k) + pnts * h
         end do
      end do

      ! calculate numerical derivatives
      ! second order formula for 1st derivative
      ! 3/5/7 point formula dependent on "points" for 2nd derivative

      ! total wavefunction
      flapl = 0d0
      do k = 1, 3 * ne
         select case (points)
         case (3)
            fx(k) = (psi0(1, k) - psi0(2, k)) / (2d0 * h)
            fxx(k) = (-2d0 * f1 + psi0(2, k) + psi0(1, k)) / h**2
         case (5)
            fx(k) = (8d0 * (psi0(1, k) - psi0(2, k)) - &
                  (psi0(3, k) - psi0(4, k))) / (12d0 * h)
            fxx(k) = (-30d0 * f1 + 16d0 * (psi0(1, k) + psi0(2, k))&
                  - (psi0(3, k) + psi0(4, k))) / (12d0 * h**2)
         case (7)
            fx(k) = (45d0 * (psi0(1, k) - psi0(2, k)) - &
                  9d0 * (psi0(3, k) - psi0(4, k)) + &
                  (psi0(5, k) - psi0(6, k))) / (60d0 * h)
            fxx(k) = (-490d0 * f1 + 270d0 * (psi0(1, k) + psi0(2, k))&
                  - 27d0 * (psi0(3, k) + psi0(4, k))&
                  + 2d0 * (psi0(5, k) + psi0(6, k))) / (180d0 * h**2)
         case (9)
            fx(k) = (672d0 * (psi0(1, k) - psi0(2, k)) - &
                  168d0 * (psi0(3, k) - psi0(4, k)) + &
                  32d0 * (psi0(5, k) - psi0(6, k)) - &
                  3d0 * (psi0(7, k) - psi0(8, k))) / (840d0 * h)
            fxx(k) = (-71750d0 * f1 + 40320d0 * (psi0(1, k) + psi0(2, k))&
                  - 5040d0 * (psi0(3, k) + psi0(4, k))&
                  + 640d0 * (psi0(5, k) + psi0(6, k))&
                  - 45d0 * (psi0(7, k) + psi0(8, k))) / (25200d0 * h**2)
         case default
            call abortp('(eloctest): wrong input for points')
         end select

         flapl = flapl + fxx(k)

      end do
      numEloc = -0.5d0 * flapl / f1 + vpots0

      maxError = 0
      do i = 1, 3 * ne
          maxError = max( maxError, abs((fx(i) - gr0(i))/gr0(i)) )
      end do

      if (maxError > tolGrad .or. abs(eloc0 - numEloc) > tolEloc) passed = .false.

      if (.not.passed .or. logmode >= 3) then

         write(iul, *)
         if (passed) then
            write(iul, *) ' numerical gradient and eloc test  *** passed *** '
         else 
            write(iul, *) ' numerical gradient and eloc test *** failed ***'
         end if
         write(iul, *)
         write(iul, *) ' numerical derivatives for Psi:'
         write(iul, *) ' computed, numerical, abs error, rel error'
         write(iul, '(A,4G17.7)') ' laplacian:', lapl0, flapl, lapl0 - flapl, &
             (lapl0 - flapl) / lapl0
         write(iul, '(A,4G17.7)') ' Eloc:', eloc0, numEloc, eloc0 - numEloc, &
             (numEloc - eloc0) / eloc0
         write(iul, *) ' gradient:'
         do i = 1, 3 * ne
            write(iul, '(i5,4g17.7)') i, gr0(i), fx(i), (fx(i) - gr0(i)), &
                     (fx(i) - gr0(i)) / gr0(i)
         end do
         write(iul, *) ' maximal relative gradient error:', maxError

         !     // Jastrow part
         flapl = 0d0
         do k = 1, 3 * ne
            select case (points)
            case (3)
               fx(k) = (psij(1, k) - psij(2, k)) / (2d0 * h)
               fxx(k) = (-2d0 * fj1 + psij(2, k) + psij(1, k)) / h**2
            case (5)
               fx(k) = (8d0 * (psij(1, k) - psij(2, k)) - &
                     (psij(3, k) - psij(4, k))) / (12d0 * h)
               fxx(k) = (-30d0 * fj1 + 16d0 * (psij(1, k) + psij(2, k))&
                     - (psij(3, k) + psij(4, k))) / (12d0 * h**2)
            case (7)
               fx(k) = (45d0 * (psij(1, k) - psij(2, k)) - &
                     9d0 * (psij(3, k) - psij(4, k)) + &
                     (psij(5, k) - psij(6, k))) / (60d0 * h)
               fxx(k) = (-490d0 * fj1 + 270d0 * (psij(1, k) + psij(2, k))&
                     - 27d0 * (psij(3, k) + psij(4, k))&
                     + 2d0 * (psij(5, k) + psij(6, k))) / (180d0 * h**2)
            case (9)
               fx(k) = (672d0 * (psij(1, k) - psij(2, k)) - &
                     168d0 * (psij(3, k) - psij(4, k)) + &
                     32d0 * (psij(5, k) - psij(6, k)) - &
                     3d0 * (psij(7, k) - psij(8, k))) / (840d0 * h)
               fxx(k) = (-71750d0 * fj1 + 40320d0 * (psij(1, k) + psij(2, k))&
                     - 5040d0 * (psij(3, k) + psij(4, k))&
                     + 640d0 * (psij(5, k) + psij(6, k))&
                     - 45d0 * (psij(7, k) + psij(8, k))) / (25200d0 * h**2)
            case default
               call abortp('(elocwf): wrong input for points')
            end select

            flapl = flapl + fxx(k)
         end do

         write(iul, *) ' '
         write(iul, *) ' JASTROW part: numerical derivatives'
         write(iul, '(A,4G17.7)') ' laplacian:', jlapl0, flapl, jlapl0 - flapl, &
             (jlapl0 - flapl) / jlapl0
         write(iul, *) ' gradient:'
         do i = 1, 3 * ne
            write(iul, '(I5,4G17.7)') i, grj0(i), fx(i), (fx(i) - grj0(i)), &
                     (fx(i) - grj0(i)) / grj0(i)
         enddo
         maxError = 0
         do i = 1, 3 * ne
             maxError = max(maxError, abs((fx(i) - grj0(i)) / grj0(i)))
         enddo
         write(iul, *) ' maximal relative gradient error:', maxError

         !     // Phi part
         flapl = 0d0
         do k = 1, 3 * ne
            select case (points)
            case (3)
               fx(k) = (psif(1, k) - psif(2, k)) / (2d0 * h)
               fxx(k) = (-2d0 * ff1 + psif(2, k) + psif(1, k)) / h**2
            case (5)
               fx(k) = (8d0 * (psif(1, k) - psif(2, k)) - &
                     (psif(3, k) - psif(4, k))) / (12d0 * h)
               fxx(k) = (-30d0 * ff1 + 16d0 * (psif(1, k) + psif(2, k))&
                     - (psif(3, k) + psif(4, k))) / (12d0 * h**2)
            case (7)
               fx(k) = (45d0 * (psif(1, k) - psif(2, k)) - &
                     9d0 * (psif(3, k) - psif(4, k)) + &
                     (psif(5, k) - psif(6, k))) / (60d0 * h)
               fxx(k) = (-490d0 * ff1 + 270d0 * (psif(1, k) + psif(2, k))&
                     - 27d0 * (psif(3, k) + psif(4, k))&
                     + 2d0 * (psif(5, k) + psif(6, k))) / (180d0 * h**2)
            case (9)
               fx(k) = (672d0 * (psif(1, k) - psif(2, k)) - &
                     168d0 * (psif(3, k) - psif(4, k)) + &
                     32d0 * (psif(5, k) - psif(6, k)) - &
                     3d0 * (psif(7, k) - psif(8, k))) / (840d0 * h)
               fxx(k) = (-71750d0 * ff1 + 40320d0 * (psif(1, k) + psif(2, k))&
                     - 5040d0 * (psif(3, k) + psif(4, k))&
                     + 640d0 * (psif(5, k) + psif(6, k))&
                     - 45d0 * (psif(7, k) + psif(8, k))) / (25200d0 * h**2)
            case default
               call abortp('(eloctest): wrong input for points')
            end select

            flapl = flapl + fxx(k)

         end do
         write(iul, *) ' '
         write(iul, *) ' DETERMINANT PART: numerical derivatives:'
         write(iul, '(A,4G17.7)') ' laplacian:', flapl0, flapl, flapl0 - flapl, &
             (flapl0 - flapl) / flapl0
         write(iul, *) ' gradient:'
         do i = 1, 3 * ne
            write(iul, '(I5,4G17.7)') i, grf0(i), fx(i), (fx(i) - grf0(i)), &
                     (fx(i) - grf0(i)) / grf0(i)
         enddo
         maxError = 0
         do i = 1, 3 * ne
             maxError = max(maxError, abs((fx(i) - grf0(i)) / grf0(i)))
         enddo
         write(iul, *) ' maximal relative gradient error:', maxError
      end if

      call eConfigArray_destroy(ec)

   end subroutine elocDerivTest



   subroutine elocDerivExtrapTest(xx, yy, zz, h0, tolGrad, tolEloc, passed)

      real(r8), intent(inout) :: xx(:), yy(:), zz(:)
      real(r8), intent(in)    :: h0                   ! initial denominator in numerical derivatives
      real(r8), intent(in)    :: tolGrad, tolEloc     ! tolerance for max(abs_grad_error) and max(abs_eloc_error)
      logical, intent(out)    :: passed               ! .true. if abs_errors smaller than tolGrad and tolEloc
      integer i, k, n, nMax
      real(r8), pointer :: x(:), y(:), z(:)
      real(r8) psiv, expU, h, d2ydx2
      real(r8) jas, phi, lapl, ulapl, flapl, jgrad(3 * ne), fgrad(3 * ne)
      real(r8) psi0(3, 3 * ne), f1
      real(r8) gradient(3 * ne), epsGradient(3 * ne),  epsLapl

      real(r8) psif(3, 3 * ne), ff1
      real(r8) psij(3, 3 * ne), fj1
      real(r8) vpots0, vpot, grad(3 * ne)
      real(r8) gr0(3 * ne), grj0(3 * ne), grf0(3 * ne)
      real(r8) numEloc, eloc0, lapl0, flapl0, jlapl0
      real(r8) maxError, maxErrorPhi, maxErrorU
      real(r8), allocatable :: t0grad(:,:), t0gradj(:,:), t0gradf(:,:)
      real(r8), allocatable :: t0lapl(:), t0laplu(:), t0laplf(:)
      real(r8), allocatable :: characteristicScale(:), characteristicScalej(:), characteristicScalef(:)
      type(eConfigArray) :: ec

      passed = .true.

      call eConfigArray_new(ec, ne, 1)
      call eConfigArray_set(ec, 1, xx, yy, zz)
      call eConfigArray_getPtr(ec, 1, x, y, z)   ! x,y,z pointer to ec positions
      call eloc(0, ec, 'none')

      maxError = 0
      maxErrorPhi = 0
      maxErrorU = 0
      expU = exp(elU(1))
      vpot = elVen(1) + elVee(1)
      jas = elU(1)
      phi = elPhi(1)
      psiv = elPhi(1) * expU
      do i = 1, ne
         grad(3 * i - 2) = elxDrift(i, 1)
         grad(3 * i - 1) = elyDrift(i, 1)
         grad(3 * i) = elzDrift(i, 1)
         jgrad(3 * i - 2) = elUgrad(3 * i - 2, 1)
         jgrad(3 * i - 1) = elUgrad(3 * i - 1, 1)
         jgrad(3 * i) = elUgrad(3 * i, 1)
         fgrad(3 * i - 2) = elFgrad(3 * i - 2, 1)
         fgrad(3 * i - 1) = elFgrad(3 * i - 1, 1)
         fgrad(3 * i) = elFgrad(3 * i, 1)
      end do


      !  analytical values
      f1 = psiv
      fj1 = jas
      ff1 = phi
      eloc0 = elEloc(1)
      jlapl0 = elUlapl(1)
      flapl0 = elFlapl(1)
      vpots0 = elVen(1) + elVee(1) + vpot0
      lapl0 = -2 * (elEloc(1) - vpots0) * psiv
      do i = 1, 3 * ne
         gr0(i) = grad(i) * psiv
         grj0(i) = jgrad(i)
         grf0(i) = fgrad(i)
      end do

      !  Calculate points for numerical derivaties

      nMax = 3
      allocate(t0grad(0 : nMax, 3 * ne), t0gradj(0 : nMax, 3 * ne), t0gradf(0 : nMax, 3 * ne))
      allocate(t0lapl(0 : nMax), t0laplu(0 : nMax), t0laplf(0 : nMax))
      allocate(characteristicScale(3 * ne), characteristicScalej(3 * ne), characteristicScalef(3 * ne))

      h = h0
      do n = 0, nMax
         lapl = 0
         ulapl = 0
         flapl = 0
         do k = 1, ne
            x(k) = x(k) + h
            call eloc(0, ec, 'none')
            psi0(3, 3 * k - 2) = elPhi(1) * exp(elU(1))
            psij(3, 3 * k - 2) = elU(1)
            psif(3, 3 * k - 2) = elPhi(1)

            psi0(2, 3 * k - 2) = psiv
            psij(2, 3 * k - 2) = jas
            psif(2, 3 * k - 2) = phi

            x(k) = x(k) - 2._r8 * h
            call eloc(0, ec, 'none')
            psi0(1, 3 * k - 2) = elPhi(1) * exp(elU(1))
            psij(1, 3 * k - 2) = elU(1)
            psif(1, 3 * k - 2) = elPhi(1)
            x(k) = x(k) + h

            t0grad(n, 3 * k - 2) = numFirstDerivative3Point(psi0(:, 3 * k - 2), h)
            d2ydx2 = numSecondDerivative3Point(psi0(:, 3 * k - 2), h)
            if (n == 0) characteristicScale(3 * k - 2) = sqrt( abs( psiv / d2ydx2 ))
            lapl = lapl + d2ydx2

            t0gradj(n, 3 * k - 2) = numFirstDerivative3Point(psij(:, 3 * k - 2), h)
            d2ydx2 = numSecondDerivative3Point(psij(:, 3 * k - 2), h)
            if (n == 0) characteristicScalej(3 * k - 2) = sqrt( abs( jas / d2ydx2 ))
            ulapl = ulapl + d2ydx2

            t0gradf(n, 3 * k - 2) = numFirstDerivative3Point(psif(:, 3 * k - 2), h)
            d2ydx2 = numSecondDerivative3Point(psif(:, 3 * k - 2), h)
            if (n == 0) characteristicScalej(3 * k - 2) = sqrt( abs( phi / d2ydx2 ))
            flapl = flapl + d2ydx2

            y(k) = y(k) + h
            call eloc(0, ec, 'none')
            psi0(3, 3 * k - 1) = elPhi(1) * exp(elU(1))
            psij(3, 3 * k - 1) = elU(1)
            psif(3, 3 * k - 1) = elPhi(1)

            psi0(2, 3 * k - 1) = psiv
            psij(2, 3 * k - 1) = jas
            psif(2, 3 * k - 1) = phi

            y(k) = y(k) - 2._r8 * h
            call eloc(0, ec, 'none')
            psi0(1, 3 * k - 1) = elPhi(1) * exp(elU(1))
            psij(1, 3 * k - 1) = elU(1)
            psif(1, 3 * k - 1) = elPhi(1)
            y(k) = y(k) + h

            t0grad(n, 3 * k - 1) = numFirstDerivative3Point(psi0(:, 3 * k - 1), h)
            d2ydx2 = numSecondDerivative3Point(psi0(:, 3 * k - 1), h)
            if (n == 0) characteristicScale(3 * k - 1) = sqrt( abs( psiv / d2ydx2 ))
            lapl = lapl + d2ydx2

            t0gradj(n, 3 * k - 1) = numFirstDerivative3Point(psij(:, 3 * k - 1), h)
            d2ydx2 = numSecondDerivative3Point(psij(:, 3 * k - 1), h)
            if (n == 0) characteristicScalej(3 * k - 1) = sqrt( abs( jas / d2ydx2 ))
            ulapl = ulapl + d2ydx2

            t0gradf(n, 3 * k - 1) = numFirstDerivative3Point(psif(:, 3 * k - 1), h)
            d2ydx2 = numSecondDerivative3Point(psif(:, 3 * k - 1), h)
            if (n == 0) characteristicScalej(3 * k - 1) = sqrt( abs( phi / d2ydx2 ))
            flapl = flapl + d2ydx2

            z(k) = z(k) +  h
            call eloc(0, ec, 'none')
            psi0(3, 3 * k) = elPhi(1) * exp(elU(1))
            psij(3, 3 * k) = elU(1)
            psif(3, 3 * k) = elPhi(1)

            psi0(2, 3 * k) = psiv
            psij(2, 3 * k) = jas
            psif(2, 3 * k) = phi

            z(k) = z(k) - 2._r8 * h
            call eloc(0, ec, 'none')
            psi0(1, 3 * k) = elPhi(1) * exp(elU(1))
            psij(1, 3 * k) = elU(1)
            psif(1, 3 * k) = elPhi(1)
            z(k) = z(k) + h

            t0grad(n, 3 * k) = numFirstDerivative3Point(psi0(:, 3 * k), h)
            d2ydx2 = numSecondDerivative3Point(psi0(:, 3 * k), h)
            if (n == 0) characteristicScale(3 * k) = sqrt( abs( psiv / d2ydx2 ))
            lapl = lapl + d2ydx2

            t0gradj(n, 3 * k) = numFirstDerivative3Point(psij(:, 3 * k), h)
            d2ydx2 = numSecondDerivative3Point(psij(:, 3 * k), h)
            if (n == 0) characteristicScalej(3 * k) = sqrt( abs( jas / d2ydx2 ))
            ulapl = ulapl + d2ydx2

            t0gradf(n, 3 * k) = numFirstDerivative3Point(psif(:, 3 * k), h)
            d2ydx2 = numSecondDerivative3Point(psif(:, 3 * k), h)
            if (n == 0) characteristicScalej(3 * k) = sqrt( abs( phi / d2ydx2 ))
            flapl = flapl + d2ydx2
         end do

         t0lapl(n) = lapl
         t0laplu(n) = ulapl
         t0laplf(n) = flapl

         h = h / 2

      end do

      call nevilleExtrapolation(t0lapl, lapl, epsLapl)
      numEloc = -0.5d0 * lapl / f1 + vpots0

      maxError = 0
      do i = 1, 3 * ne
         call nevilleExtrapolation(t0grad(0:, i), gradient(i), epsGradient(i))
         maxError = max( maxError, abs((gradient(i) - gr0(i))/gr0(i)) )
      end do

      if (maxError > tolGrad .or. abs(eloc0 - numEloc) > tolEloc) passed = .false.

      if (logmode >= 2) then
         write(iul, *)
         if (passed) then
            write(iul, *) ' numerical gradient and eloc test (Neville) *** passed *** '
         else
            write(iul, *) ' numerical gradient and eloc test (Neville) *** failed ***'
            write(iul, *) maxError, tolGrad, abs(eloc0 - numEloc), tolEloc
         end if
      end if

      if (.not.passed .or. logmode >= 3) then

         write(iul, *)
         write(iul, *) ' numerical derivatives for Psi:'
         write(iul, *) ' computed, numerical, abs error, rel error, extrap error, char scale'
         write(iul, '(A,5g20.10)') ' laplacian:', lapl0, lapl, lapl0 - lapl, &
                 (lapl0 - lapl) / lapl0, epsLapl
         write(iul, '(A,4g20.10)') ' Eloc:', eloc0, numEloc, eloc0 - numEloc, &
                 (numEloc - eloc0) / eloc0
         write(iul, *) ' gradient:'
         do i = 1, 3 * ne
            write(iul, '(i5,6g20.10)') i, gr0(i), gradient(i), (gradient(i) - gr0(i)), &
                    (gradient(i) - gr0(i)) / gr0(i), epsGradient(i), characteristicScale(i)
         end do
         write(iul, *) ' maximal relative gradient error:', maxError


         ! Jastrow U errors

         call nevilleExtrapolation(t0laplu, ulapl, epsLapl)

         maxError = 0
         do i = 1, 3 * ne
            call nevilleExtrapolation(t0gradj(0:, i), gradient(i), epsGradient(i))
            maxError = max( maxError, abs((gradient(i) - grj0(i))/grj0(i)) )
         end do

         write(iul, *) ' '
         write(iul, *) ' JASTROW part: numerical derivatives'
         write(iul, '(a,5g20.10)') ' laplacian:', jlapl0, ulapl, jlapl0 - ulapl, &
                 (jlapl0 - ulapl) / jlapl0, epsLapl
         write(iul, *) ' gradient:'
         do i = 1, 3 * ne
            write(iul, '(i5,6g20.10)') i, grj0(i), gradient(i), (gradient(i) - grj0(i)), &
                    (gradient(i) - grj0(i)) / grj0(i), epsGradient(i), characteristicScalej(i)
         enddo
         write(iul, *) ' maximal relative U gradient error:', maxError

         ! Phi part

         call nevilleExtrapolation(t0laplf, flapl, epsLapl)

         maxError = 0
         do i = 1, 3 * ne
            call nevilleExtrapolation(t0gradf(0:, i), gradient(i), epsGradient(i))
            maxError = max( maxError, abs((gradient(i) - grf0(i))/grf0(i)) )
         end do

         write(iul, *) ' '
         write(iul, *) ' DETERMINANT PART: numerical derivatives:'
         write(iul, '(A,5g20.10)') ' laplacian:', flapl0, flapl, flapl0 - flapl, &
                 (flapl0 - flapl) / flapl0, epsLapl
         write(iul, *) ' gradient:'
         do i = 1, 3 * ne
            write(iul, '(i5,6g20.10)') i, grf0(i), gradient(i), (gradient(i) - grf0(i)), &
                    (gradient(i) - grf0(i)) / grf0(i), epsGradient(i), characteristicScalef(i)
         enddo
         write(iul, *) ' maximal relative det gradient error:', maxError
      end if

      call eConfigArray_destroy(ec)

   end subroutine elocDerivExtrapTest


   subroutine elocDerivExtrap1Test(xx, yy, zz, h0, tolGrad, tolEloc, passed)

      real(r8), intent(inout) :: xx(:), yy(:), zz(:)
      real(r8), intent(in)    :: h0                   ! initial denominator in numerical derivatives
      real(r8), intent(in)    :: tolGrad, tolEloc     ! tolerance for max(abs_grad_error) and max(abs_eloc_error)
      logical, intent(out)    :: passed               ! .true. if abs_errors smaller than tolGrad and tolEloc
      integer i, j, k, n, nMax, iter, iteru, iterf, iterl, iterlu, iterlf
      real(r8), pointer :: x(:), y(:), z(:)
      real(r8) psiv, expU, h
      real(r8) jas, phi, lapl, ulapl, flapl, jgrad(3 * ne), fgrad(3 * ne), lapli(ne), flapli(ne), ulapli(ne)
      real(r8) psi0(3)
      real(r8) gradient, epsGradient, ugradient, epsUGradient, fgradient, epsFGradient
      real(r8) epsLapl, epsULapl, epsFLapl

      real(r8) psif(3)
      real(r8) psij(3)
      real(r8) vpots0, vpot, grad(3 * ne)
      real(r8) numEloc, eloc0, lapl0, flapl0, jlapl0, lapl1, ulapl1, flapl1, total_lapl, total_ulapl, total_flapl
      real(r8) maxError, maxFError, maxUError, loss, loss0
      real(r8), allocatable :: t0grad(:), t0gradj(:), t0gradf(:)
      real(r8), allocatable :: t0lapl(:), t0laplu(:), t0laplf(:)
      real(r8) :: characteristicScale, characteristicScaleJ, characteristicScaleF
      type(eConfigArray) :: ec

      passed = .true.

      call eConfigArray_new(ec, ne, 1)
      call eConfigArray_set(ec, 1, xx, yy, zz)
      call eConfigArray_getPtr(ec, 1, x, y, z)   ! x,y,z pointer to ec positions
      call eloc(0, ec, 'none')

      expU = exp(elU(1))
      vpot = elVen(1) + elVee(1)
      jas = elU(1)
      phi = elPhi(1)
      psiv = elPhi(1) * expU
      do i = 1, ne
         grad(3 * i - 2) = elxDrift(i, 1)        ! drift = grad / psi
         grad(3 * i - 1) = elyDrift(i, 1)
         grad(3 * i) = elzDrift(i, 1)
         jgrad(3 * i - 2) = elUgrad(3 * i - 2, 1)
         jgrad(3 * i - 1) = elUgrad(3 * i - 1, 1)
         jgrad(3 * i) = elUgrad(3 * i, 1)
         fgrad(3 * i - 2) = elFgrad(3 * i - 2, 1) / phi
         fgrad(3 * i - 1) = elFgrad(3 * i - 1, 1) / phi
         fgrad(3 * i) = elFgrad(3 * i, 1) / phi
         lapli(i) = -2 * elEkini(i, 1)          ! lapl / psi
         flapli(i) = elFLapli(i, 1) / phi
         ulapli(i) = elULapli(i, 1)
      end do


      !  analytical values
      eloc0 = elEloc(1)
      jlapl0 = elUlapl(1)
      flapl0 = elFlapl(1) / phi
      vpots0 = elVen(1) + elVee(1) + vpot0
      lapl0 = -2 * (elEloc(1) - vpots0)

      !  Calculate points for numerical derivaties

      nMax = 6
      allocate(t0grad(0 : nMax), t0gradj(0 : nMax), t0gradf(0 : nMax))
      allocate(t0lapl(0 : nMax), t0laplu(0 : nMax), t0laplf(0 : nMax))
      maxError = 0
      maxUError = 0
      maxFError = 0
      lapl = 0
      ulapl = 0
      flapl = 0
      total_lapl = 0
      total_ulapl = 0
      total_flapl = 0

      if (logmode >= 3) then
         write(iul, *)
         write(iul, *) ' numerical derivatives:'
         write(iul, *) ' computed, numerical, abs error, rel error, extrap error, char scale, loss'
      end if

      do k = 1, 3 * ne

         h = h0
         call internal_getElocValues(k, h, psi0, psij, psif)
         loss = lossOfDigits(psi0(1), psi0(3))
         loss0 = loss
         if (loss > 2._r8) then
            h = 10 * h
            call internal_getElocValues(k, h, psi0, psij, psif)
            loss = lossOfDigits(psi0(1), psi0(3))
         end if
         t0grad(0) = numFirstDerivative3Point(psi0(:), h)
         t0lapl(0) = numSecondDerivative3Point(psi0(:), h)
         t0gradj(0) = numFirstDerivative3Point(psij(:), h)
         t0laplu(0) = numSecondDerivative3Point(psij(:), h)
         t0gradf(0) = numFirstDerivative3Point(psif(:), h)
         t0laplf(0) = numSecondDerivative3Point(psif(:), h)
         characteristicScale = sqrt( abs( 1._r8 / t0lapl(0) ))
         characteristicScaleJ = sqrt( abs( jas / t0laplu(0) ))
         characteristicScaleF = sqrt( abs( 1._r8 / t0laplf(0) ))
         do n = 1, nMax
            h = h / 2
            call internal_getElocValues(k, h, psi0, psij, psif)
            t0grad(n) = numFirstDerivative3Point(psi0(:), h)
            t0lapl(n) = numSecondDerivative3Point(psi0(:), h)
            t0gradj(n) = numFirstDerivative3Point(psij(:), h)
            t0laplu(n) = numSecondDerivative3Point(psij(:), h)
            t0gradf(n) = numFirstDerivative3Point(psif(:), h)
            t0laplf(n) = numSecondDerivative3Point(psif(:), h)
         end do

         call nevilleExtrapolation(t0grad(0:), gradient, epsGradient, iter)
         call nevilleExtrapolation(t0gradj(0:), ugradient, epsUGradient, iteru)
         call nevilleExtrapolation(t0gradf(0:), fgradient, epsFGradient, iterf)
         call nevilleExtrapolation(t0lapl(0:), lapl1, epsLapl, iterl)
         call nevilleExtrapolation(t0laplu(0:), ulapl1, epsULapl, iterlu)
         call nevilleExtrapolation(t0laplf(0:), flapl1, epsFLapl, iterlf)

         maxError = max( maxError, abs((gradient - grad(k))/grad(k)) )
         maxUError = max( maxUError, abs((ugradient - jgrad(k))/jgrad(k)) )
         maxFError = max( maxFError, abs((fgradient - fgrad(k))/fgrad(k)) )

         lapl = lapl + lapl1
         ulapl = ulapl + ulapl1
         flapl = flapl + flapl1

         if (logmode >= 3) then
            write(iul, '(i3,4g19.10,3g10.3,i3)') k, grad(k), gradient, gradient - grad(k), &
                    abs((gradient - grad(k))/grad(k)), epsGradient, characteristicScale, loss0, iter
            write(iul, '(i3,4g19.10,2g10.3,i13)') k, jgrad(k), ugradient, ugradient - jgrad(k), &
                    abs((ugradient - jgrad(k))/jgrad(k)), epsUGradient, characteristicScaleJ, iteru
            write(iul, '(i3,4g19.10,2g10.3,i13)') k, fgrad(k), fgradient, fgradient - fgrad(k),  &
                    abs((fgradient - fgrad(k))/fgrad(k)), epsFGradient, characteristicScaleF, iterf
         end if

         if (mod(k, 3) == 0) then
            i = k / 3
            total_lapl = total_lapl + lapl
            total_ulapl = total_ulapl + ulapl
            total_flapl = total_flapl + flapl

            if (logmode >= 3) then
               write(iul, '(i3,4g19.10,30x,i3)') i, lapli(i), lapl, lapli(i) - lapl, &
                       (lapli(i) - lapl) / lapli(i), iterl
               write(iul, '(i3,4g19.10,30x,i3)') i, ulapli(i), ulapl, ulapli(i) - ulapl, &
                       (ulapli(i) - ulapl) / jlapl0, iterlu
               write(iul, '(i3,4g19.10,30x,i3)') i, flapli(i), flapl, flapli(i) - flapl, &
                       (flapli(i) - flapl) / flapli(i), iterlf
            end if

            lapl = 0
            ulapl = 0
            flapl = 0
         end if

      end do

      numEloc = -0.5d0 * total_lapl + vpots0

      if (maxError > tolGrad .or. abs(eloc0 - numEloc) > tolEloc) passed = .false.

      if (logmode >= 3) then
         write(iul, '(a,5g20.10)') ' laplacian(Psi):', lapl0, total_lapl, lapl0 - total_lapl, &
                 (lapl0 - total_lapl) / lapl0
         write(iul, '(a,4g20.10)') ' Eloc:', eloc0, numEloc, eloc0 - numEloc, &
                 (numEloc - eloc0) / eloc0
         write(iul, '(a,5g20.10)') ' laplacian(U):', jlapl0, total_ulapl, jlapl0 - total_ulapl, &
                 (jlapl0 - total_ulapl) / jlapl0
         write(iul, '(a,5g20.10)') ' laplacian(Phi):', flapl0, total_flapl, flapl0 - total_flapl, &
                 (flapl0 - total_flapl) / flapl0
         write(iul, '(a,6g12.3)') ' grad / lapl errors:', maxError, maxUError, maxFError, &
                 (lapl0 - total_lapl) / lapl0, (jlapl0 - total_ulapl) / jlapl0, (flapl0 - total_flapl) / flapl0
      end if

      if (logmode >= 2) then
         write(iul, *)
         if (passed) then
            write(iul, *) ' numerical gradient and eloc test (Neville) *** passed *** '
         else
            write(iul, *) ' numerical gradient and eloc test (Neville) *** failed ***'
            write(iul, *) maxError, tolGrad, abs(eloc0 - numEloc), tolEloc
            write(iul, *) maxUError, maxFError
         end if
      end if


      k = 3
      write(iul,*) z(1)
      do i = 1, 10
         h = i * epsilon(1._r8)
         call internal_getElocValues(k, h, psi0, psij, psif)
         write(iul,*) i, z(1)
         do j = 1, 3
            write(iul,'(4g25.16)') h, psi0(j), psif(j), psif(j)
         end do
      end do

      call eConfigArray_destroy(ec)

   contains


      subroutine internal_getElocValues(k, h, psi0, psij, psif)
         integer, intent(in) :: k
         real(r8), intent(in) :: h
         real(r8), intent(out) :: psi0(:), psij(:), psif(:)
         integer kk, krmd

         krmd = mod(k + 2, 3)
         kk = (k + 2) / 3
         select case (krmd)
         case (0)  ! x
            x(kk) = x(kk) + h
         case (1)  ! y
            y(kk) = y(kk) + h
         case (2)  ! z
            z(kk) = z(kk) + h
         end select
         call eloc(0, ec, 'none')
         psi0(3) = elPhi(1) / psiv * exp(elU(1))
         psij(3) = elU(1)
         psif(3) = elPhi(1) / phi

         psi0(2) = 1.0_r8
         psij(2) = jas
         psif(2) = 1.0_r8

         select case (krmd)
         case (0)  ! x
            x(kk) = x(kk) - 2._r8 * h
         case (1)  ! y
            y(kk) = y(kk) - 2._r8 * h
         case (2)  ! z
            z(kk) = z(kk) - 2._r8 * h
         end select
         call eloc(0, ec, 'none')
         psi0(1) = elPhi(1) / psiv * exp(elU(1))
         psij(1) = elU(1)
         psif(1) = elPhi(1) / phi

         select case (krmd)
         case (0)  ! x
            x(kk) = x(kk) + h
         case (1)  ! y
            y(kk) = y(kk) + h
         case (2)  ! z
            z(kk) = z(kk) + h
         end select

      end subroutine internal_getElocValues


   end subroutine elocDerivExtrap1Test


   subroutine updatetest(xx, yy, zz)

      real(r8), intent(inout) :: xx(:), yy(:), zz(:)
      real(r8), parameter :: dt = 0.01    ! time step for move
      real(r8) :: x1(ne), y1(ne), z1(ne)
      real(r8), pointer :: x(:), y(:), z(:)
      real(r8) phi, u
      integer ie
      type(eConfigArray) :: ec

      call eConfigArray_new(ec, ne, 1)
      call eConfigArray_set(ec, 1, xx, yy, zz)
      call eConfigArray_getPtr(ec, 1, x, y, z)   ! x,y,z pointer to ec positions

      write(iul, *)
      write(iul, *) ' test of psi update with psit:: currently disabled'
      call abortp("UPDATE CURRENTLY NOT POSSIBLE")

      call eloc(0, ec, 'none')
      write(iul, '(2(A8,G14.8))') ' Phi = ', elPhi(1), ' U = ', elU(1)

      x1 = elxDrift(1:ne, 1) * dt
      y1 = elyDrift(1:ne, 1) * dt
      z1 = elzDrift(1:ne, 1) * dt

      x(1) = x(1) + x1(1)
      y(1) = y(1) + y1(1)
      z(1) = y(1) + z1(1)

      ie = 1
      !     !!!!!!loc_method = "SD*J"

      !     !!!!!!call psit(.true.,ie,x,y,z,phi,u)

      write(iul, *) ' update with psit'
      write(iul, '(2(A8,G14.8))') ' Phi = ', phi, ' U = ', u

      call eloc(0, ec, 'none')
      write(iul, *) ' check with eloc'
      write(iul, '(2(A8,G14.8))') ' Phi = ', elPhi(1), ' U = ', elU(1)
      call assertEqualRelative(phi, elPhi(1), 1.d-5, &
             '(updatetest): phi not equal')
      call assertEqualRelative(u, elU(1), 1.d-5, &
             '(updatetest): u not equal')

      do ie = 2, ne
         x(ie) = x(ie) + x1(ie)
         y(ie) = y(ie) + y1(ie)
         z(ie) = y(ie) + z1(ie)
         !        !!!!!!!!!!!!!call psit(.false.,ie,x,y,z,phi,u)
         write(iul, *) ' update with psit ie=', ie
         write(iul, '(2(A8,G14.8))') ' Phi = ', phi, ' U = ', u
      enddo

      call eloc(0, ec, 'none')
      write(iul, *) ' check with eloc'
      write(iul, '(2(A8,G14.8))') ' Phi = ', elPhi(1), ' U = ', elU(1)
      call assertEqualRelative(phi, elPhi(1), 1.d-5, &
             '(updatetest): phi not equal')
      call assertEqualRelative(u, elU(1), 1.d-5, &
             '(updatetest): u not equal')

      write(iul, *) ' PASSED: eloctest -- updatetest'

      call eConfigArray_destroy(ec)

   end subroutine


   function lossOfDigits(y1, y2) result(res)
      real(r8), intent(in) :: y1
      real(r8), intent(in) :: y2
      real(r8) :: res
      res = 0._r8
      if ( abs(y1 - y2) < abs(y1) ) then
         res = log10( abs( y1 / (y1 - y2) ) )
      end if

   end function lossOfDigits


END MODULE elocTest_m
