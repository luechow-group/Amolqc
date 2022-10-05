! Copyright (C) 1996, 2002, 2005, 2007-2008, 2015-2016, 2018 Arne Luechow
! Copyright (C) 2003-2005 Christian Diedrich
! Copyright (C) 2005 Annika Bande
! Copyright (C) 2013, 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later


! f90 module eloc-- local energy calculation

MODULE eloc_m
    use kinds_m, only : r8
    use OMP_LIB
    use wfData_m
    use waveFunction_m, only : readwriteWF
    use ecp_m, only: EcpType
    use wfParameters_m
    use jastrowParamData_m
    use elocData_m
    use eConfigs_m, only: eConfigArray, eConfigArray_get, eConfigArray_size, eConfigArray_set
    use jastrow_m, only : jasCalcWDerivs, JasCalc
    use aos_m, only : aocalc, aosplcalc
    use mos_m, only : mocalc
    use multiDet_m, only : mdetcalc, MDET_NONE, MDET_LU_ERR, MDET_INV_ERR
    use multiDetParam_m
    use moParam_m
    use ecp_m
    use ecpIo_m, only : initecpparams
    use aoMo_m
    use rdataUpdate_m
    use coulombDensity_m, only: get_density_repulsion, is_coulomb_density_initialized

    implicit none

    type(EcpType) :: ecp

    integer, parameter :: ELOC_NONE=0, ELOC_X_INF_ERR=1, ELOC_LU_ERR=2, ELOC_INV_ERR=3

    private
    public :: wf_init, wf_ecp_init, wf_getNCoreElecs, eloc, vpotcalc, elCutOff, elocUpdateDistAndPot, elocCalcPhi, &
            resetEloc
    public :: ELOC_NONE, ELOC_X_INF_ERR, ELOC_LU_ERR, ELOC_INV_ERR


CONTAINS

    !     ! first new routines for wf objects, currently only one, hidden from QMC
    !     ! stored here


    subroutine wf_init(lines, nl)
        integer, intent(in) :: nl
        character(len = *), intent(in) :: lines(nl)

        call readwriteWF(lines, nl, ecp)
    end subroutine wf_init


    subroutine wf_ecp_init(lines, nl)
        integer, intent(in) :: nl
        character(len = *), intent(in) :: lines(nl)
        call initecpparams(lines, nl, ecp)
    end subroutine wf_ecp_init


    function wf_getNCoreElecs(a) result(res)
        integer :: res
        integer, optional, intent(in) :: a
        if (present(a)) then
            res = ecp%NCoreElecs(atom = a)
        else
            res = ecp%NCoreElecs()
        endif
    end function wf_getNCoreElecs


    subroutine eloc(ie, ec, optType, error_code, wfpDef, wfpDT, twoLevelStep, doCalc, &
            tmove, tau, isMoved)

        ! TODO:
        ! free format!
        ! remove optType (is part of wfpDef)
        ! check two-level propagator (requires saving intermediate results)
        ! remove EConfigArray! replace by EConfig




        ! ELOC calculates psi, grad(psi), kin. energy ekin,
        ! and pot. energy Vpot, for given electron configurations (x,y,z) in ec
        ! one-electron update for ie>0. Results are stored in module
        ! data structure elocdata (eloc_d.f90) for fast and flexible access.

        ! Modifications/Versions:
        ! This version allows efficient calculation of several electron
        ! configurations simultaneously. ec contains several configurations
        ! (max is mMaxElecConfigs), AL, 2011
        ! Version 1.5 (28.4.98) "module" version. No reference to walker, Just
        !                       evaluate E_local and other stuff (stored in
        !                       local data structure in eloc.h) for some position
        !                       vector.
        ! Version 1.4 (9.3.98)  include spline evaluation of radial part of AO's
        ! Version 1.3 (sept 97) module version: use aos_m,mos,jastrow,mdet modules
        ! Version 1.2 (5/31/96) include ECP's
        ! Version 1.1 (5/20/96) allow single-electron moves; calc local energy
        !                       for the orbital part only (COMMON eltest)
        ! Version 1.0 (1/29/96)
        integer, intent(in) :: ie       ! 0: update all electron; i: update only electron i
        type(eConfigArray), intent(inout) :: ec       ! contains electron configurations: T move on exit
        character(len = *), intent(in) :: optType  ! transferred to call to JasWDerivs (REMOVE! REPLACED BY WFPDEF)
        integer, intent(out), optional :: error_code 
        type(WFParamDef), intent(in), optional :: wfpDef   ! definition of wf parameters; triggers calculation of
        !                                                                  !                                    param derivatives
        type(WFParamDerivTerms), intent(inout), optional :: wfpDT    ! contains on output parameter derivatives (of Psi and Elocal)
        integer, intent(in), optional :: twoLevelStep     ! two Level propagator step
        logical, intent(in), optional :: doCalc(eConfigArray_size(ec)) ! for twoLevel
        integer, intent(in), optional :: tmove    ! triggers T move calculation
        real(r8), intent(in), optional :: tau      ! time step for T moves
        logical, intent(inout), optional :: isMoved(:) ! T move accepted? for each EConfig

        integer i, ii, iull, n, ierr, nElecConfigs
        real(r8) ekin
        real(r8) wtimerPhi(4)
#ifdef WTIMER
        real(r8) wtimer_i(8), wtimer_f(8)
#endif
        !     ! automatic arrays for actual electron number!
        real(r8) x(ne), y(ne), z(ne)
        real(r8) rai(ncenter, ne), rij(ne, ne)
        real(r8) rrai(ncenter, ne, eConfigArray_size(ec))    ! rai for all elec configs
        real(r8) rrij(ne, ne, eConfigArray_size(ec))         ! rij for all elec configs
        real(r8) uu(eConfigArray_size(ec))
        real(r8) uugrad(3 * ne, eConfigArray_size(ec))
        real(r8) uulapl(eConfigArray_size(ec))
        real(r8) uulapli(ne, eConfigArray_size(ec))
        real(r8) phi(eConfigArray_size(ec)), &
                fgrad(3 * ne, eConfigArray_size(ec)), &
                flapl(eConfigArray_size(ec)), &
                flapli(ne, eConfigArray_size(ec))
        real(r8) ecpLocal, ecpNonlocal
        real(r8), allocatable :: ecpNLocalk(:)
        type(TMoveDataType) :: TMoveData


#ifdef WTIMER
        if (MASTER) wtimer_i(WTIMER_ELOC) = omp_get_wtime()
#endif

        if (present(error_code)) error_code = ELOC_NONE
        nElecConfigs = eConfigArray_size(ec)

        if (asserts) call assert(nElecConfigs==1, "eloc: many EConfigs disabled")
        if (debug) then
            call eConfigArray_get(ec, 1, x, y, z)
            ! check input data for NaN or Inf
            if (.not.(all(abs(x)<huge(1._r8)) .and. all(abs(y)<huge(1._r8)) .and. all(abs(z)<huge(1._r8)))) then
                if (present(error_code)) then
                    error_code = ELOC_X_INF_ERR
                    return
                end if
            end if
        end if

        iull = mytid + 80  ! for write output for each node

        rrai = 0; rrij = 0
        rai = 0; rij = 0

        if (present(twoLevelStep)) then
            call elocUpdateDistAndPot(ec, ie, nElecConfigs, rai, rij, rrai, rrij)
            if (twoLevelStep == 1) then
                call assert(separated_electrons == 0, 'eloc: twoLevelStep NOT IMPLEMENTED with separated_electrons')
#ifdef WTIMER
                if (MASTER) wtimer_i(WTIMER_PHI) = omp_get_wtime()
#endif

                call resetEloc()
                call elocCalcPhi(ie, nElecconfigs, rrai, ec, wtimerPhi, ierr)
                if (present(error_code)) then
                    if (ierr /= ELOC_NONE) then
                        error_code = ierr
                        return
                    end if
                end if

#ifdef WTIMER
                if (MASTER) wtimer_f(WTIMER_PHI) = omp_get_wtime()
#endif

            else

#ifdef WTIMER
                if (MASTER) wtimer_i(WTIMER_JAS) = omp_get_wtime()
#endif

                call elocCalcJas(ie, ec, optType, nElecConfigs, uu, uugrad, uulapl, uulapli, rrai, rrij, twoLevelStep, doCalc)

#ifdef WTIMER
                if (MASTER) wtimer_f(WTIMER_JAS) = omp_get_wtime()
#endif
                if (twoLevelStep == 2) return

            endif

        else

            call resetEloc()
            call elocUpdateDistAndPot(ec, ie, nElecConfigs, rai, rij, rrai, rrij)

#ifdef WTIMER
            if (MASTER) wtimer_i(WTIMER_PHI) = omp_get_wtime()
#endif

            call elocCalcPhi(ie, nElecconfigs, rrai, ec, wtimerPhi, ierr)
            if (present(error_code)) then
                if (ierr /= ELOC_NONE) then
                    error_code = ierr
                    return
                end if
            end if

            ! this line has to be after elocUpdateDistandPot and after elocCalcPhi
            if (separated_electrons /= 0) then
                call assert(ie == 0, 'eloc: one electron update NOT IMPLEMENTED with separated electrons')
                call assert(is_coulomb_density_initialized(), &
                        'eloc: $coulomb_density has to be called before qmc with separated electrons')
                call add_Vee_density_repulsion(nElecConfigs, ec)
            end if

#ifdef WTIMER
            if (MASTER) wtimer_f(WTIMER_PHI) = omp_get_wtime()
            if (MASTER) wtimer_i(WTIMER_JAS) = omp_get_wtime()
#endif

            call elocCalcJas(ie, ec, optType, nElecConfigs, uu, uugrad, uulapl, uulapli, rrai, rrij)

#ifdef WTIMER
            if (MASTER) wtimer_f(WTIMER_JAS) = omp_get_wtime()
#endif
        endif

        if (ecp%isInitialised()) then

#ifdef WTIMER
            if (MASTER) wtimer_i(WTIMER_PP) = omp_get_wtime()
#endif
            if (present(wfpDef) .and. present(wfpDT)) then
                if (wfpDT%ELiCalc .and. (wfpDT%noCalc.eqv. .false.))then
                    call elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal, &
                            ecpNonlocal, ecpNLocalk = ecpNLocalk, wfpDef = wfpDef)
                else
                    call elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal, &
                            ecpNonlocal)
                end if
            else if (present(tmove)) then
                TMoveData%tmove = tmove
                TMoveData%tau = tau
                call elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal, &
                        ecpNonlocal, TMoveData = TMoveData)
                isMoved(1) = TMoveData%isMoved
            else
                call elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal, ecpNonlocal)
            endif

            elVen(1) = elVen(1) + ecpLocal + ecpNonlocal
            elECPPot = ecpLocal
            elECPPotNl = ecpNonlocal

#ifdef WTIMER
            if (MASTER) wtimer_f(WTIMER_PP) = omp_get_wtime()
#endif
        endif

        n = 1   !!! only one EConfig !!!

        if (present(doCalc)) then
            if (.not. doCalc(1)) return
        endif

        !     // Grad(Psi_G) / psi, Laplacian and kinetic energy
        phi(n) = elPhi(n)
        fgrad(:, n) = elFgrad(:, n)
        flapl(n) = elFlapl(n)
        flapli(:, n) = elFlapli(:, n)
        ii = 1
        do i = 1, ne
            elxDrift(i, n) = fgrad(ii, n) / phi(n) + uugrad(ii, n)
            ii = ii + 1
            elyDrift(i, n) = fgrad(ii, n) / phi(n) + uugrad(ii, n)
            ii = ii + 1
            elzDrift(i, n) = fgrad(ii, n) / phi(n) + uugrad(ii, n)
            ii = ii + 1
        enddo
        !     // kinetic energy for Psi_G
        ii = 1
        do i = 1, ne
            elEkini(i, n) = -0.5_r8 * (flapli(i, n) / phi(n) + uulapli(i, n)&
                    + 2._r8 * fgrad(ii, n) / phi(n) * uugrad(ii, n)&
                    + 2._r8 * fgrad(ii + 1, n) / phi(n) * uugrad(ii + 1, n)&
                    + 2._r8 * fgrad(ii + 2, n) / phi(n) * uugrad(ii + 2, n)&
                    + uugrad(ii, n) * uugrad(ii, n)&
                    + uugrad(ii + 1, n) * uugrad(ii + 1, n)&
                    + uugrad(ii + 2, n) * uugrad(ii + 2, n))
            ii = ii + 3
        enddo

        if(do_epart) elEkin_epart(1:ne, n) = elEkini(1:ne, n)

        ekin = 0._r8                ! local kinetic energy with Jastrow
        do i = 1, ne
            ekin = ekin + elEkini(i, n)
        enddo

        !      // local energies
        elEloc(n) = ekin + elVee(n) + elVen(n) + vpot0

        if (logmode >= 5) then
            write(iul, '(A,5G20.10)') "eloc:", elEloc(n), ekin, elVen(n) + elVee(n)+ vpot0, elECPPot, elECPPotNl
            write(iul, '(A,2G20.10)') "phi,U:", elPhi(n), uu
        endif

        if (mElocCutOff) then
            call elCutOff(n)
            if (logmode >= 5) write(iul, '(A,4G15.6)') "eloc (after cutoff):", elEloc(n), ekin, &
                    elVen(n) + elVee(n) + vpot0, elECPPot
        endif

        if (present(wfpDef) .and. present(wfpDT)) then
            wfpDT%eloc = elEloc(1)
            wfpDT%phi = phi(1)
            wfpDT%U = uu(1)
            !        ! calculate contributions to wf parameter derivative terms
            if(wfpDT%noCalc .eqv. .false.) call internal_calcDerivContribs()
        endif



#ifdef WTIMER
        if (MASTER) then
            wtimer_f(WTIMER_ELOC) = omp_get_wtime()
            wtimer(WTIMER_ELOC) = wtimer(WTIMER_ELOC) + (wtimer_f(WTIMER_ELOC) - wtimer_i(WTIMER_ELOC))
            wtimer(WTIMER_PHI) = wtimer(WTIMER_PHI) + (wtimer_f(WTIMER_PHI) - wtimer_i(WTIMER_PHI))
            wtimer(WTIMER_PP) = wtimer(WTIMER_PP) + (wtimer_f(WTIMER_PP) - wtimer_i(WTIMER_PP))
            wtimer(WTIMER_JAS) = wtimer(WTIMER_JAS) + (wtimer_f(WTIMER_JAS) - wtimer_i(WTIMER_JAS))
            wtimer(WTIMER_AO:WTIMER_MDET) = wtimer(WTIMER_AO:WTIMER_MDET) + wtimerPhi
        endif
#endif


    contains

        subroutine internal_calcDerivContribs()

            integer np, npJ, npMO, npCI, k, l, i
            real(r8) tmp

            call assert(separated_electrons == 0,&
                    'internal_calcDerivContribs: eloc derivs NOT IMPLEMENTED with separated_electrons')

            call wfparams_getNJastrowParams(wfpDef, npJ)
            call wfparams_getNDetParams(wfpDef, npCI, npMO)
            np = npJ + npCI + npMO

            if (.not.allocated(wfpDT%Eli)) then
                allocate(wfpDT%fi(np), wfpDT%Eli(np), wfpDT%fij(np, np))
            else
                if (size(wfpDT%Eli) /= np) call abortp("eloc:internal_calcDerivContribs: illegal sizes")
            endif

            !          wfpDT%eloc = elEloc(1)
            !          wfpDT%phi  = phi(1)
            !          wfpDT%U    = uu(1)

            !        ! parameter order: Jastrow, MO, CI

            !        ! Jastrow block
            if (npJ > 0) then

                if (wfpDT%fiCalc) then
                    do k = 1, npJ
                        wfpDT%fi(k) = uk(k)
                    end do
                end if

                if (wfpDT%ELiCalc) then
                    do k = 1, npJ
                        tmp = 0._r8
                        do i = 1, ne
                            tmp = tmp + elxDrift(i, 1) * ukgrad(3 * i - 2, k) + elyDrift(i, 1) * ukgrad(3 * i - 1, k)&
                                    + elzDrift(i, 1) * ukgrad(3 * i, k)
                        end do
                        wfpDT%ELi(k) = -0.5_r8 * uklapl(k) - tmp
                    end do
                    if (ecp%isInitialised()) then
                        wfpDT%ELi(1:npJ) = wfpDT%ELi(1:npJ)&
                                + ecpNLocalk(1:npJ) - ecpNonlocal * uk(1:npJ)
                    end if
                end if

                if (wfpDT%fijCalc) then
                    do k = 1, npJ
                        do l = 1, k
                            wfpDT%fij(k, l) = uk(l) * uk(k)
                        end do
                    end do
                end if

            endif

            !        ! MO block
            if (npMO > 0) then
                call moparam_calcderivs(wfpDT%ELiCalc)

                if (wfpDT%fiCalc) then
                    do i = 1, npMO
                        wfpDT%fi(npJ + i) = mok(i) / elPhi(1)
                    end do
                end if

                if (wfpDT%ELiCalc) then
                    do i = 1, npMO
                        tmp = moklapl(i) + 2 * dot_product(mokgrad(:, i), elUGrad(1:3 * ne, 1))&
                                + mok(i) * dot_product(elUGrad(1:3 * ne, 1), elUGrad(1:3 * ne, 1)) + mok(i) * elULapl(1)
                        tmp = tmp * (-0.5_r8) / elPhi(1)
                        wfpDT%Eli(npJ + i) = tmp + (elVee(1) + elVen(1) + vpot0 - elEloc(1)) * mok(i) / elPhi(1)
                    end do
                    if (ecp%isInitialised()) then
                        wfpDT%ELi(npJ + 1:npJ + npMO) = wfpDT%ELi(npJ + 1:npJ + npMO) + &
                                ecpNLocalk(npJ + 1:npJ + npMO) - ecpNonlocal * mok(1:npMO) / elPhi(1)
                    end if
                end if

                if (wfpDT%fijCalc) then
                    wfpDT%fij(npJ + 1:npJ + npMO, npJ + 1:npJ + npMO) = 0._r8
                end if

            endif

            !        ! CI block
            if (npCI > 0) then
                call mdetcalcderivs()

                if (wfpDT%fiCalc) then
                    do i = 1, npCI
                        wfpDT%fi(npMO + npJ + i) = fk(i) / elPhi(1)
                    end do
                end if

                if (wfpDT%ELiCalc) then
                    do i = 1, npCI
                        tmp = fklapl(i) + 2 * dot_product(fkgrad(:, i), elUGrad(1:3 * ne, 1))&
                                + fk(i) * dot_product(elUGrad(1:3 * ne, 1), elUGrad(1:3 * ne, 1)) + fk(i) * elULapl(1)
                        tmp = tmp * (-0.5_r8) / elPhi(1)
                        wfpDT%ELi(npMO + npJ + i) = tmp + (elVee(1) + elVen(1) + vpot0 - elEloc(1)) * fk(i) / elPhi(1)
                    end do
                    if (ecp%isInitialised()) then
                        wfpDT%ELi(npJ + npMO + 1:np) = wfpDT%ELi(npJ + npMO + 1:np) + &
                                ecpNLocalk(npJ + npMO + 1:np) - ecpNonlocal * fk(1:npCI) / elPhi(1)
                    end if
                end if

                if (wfpDT%fijCalc) then
                    wfpDT%fij(npMO + npJ + 1:np, npMO + npJ + 1:np) = 0._r8
                end if
            end if


            !        ! mixed blocks
            if (wfpDT%fijCalc) then
                !           ! jas-MO
                do k = 1, npJ
                    do i = 1, npMO
                        wfpDT%fij(npJ + i, k) = mok(i) / elPhi(1) * uk(k)
                    end do
                end do
                !           ! jas-CI
                do k = 1, npJ
                    do i = 1, npCI
                        wfpDT%fij(npJ + npMO + i, k) = fk(i) / elPhi(1) * uk(k)
                    end do
                end do
                !           ! MO-CI
                do k = 1, npMO
                    do i = 1, npCI
                        wfpDT%fij(npJ + npMO + i, npJ + k) = fk(i) / elPhi(1) * mok(k)
                    end do
                end do
            end if

        end subroutine internal_calcDerivContribs
    end subroutine eloc


    !===========================================================

    !     --------------------------!
    subroutine resetEloc()
        !     --------------------------!
        elPhi = 0
        elFgrad = 0
        elFlapl = 0
        elU = 0
        elUgrad = 0
        elUlapl = 0
        elxDrift = 0; elyDrift = 0; elzDrift = 0
        elECPPot = 0
        elECPPotNl = 0
        elVen = 0
        elVee = 0
        elEkini = 0
        elEloc = 0
        elSingularity = .false.
    end subroutine resetEloc

    !===========================================================

    !     -----------------------------------------!
    subroutine vpotcalc(ie, x, y, z, rai, rij, ven, vee, n)
        !     -----------------------------------------!

        ! vpotcalc calculates the potential energies for given position vector
        ! NOTE: rij is calculated ONLY FOR i<j
        ! ie == 0: new calculation of vpot,rai and rij
        ! ie == n: update only for n-th electron
        ! on entry: for one-electron update the vectors x,y,z must be changed compared
        !           to the last call to this routine only at electron ie.
        !           rai and rij contain the old distances.
        !           vpot is then required have 'old' value of vpot

        integer, intent(in), optional :: n
        integer, intent(in) :: ie    ! ie==0: recalculate rai,rij,vpot,
        !                                    !        else update electron ie
        real(r8), intent(in) :: x(:), y(:), z(:)
        real(r8), intent(inout) :: rai(:, :), rij(:, :)
        real(r8), intent(inout) :: ven  ! input: for ie>0 old potential
        real(r8), intent(inout) :: vee  ! input: for ie>0 old potential
        !                                      ! output: vpot for x,y,z
        integer a, i, j
        real(r8) veni(nmax)             ! individual potential energy
        save veni

        !     // Potential energy vpot and distances
        if (do_epart) then
            call assert(separated_electrons == 0,&
                    'vpotcalc: epart NOT IMPLEMENTED with separated_electrons')
            if (ie == 0) then
                !        // Recalculate all distances and potential energy
                ven = 0._r8
                vee = 0._r8
                do i = 1, ne
                    veni(i) = 0._r8
                    do a = 1, ncenter
                        rai(a, i) = sqrt((x(i) - atoms(a)%cx)**2&
                                + (y(i) - atoms(a)%cy)**2 + (z(i) - atoms(a)%cz)**2)
                        veni(i) = veni(i) - atoms(a)%za / rai(a, i)
                        elVne_epart(a, i, n) = -atoms(a)%za / rai(a, i)
                    enddo
                    ven = ven + veni(i)

                    do j = i + 1, ne
                        rij(i, j) = sqrt((x(i) - x(j))**2 + &
                                (y(i) - y(j))**2 + (z(i) - z(j))**2)
                        vee = vee + 1._r8 / rij(i, j)
                        elVee_epart(i, j, n) = 1._r8 / rij(i, j)
                    enddo
                enddo
            else
                !        // Update distances and potential energy for electron ie
                i = ie
                ven = ven - veni(i)
                veni(i) = 0._r8
                do a = 1, ncenter
                    rai(a, i) = sqrt((x(i) - atoms(a)%cx)**2&
                            + (y(i) - atoms(a)%cy)**2 + (z(i) - atoms(a)%cz)**2)
                    veni(i) = veni(i) - atoms(a)%za / rai(a, i)
                    elVne_epart(a, i, n) = -atoms(a)%za / rai(a, i)
                enddo
                ven = ven + veni(i)

                do j = 1, i - 1
                    vee = vee - 1._r8 / rij(j, i)
                    rij(j, i) = sqrt((x(i) - x(j))**2 + &
                            (y(i) - y(j))**2 + (z(i) - z(j))**2)
                    vee = vee + 1._r8 / rij(j, i)
                    elVee_epart(j, i, n) = 1._r8 / rij(j, i) !is this correct?
                enddo
                do j = i + 1, ne
                    vee = vee - 1._r8 / rij(i, j)
                    rij(i, j) = sqrt((x(i) - x(j))**2 + &
                            (y(i) - y(j))**2 + (z(i) - z(j))**2)
                    vee = vee + 1._r8 / rij(i, j)
                    elVee_epart(i, j, n) = 1._r8 / rij(i, j) !is this correct?
                enddo
            endif !! (ie==0)
        else  !! (.not.do_epart)
            if (ie == 0) then
                !        // Recalculate all distances and potential energy
                ven = 0._r8
                vee = 0._r8
                do i = 1, ne
                    veni(i) = 0._r8
                    do a = 1, ncenter
                        rai(a, i) = sqrt((x(i) - atoms(a)%cx)**2&
                                + (y(i) - atoms(a)%cy)**2 + (z(i) - atoms(a)%cz)**2)
                        veni(i) = veni(i) - atoms(a)%za / rai(a, i)
                    enddo
                    ven = ven + veni(i)

                    do j = i + 1, ne
                        rij(i, j) = sqrt((x(i) - x(j))**2 + &
                                (y(i) - y(j))**2 + (z(i) - z(j))**2)
                        vee = vee + 1._r8 / rij(i, j)
                    enddo
                enddo
            else
                call assert(separated_electrons == 0,&
                        'vpotcalc: 1e update NOT IMPLEMENTED with separated_electrons')
                !        // Update distances and potential energy for electron ie
                i = ie
                ven = ven - veni(i)
                veni(i) = 0._r8
                do a = 1, ncenter
                    rai(a, i) = sqrt((x(i) - atoms(a)%cx)**2&
                            + (y(i) - atoms(a)%cy)**2 + (z(i) - atoms(a)%cz)**2)
                    veni(i) = veni(i) - atoms(a)%za / rai(a, i)
                enddo
                ven = ven + veni(i)

                do j = 1, i - 1
                    vee = vee - 1._r8 / rij(j, i)
                    rij(j, i) = sqrt((x(i) - x(j))**2 + &
                            (y(i) - y(j))**2 + (z(i) - z(j))**2)
                    vee = vee + 1._r8 / rij(j, i)
                enddo
                do j = i + 1, ne
                    vee = vee - 1._r8 / rij(i, j)
                    rij(i, j) = sqrt((x(i) - x(j))**2 + &
                            (y(i) - y(j))**2 + (z(i) - z(j))**2)
                    vee = vee + 1._r8 / rij(i, j)
                enddo
            endif !! (ie==0)
        endif !! (do_epart)

    end subroutine vpotcalc


    !     ------------------------
    subroutine elCutOff(n)
        !     ------------------------

        !  Cutoff of Rothstein/Vrbik (see Umrigar 1993 paper)
        !  cutoff for both local energy and drift. Cutoff vanishes
        !  for tau->0. Reference for local energy cutoff is 'evar'.
        !  Careful: 'evar' must have an appropriate value
        !  Alternatively: use best current estimate for <E>
        !  Now: CUTOFF for E_local(Psi_T) only, not for Psi_G!
        !  but: for Drift(Psi_G)

        integer, intent(in) :: n    ! electron configuration
        integer i
        real(r8) frac
        logical driftCut

        frac = mCutOffFactor / sqrt(mTauCutOff)

        if (abs(elEloc(n) - E_trial) > 2._r8 * frac) then
            if (elEloc(n) >= E_trial) then
                elEloc = E_trial + 2._r8 * frac
            else
                elEloc = E_trial - 2._r8 * frac
            endif
            mElocCut = mElocCut + 1
        endif
        mElocCutCount = mElocCutCount + 1

        frac = mCutOffFactor / mTauCutOff
        driftCut = .false.
        do i = 1, ne
            if (abs(elxDrift(i, n)) > frac) then
                driftCut = .true.
                if (elxDrift(i, n) >= 0._r8) then
                    elxDrift(i, n) = frac
                else
                    elxDrift(i, n) = -frac
                endif
            endif
            if (abs(elyDrift(i, n)) > frac) then
                driftCut = .true.
                if (elyDrift(i, n) >= 0._r8) then
                    elyDrift(i, n) = frac
                else
                    elyDrift(i, n) = -frac
                endif
            endif
            if (abs(elzDrift(i, n)) > frac) then
                driftCut = .true.
                if (elzDrift(i, n) >= 0._r8) then
                    elzDrift(i, n) = frac
                else
                    elzDrift(i, n) = -frac
                endif
            endif
        enddo
        if (driftCut) then
            mDriftCut = mDriftCut + 1
        end if
        mDriftCutCount = mDriftCutCount + 1

    end subroutine elCutOff


    !     ------------------------
    subroutine elocUpdateDistAndPot(ec, ie, nElecConfigs, rai, rij, rrai, rrij, doCalc)
        !     ------------------------
        integer, intent(in) :: ie       ! 0: update all electron; i: update only electron i
        type(eConfigArray), intent(inout) :: ec  ! contains electron configurations
        real(r8) x(ne), y(ne), z(ne), ven, vee
        integer :: n
        integer, intent(in) :: nElecConfigs
        real(r8), intent(inout) :: rai(ncenter, ne), rij(ne, ne)
        real(r8), intent(inout) :: rrai(ncenter, ne, eConfigArray_size(ec))    ! rai for all elec configs
        real(r8), intent(inout) :: rrij(ne, ne, eConfigArray_size(ec))         ! rij for all elec configs
        logical, intent(in), optional :: doCalc(eConfigArray_size(ec))

        do n = 1, nElecConfigs
            if(present(doCalc)) then
                if(.not. doCalc(n)) cycle
            endif
            call eConfigArray_get(ec, n, x, y, z)
            !        ! check input data for NaN or Inf
            call assert(all(abs(x)<huge(1._r8)) .and.&
                    all(abs(y)<huge(1._r8)) .and. all(abs(z)<huge(1._r8)), &
                    "eloc: illegal x,y,z coords in ec")
            !        // calculate/update particle distances and potential energy
            call vpotcalc(0, x, y, z, rai, rij, ven, vee, n)
            call assert(all(rai<huge(1._r8)), "eloc: illegal rai values")
            call assert(all(rij<huge(1._r8)), "eloc: illegal rij values")
            elVen(n) = ven
            elVee(n) = vee
            rrai(:, :, n) = rai
            rrij(:, :, n) = rij
        end do
        return
    end subroutine elocUpdateDistAndPot


    subroutine elocCalcPhi(ie, nElecconfigs, rrai, ec, wtimerPhi, error_code)
    !------------------------------------------------------------------------
        integer, intent(in)       :: ie       ! 0: update all electron; i: update only electron i
        integer, intent(in)       :: nElecConfigs
        real(r8), intent(in)  :: rrai(:, :, :)    ! rai for all elec configs
        real(r8), intent(out) :: wtimerPhi(4)
        type(eConfigArray), intent(inout) :: ec  ! contains electron configurations
        integer, intent(out)      :: error_code 

        real(r8) x(ne), y(ne), z(ne)
        integer :: n, ierr, i
        real(r8) phi(eConfigArray_size(ec)), &
                fgrad(3 * ne, eConfigArray_size(ec)), &
                flapl(eConfigArray_size(ec)), &
                flapli(ne, eConfigArray_size(ec))
#ifdef WTIMER
        real(r8) wtimer1, wtimer2, wtimer3, wtimer4
#endif

        wtimerPhi = 0._r8
        error_code = ELOC_NONE

        if (aomocomb) then ! use combined AO/MO calculation

            ECLOOP : do n = 1, nElecConfigs

                call eConfigArray_get(ec, n, x, y, z)
                !        !rai = rrai(:,:,n)
#ifdef WTIMER
                if (MASTER) wtimer1 = omp_get_wtime()
#endif
                if (spline) then
                    if (cutmo) then
                        call aomocutspl_calc(ie, x, y, z, rrai(:, :, n))
                    else
                        call aomospl_calc(ie, x, y, z, rrai(:, :, n))
                    endif
                else
                    if (cutmo) then
                        call aomocut_calc(ie, x, y, z, rrai(:, :, n))
                    else
                        call aomo_calc(ie, x, y, z, rrai(:, :, n))
                    endif
                endif
#ifdef WTIMER
                if (MASTER) wtimer2 = omp_get_wtime()
#endif
                phi(n) = 0
                fgrad(:, n) = 0
                flapli(:, n) = 0
                flapl(n) = 0
                call mdetcalc(ie, 1, phi(n:n), fgrad(:, n:n), flapli(:, n:n), &
                        flapl(n:n), ierr)
#ifdef WTIMER
                if (MASTER) wtimer3 = omp_get_wtime()
#endif

                if (ierr /= MDET_NONE) then
                    if (mStopAtSingularity) then
                        write(iull, *) 'mdetcalc (AOMO) failed for:'
                        write(iull, '(i5,g20.10)') ierr, phi(1)
                        do i = 1, ne
                            write(iull, '(i5,3f14.7)') i, x(i), y(i), z(i)
                        end do
                        call abortp("eloc: singularity in determinant")
                    else
                        if (ierr == MDET_INV_ERR) then
                            error_code = ELOC_INV_ERR
                        else if (ierr == MDET_LU_ERR) then
                            error_code = ELOC_LU_ERR
                        end if
                        elSingularity = .true.
                        return
                    end if
                end if
                elPhi(n) = phi(n)
                elFgrad(:, n) = fgrad(:, n)
                elFlapl(n) = flapl(n)
                elFlapli(:, n) = flapli(:, n)

#ifdef WTIMER
                if (MASTER) then
                    wtimerPhi(3) = wtimerPhi(3) + (wtimer2 - wtimer1)
                    wtimerPhi(4) = wtimerPhi(4) + (wtimer3 - wtimer2)
                endif
#endif

            enddo ECLOOP

        else

#ifdef WTIMER
            if (MASTER) wtimer1 = omp_get_wtime()
#endif
            if (spline) then
                call aosplcalc(ie, ec, rrai)
            else
                call aocalc(ie, ec, rrai)
            endif
            !       // calculate MO's for current AO's
#ifdef WTIMER
            if (MASTER) wtimer2 = omp_get_wtime()
#endif
            call mocalc(ie)
#ifdef WTIMER
            if (MASTER) wtimer3 = omp_get_wtime()
#endif

            phi = 0
            fgrad = 0
            flapli = 0
            flapl = 0
            call mdetcalc(ie, nElecconfigs, phi, fgrad, flapli, flapl, ierr)
#ifdef WTIMER
            if (MASTER) wtimer4 = omp_get_wtime()
#endif
            if (ierr /= MDET_NONE) then
                if (mStopAtSingularity) then
                    write(iull, *) 'mdetcalc failed for:'
                    write(iull, '(i5,g20.10)') ierr, phi(1)
                    call eConfigArray_get(ec, 1, x, y, z)
                    do i = 1, ne
                        write(iull, '(3f14.7)') x(i), y(i), z(i)
                    end do
                    call abortp("eloc: singularity in determinant")
                else
                    if (ierr == MDET_INV_ERR) then
                        error_code = ELOC_INV_ERR
                    else if (ierr == MDET_LU_ERR) then
                        error_code = ELOC_LU_ERR
                    end if
                    elSingularity = .true.
                    return
                end if
            end if

            do n = 1, nElecConfigs
                elPhi(n) = phi(n)
                elFgrad(:, n) = fgrad(:, n)
                elFlapl(n) = flapl(n)
                elFlapli(:, n) = flapli(:, n)
            enddo

#ifdef WTIMER
            if (MASTER) then
                wtimerPhi(1) = (wtimer2 - wtimer1)
                wtimerPhi(2) = (wtimer3 - wtimer2)
                wtimerPhi(4) = (wtimer4 - wtimer3)
            endif
#endif

        endif
    end subroutine elocCalcPhi


    !     ----------------------------------------------------------------------------------------------------------
    subroutine elocCalcJas(ie, ec, optType, nElecConfigs, uu, uugrad, uulapl, uulapli, rrai, rrij, twoLevelStep, calcJas)
        !     ----------------------------------------------------------------------------------------------------------
        integer, intent(in) :: ie, nElecConfigs       ! 0: update all electron; i: update only electron i
        type(eConfigArray), intent(inout) :: ec  ! contains electron configurations
        character(len = *), intent(in) :: optType ! transferred to call to JasWDerivs

        integer :: n
        !     ! automatic arrays for actual electron number!
        real(r8) :: x(ne), y(ne), z(ne)
        real(r8), intent(in) :: rrai(ncenter, ne, eConfigArray_size(ec))
        real(r8), intent(in) :: rrij(ne, ne, eConfigArray_size(ec))
        real(r8) :: u, ugrad(3 * ne), ulapl, ulapli(ne)
        real(r8), intent(out) :: uu(eConfigArray_size(ec))
        real(r8), intent(out) :: uugrad(3 * ne, eConfigArray_size(ec))
        real(r8), intent(out) :: uulapl(eConfigArray_size(ec))
        real(r8), intent(out) :: uulapli(ne, eConfigArray_size(ec))
        integer, intent(in), optional :: twoLevelStep
        logical, intent(in), optional :: calcJas(eConfigArray_size(ec))
        integer :: step
        logical :: calc

        step = 3
        if(present(twoLevelStep)) step = twoLevelStep

        do n = 1, nElecConfigs
            if (step == 2 .and. jastype /= 'none') then
                !       ! calculate only jastrow value without derivatives
                u = 0
                ugrad = 0
                ulapl = 0
                ulapli = 0
                call jasCalc(.false., ie, rrai(:, :, n), rrij(:, :, n), u)
            else if (step == 3 .and. jastype /= 'none') then
                !       ! calculate jastrow
                call eConfigArray_get(ec, n, x, y, z)
                u = 0
                ugrad = 0
                ulapl = 0
                ulapli = 0
                calc = .true.
                if(present(calcJas)) calc = calcJas(n)

                if(calc) then
                    call jasCalcWDerivs(ie, x, y, z, rrai(:, :, n), rrij(:, :, n), optType, u, ugrad, ulapl, ulapli, n)
                endif
            else
                !       ! assign default values
                u = 0
                ugrad = 0
                ulapl = 0
                ulapli = 0
            endif
            uu(n) = u
            uugrad(:, n) = ugrad
            uulapl(n) = ulapl
            uulapli(:, n) = ulapli
            elU(n) = u
            elUGrad(:, n) = ugrad(:)
            elULapl(n) = ulapl
            elULapli(:, n) = ulapli
        enddo ! nElecConfigs

    end subroutine elocCalcJas

    !     -----------------------------------------------------------------------------------------------------
    subroutine elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal, &
            ecpNonlocal, ecpNLocalk, wfpDef, TMoveData)
        !     -----------------------------------------------------------------------------------------------------
        integer, intent(in) :: ie            ! ie==0: all electron, ie>0: update ie only
        integer, intent(in) :: nElecConfigs
        type(eConfigArray), intent(inout) :: ec            ! contains electron configurations
        real(r8), intent(in) :: rrai(:, :, :)   ! nuc-elec distances for all EConfigs
        real(r8), intent(in) :: rrij(:, :, :)   ! elec-elec distances for all EConfigs
        real(r8), intent(out) :: ecpLocal      ! local ECP contribution
        real(r8), intent(out) :: ecpNonlocal   ! localized nonlocal ECP contribution
        real(r8), intent(out), allocatable, optional :: ecpNLocalk(:) ! parameter derivative of nonlocal contrib
        type(WFParamDef), intent(in), optional :: wfpDef    ! wfParamDefinition with info about Jastrow terms
        type(TMoveDataType), intent(inout), optional :: TMoveData ! triggers T move calculation

        real(r8) :: x(size(rrij, 1)), y(size(rrij, 1)), z(size(rrij, 1))
        integer :: np, npJ, npJ1, npJ2, npJnl, npCI, npMO, i
        type(RdataUpdate) :: Rdu

        !        !
        !        ! check if  arg 'ie' makes any sense. One-electron moves? Really implemented?
        !        !

        call assert(nElecConfigs==1, "elocCalcECP: fast ECP Jastrow currently only for walker_block==1")
        !        ! careful: ec used for aos access for anisotropic Jastrow
        !        ! also: initialization of Jastrow terms for ECP updates relying on walker_block==1
        !        !
        ecpLocal = 0._r8
        ecpNonlocal = 0._r8
        call eConfigArray_get(ec, 1, x, y, z)

        !        !!!write(iul,'(a)') 'DBG:elocCalcECP:x,y,z:'
        !        !!!do i=1,ne
        !        !!!   write(iul,'(i3,3g20.10)') i,x(i),y(i),z(i)
        !        !!!enddo

        if (present(ecpNLocalk) .and. present(wfpDef)) then
            call wfparams_getNJastrowParams(wfpDef, npJ, npJ1, npJ2, npJnl)
            call wfparams_getNDetParams(wfpDef, npCI, npMO)
            np = npJ + npCI + npMO
            call Rdu%initWParamDerivs(x, y, z, npJ1, npJ2, npCI, npMO, rrai(:, :, 1), rrij(:, :, 1))  ! check if new (i.e. allocation) here OK
            Rdu%phi0 = elPhi(1)
            if (.not.allocated(ecpNLocalk)) then
                allocate(ecpNLocalk(np))
            else
                if (size(ecpNLocalk) /= np) call abortp("elocCalcECP: inconsistent sizes")
            end if
            ecpNLocalk = 0._r8
            call ecp%calculate(Rdu, ecpLocal, ecpNonlocal, ecpNonlocalk = ecpNLocalk)
        else if (present(TMoveData)) then
            !           ! possible T move
            call Rdu%init(x, y, z, rrai(:, :, 1), rrij(:, :, 1))                   !!! check if new (i.e. allocation) here OK
            Rdu%phi0 = elPhi(1)
            call ecp%calculate(Rdu, ecpLocal, ecpNonlocal, TMoveData = TMoveData)
            if (TMoveData%isMoved) then
                call eConfigArray_set(ec, 1, Rdu%x, Rdu%y, Rdu%z)
            end if
        else
            call Rdu%init(x, y, z, rrai(:, :, 1), rrij(:, :, 1))                   !!! check if new (i.e. allocation) here OK
            Rdu%phi0 = elPhi(1)
            call ecp%calculate(Rdu, ecpLocal, ecpNonlocal)
        end if

        call Rdu%delete()

    end subroutine elocCalcECP

    subroutine add_Vee_density_repulsion(nElecconfigs, ec)
        integer, intent(in) :: nElecConfigs
        type(eConfigArray), intent(inout) :: ec
        real(r8) :: x(ne), y(ne), z(ne)
        integer :: n, i

        do n = 1, nElecConfigs
            call eConfigArray_get(ec, n, x, y, z)
            do i = 1, ne
                elVee(n) = elVee(n) + get_density_repulsion([x(i), y(i), z(i)])
            end do
        end do
    end subroutine add_Vee_density_repulsion

END MODULE eloc_m