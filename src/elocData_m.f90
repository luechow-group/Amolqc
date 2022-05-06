! Copyright (C) 1999, 2015, 2018 Arne Luechow
! Copyright (C) 2013 Alexander Sturm
! Copyright (C) 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later


! elocData_m contains values pertaining to local energy
! that are accessed throughout the code
! All "elXXX" data contain just the "current" values
! whatever that means.
! The "elXXX" values are correct *after* a call to eloc(x,y,z)
! for the positions x,y,z.

! Currently: global data
! TODO: access functions get/set
!


!     -----------------
MODULE elocData_m
    !     -----------------

    use kinds_m, only : r8, i8
    use error_m
    use wfData_m
    implicit none

    ! main return data from eloc, for several elec configurations
    real(r8), allocatable :: elPhi(:)        ! Phi_G (determinantal part)
    real(r8), allocatable :: elU(:)          ! Jastrow exponent U
    real(r8), allocatable :: elEloc(:)       ! E_local(Psi_G)
    real(r8), allocatable :: elVen(:)        ! Ven potential energy
    real(r8), allocatable :: elVee(:)        ! Vee potential energy
    real(r8), allocatable :: elxDrift(:, :)   ! nabla(Psi_G) / Psi_G
    real(r8), allocatable :: elyDrift(:, :)   ! nabla(Psi_G) / Psi_G
    real(r8), allocatable :: elzDrift(:, :)   ! nabla(Psi_G) / Psi_G

    ! components of Eloc and Psi (only for one configuration)
    real(r8), allocatable :: elEkini(:, :)  ! -0.5*nabla_i**2 Psi_G / Psi_G
    real(r8), allocatable :: elUgrad(:, :)  ! grad(U)
    real(r8), allocatable :: elULapl(:)    ! laplacian(U)
    real(r8), allocatable :: elULapli(:,:)    ! laplacian(U(i)
    real(r8), allocatable :: elFgrad(:, :)  ! grad(Phi_G)
    real(r8), allocatable :: elFLapl(:)    ! laplacian(Phi_G)
    real(r8), allocatable :: elFLapli(:, :) ! laplacian(Phi_G(i,n))
    logical :: elSingularity ! singular determinant

    ! additional stuff
    real(r8) :: elECPPot, elECPPotNl ! Potential,Non-local part

    ! E_local cutoff (Rothstein/Vrbik)
    real(r8) :: mTauCutOff = 0
    real(r8) :: mCutOffFactor = 0
    logical :: mElocCutOff = .false.

    ! max and actual number of simultaneous electron configurations
    integer :: mMaxElecConfigs = 0
    integer :: mElecConfigs = 0
    ! Arrays for energy partitioning
    real(r8), allocatable :: elEkin_epart(:, :)
    real(r8), allocatable :: elVne_epart(:, :, :)
    real(r8), allocatable :: elVee_epart(:, :, :)

    ! stop when singularity in det is hit?
    logical :: mStopAtSingularity = .true.

    integer(i8) :: mElocCut, mElocCutCount
    integer(i8) :: mDriftCut, mDriftCutCount


CONTAINS


    !     -------------------------------
    subroutine eloc_initialize(nec)
        !     -------------------------------

        integer, intent(in), optional :: nec
        integer nn, alstat
        if (present(nec)) then
            mMaxElecConfigs = nec
        else
            mMaxElecConfigs = 1
        endif
        nn = mMaxElecConfigs

        call assert(ne>0 .and. ncenter>0, &
                'eloc_initialize: electron and nucleus number not set')
        if (allocated(elEloc)) then
            if (size(elEloc) /= nn) call eloc_deallocate()
        endif

        if (.not. allocated(elEloc)) then
            allocate(elPhi(nn), elU(nn), elEloc(nn), elVen(nn), elVee(nn),&
                    elxDrift(ne, nn), elyDrift(ne, nn), elzDrift(ne, nn), &
                    elEkini(ne, nn), elUgrad(3 * ne, nn), elULapl(nn), &
                    elFgrad(3 * ne, nn), elFLapl(nn), elFLapli(ne, nn), &
                    elULapli(ne, nn), stat = alstat)
            if(do_epart) allocate(elEkin_epart(ne, nn), elVne_epart(ncenter, ne, nn)&
                    , elVee_epart(ne, ne, nn))
            call assert(alstat==0, 'eloc_initialize: allocation3 failed')
            elPhi = 0; elU = 0; elEloc = 0; elVen = 0; elVee = 0; elxDrift = 0
            elyDrift = 0; elzDrift = 0; elEkini = 0; elUgrad = 0
            elFgrad = 0; elFLapli = 0; elULapli = 0; elECPPotNl = 0
        endif
    end subroutine eloc_initialize

    logical function  eloc_getStopAtSingularity()
        eloc_getStopAtSingularity = mStopAtSingularity
    end function  eloc_getStopAtSingularity

    subroutine eloc_setStopAtSingularity(l)
        logical, intent(in) :: l
        mStopAtSingularity = l
    end subroutine eloc_setStopAtSingularity

    logical function eloc_isAtSingularity()
        eloc_isAtSingularity = elSingularity
    end function eloc_isAtSingularity

    !     -----------------------------
    subroutine eloc_deallocate()
        !     -----------------------------

        call assert(allocated(elEloc), &
                'eloc: deallocation before allocation')
        deallocate(elPhi, elU, elEloc, elVen, elVee, elxDrift, elyDrift, elzDrift, &
                elEkini, elUgrad, elULapl, elFgrad, elFLapl, elFLapli)
        if(do_epart) deallocate(elEkin_epart, elVne_epart, elVee_epart)
    end subroutine eloc_deallocate

    !     -------------------------------------------------------
    subroutine eloc_getCurrentElocData(nec, phi, u, eloc, vpot, &
            xdrift, ydrift, zdrift)
        !     -------------------------------------------------------

        !     ! CHANGE to pointer!
        integer, intent(inout) :: nec
        real(r8), intent(inout) :: phi(:)
        real(r8), intent(inout) :: u(:)
        real(r8), intent(inout) :: eloc(:)
        real(r8), intent(inout) :: vpot(:)
        real(r8), intent(inout) :: xdrift(:, :), ydrift(:, :), zdrift(:, :)
        nec = mElecConfigs
        phi = elPhi
        u = elU
        eloc = elEloc
        vpot = elVen + elVee
        xdrift = elxDrift; ydrift = elyDrift; zdrift = elzDrift
    end subroutine eloc_getCurrentElocData

    !     -------------------------------------------------------
    subroutine eloc_getCurrentElocData1(phi, u, eloc, vpot, &
            xdrift, ydrift, zdrift)
        !     -------------------------------------------------------

        !     ! this returns non-array data containing the first elec config
        !     ! CHANGE to returning pointers to avoid copying
        real(r8), intent(inout) :: phi
        real(r8), intent(inout) :: u
        real(r8), intent(inout) :: eloc
        real(r8), intent(inout) :: vpot
        real(r8), intent(inout) :: xdrift(:), ydrift(:), zdrift(:)
        phi = elPhi(1)
        u = elU(1)
        eloc = elEloc(1)
        vpot = elVen(1) + elVee(1)
        xdrift = elxDrift(:, 1)
        ydrift = elyDrift(:, 1)
        zdrift = elzDrift(:, 1)
    end subroutine eloc_getCurrentElocData1

    subroutine eloc_getCurrentElocEpart(nec, kin, vee, vne)
        real(r8), intent(inout) :: kin(:)
        real(r8), intent(inout) :: vee(:, :)
        real(r8), intent(inout) :: vne(:, :)
        integer, intent(inout) :: nec

        kin = elEkin_epart(nec, :)
        vee = elVee_epart(nec, :, :)
        vne = elVne_epart(nec, :, :)

    end subroutine eloc_getCurrentElocEpart

    subroutine eloc_getCurrentElocEpart1(kin, vee, vne)
        real(r8), intent(inout) :: kin(:)
        real(r8), intent(inout) :: vee(:, :)
        real(r8), intent(inout) :: vne(:, :)

        kin = elEkin_epart(1, :)
        vee = elVee_epart(1, :, :)
        vne = elVne_epart(1, :, :)

    end subroutine eloc_getCurrentElocEpart1

    subroutine eloc_getECPContribs1(EECPlocal, EECPnonlocal)
        real(r8), intent(inout) :: EECPlocal, EECPnonlocal
        EECPlocal = elECPPot
        EECPnonlocal = elECPPotNl
    end subroutine


    !     ---------------------------------------
    subroutine setElocCutOff(cutoff, tau, cf)
        !     ---------------------------------------

        logical, intent(in) :: cutoff
        real(r8), intent(in) :: tau
        real(r8), intent(in) :: cf

        mElocCutOff = cutoff
        mTauCutOff = tau
        mCutOffFactor = cf
    end subroutine setElocCutOff


    !     ---------------------------------------
    subroutine getElocCutOff(cutoff, tau, cf)
        !     ---------------------------------------

        logical, intent(out) :: cutoff
        real(r8), intent(out), optional :: tau
        real(r8), intent(out), optional :: cf

        cutoff = mElocCutOff
        if (present(tau)) tau = mTauCutOff
        if (present(cf)) cf = mCutOffFactor
    end subroutine getElocCutOff


END MODULE elocData_m

