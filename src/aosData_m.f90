! Copyright (C) 1998, 2012, 2014, 2018 Arne Luechow
! Copyright (C) 1999 Sebastian Manten
! Copyright (C) 1999 Tony C. Scott
! Copyright (C) 2013 Alexander Sturm
! Copyright (C) 2015 Christoph Schulte
!
! SPDX-License-Identifier: GPL-3.0-or-later

!     --------------
MODULE aosData_m
    !     --------------
    !
    ! contains the data of AO module including only constructors and
    ! desctructors. AO calculation is in module "aos"

      use kinds_m, only: r8
      use globalUtils_m
      use error_m
      use wfData_m
      use cspline_m
      use cuspOpt_m
      use utils_m
      use atom_m, only: getPSEIdx, getPSESymbol
      use wfData_m
#ifdef MPI
      use MPI_F08
#endif
      implicit none

    public
    private :: readGaussianBasis, readAmolqcBasis, &
            calcNormFactors, getso

    ! definition of the basis functions

    integer :: ngto(basmax)            ! # GTOs in contraction
    integer :: so(basmax)              ! spline function # of basis functions
    integer :: somax = 0                 ! total # of functions to spline in WF

    real(r8) bzet(basmax)                ! orbital exponents (zeta)
    real(r8) norm(basmax)                ! norm factors
    real(r8) cntrctn(2, cnmax, basmax)     ! GTO contraction (alp_i,c_i!)
    real(r8) stoc(4, basmax)              ! coeff. for correction of cusp
    !                                      ! 1. exp. Koeffizient ,< 0.2 fuer Opt
    !                                      ! 2. Vorfaktor
    !                                      ! 3. i Anfang Uebergangsbereich in Bohr
    !                                      ! 4. i Ende Uebergangsbereich in Bohr

    character bl(basmax)*1             ! l of orbs as 'S','P','D','F','G'
    character typ(basmax)*3           ! 'STO' or 'GTO'

    !     each basis function (STO, GTO, or contracted GTO) evaluated
    !     for each electron position, of each electron configuration
    real(r8), allocatable :: uao(:, :, :)
    real(r8), allocatable :: uxao(:, :, :), uyao(:, :, :), uzao(:, :, :)
    real(r8), allocatable :: u2ao(:, :, :)

    !cc stuff for AO - Cutoff
    real(r8), allocatable :: aocuts(:)
    !     real(r8), aocutsp(basmax)  ! AO Cutoff array (distances) for
    !                                        ! laplacian       : aocuts
    !                                        ! and the orbitals: aucutsp
    ! variables for LINSCAL only
    !      integer AOSofCenters(2,basmax)     ! first + last AOS of centers
    !      integer AOSListsI(basmax,nmax)     ! lists of non-zero diff. AO's of the electrons
    !      integer AOSListsIEnd(nmax)       ! end  of the lists
    !      integer AOSListsf(basmax)          !
    !      integer AOSListsfEnd             !
    !      integer CenterofAO(basmax)         ! index of the center of the AO
    !      integer basftobas(basmax)          ! first bas of basf , basf index diff. b.function
    !                                      !                     bas  index indi. b.function
    !      real(r8) AORadius(basmax)            ! cutoff-radii of the AOS, distinct AO types
    !      real(r8) AOSRadius(3,basmax)         ! cutoff-radii of the AOS, individual bf

    integer :: mAOElecConfigs = 0   ! arrays allocated for n electron configurations

    type AOIndex
        integer :: nums(amax) = 0
        integer :: nump(amax) = 0
        integer :: numd(amax) = 0
        integer :: numf(amax) = 0         ! number of s,p,d,f basis functions at each nuc
        integer :: numg(amax) = 0
        integer :: s(amax, basmax) = 0    ! index of basis function (basf)
        integer :: p(amax, basmax) = 0    ! for each nucleus.
        integer :: d(amax, basmax) = 0
        integer :: f(amax, basmax) = 0
        integer :: g(amax, basmax) = 0
    end type AOIndex

    type(AOIndex) :: AOIdx

    !     !!!save cntrctn


CONTAINS

    !     ------------------------------
    subroutine aos_initialize(nec)
        !     ------------------------------
        integer, intent(in), optional :: nec
        integer nn, alstat
        if (present(nec)) then
            mAOElecConfigs = nec
        else
            mAOElecConfigs = 1
        endif
        nn = mAOElecConfigs

        call assert(nbas>0 .and. ne>0, &
                'aos_initialize: nbas or ne not yet set')

        if (allocated(uao)) then
            call aos_deallocate()
        endif

        allocate(uao(nbas, ne, nn), uxao(nbas, ne, nn), uyao(nbas, ne, nn), &
                uzao(nbas, ne, nn), u2ao(nbas, ne, nn), stat = alstat)
        call assert(alstat==0, 'aos_initialize: allocation failed')

    end subroutine aos_initialize

    !     ---------------------------
    subroutine aos_deallocate()
        !     ---------------------------

        integer alstat
        if (allocated(uao)) then
            deallocate(uao, uxao, uyao, uzao, u2ao, stat = alstat)
            call assert(alstat==0, 'aos_deallocate: failed')
        endif
    end subroutine aos_deallocate


    !     -----------------------------
    subroutine aoinputg(lines)
        !     -----------------------------

        !     basis set input in Gaussian format to be read string array lines(nl)
        !     basis input assumed to start in line 2

        character(len = *), intent(in) :: lines(:)! lines array
        integer k, socounter, iflag, io, idx, iu

        nbas = 0                        ! counts AO's, individual basis functions
        nbasf = 0                       ! counts distinct AO types
        socounter = 0                   ! counts the contracted basis functions
        idx = 2
        do k = 1, ncenter
            read(lines(idx), '(A)', iostat = io) atoms(k)%ba           ! name of atom(!) instead of basis name
            call assert(io==0, '(aoinputg) expecting basis set name')
            idx = idx + 1
            call readGaussianBasis(lines, idx, k, socounter, iflag)
        enddo
        somax = socounter

        call calcNormFactors()

        !     ! initialize spline tables if spline option set
        if (spline) then
            call aospline(splinpnts)
        endif
        if (aosopt) then
           iu = iul
           call aooutputg(iu)
           call stopp('Normal termination after AO spline optimization')
        end if   
        if (cutao) then
            call ao_cut(aocutoff)
        endif

        call aos_initialize()

    end subroutine aoinputg


    !     ----------------------
    subroutine aoinputex()
        !     ----------------------

        ! read basis set from basisset library

        integer, parameter :: MAXLINES = 5000
        integer, parameter :: MAXLEN = 120
        character(len = MAXLEN) :: lines(MAXLINES) ! lines array
        integer :: nl      ! actual # of lines
        integer k, socounter, cha, io, iflag, idx, iu
        character(len=20) basislist
        character(len=180) basispath, basisname
        logical fileExists

        nbas = 0                        ! counts AO's, individual basis functions
        nbasf = 0                       ! counts distinct AO types
        socounter = 0                   ! counts the contracted basis functions
        call getAmolqcPath(basispath)
        call assert(len(trim(basispath))<176, &
                "aoinputex: amolqc path length exceeds definition")
        basispath = trim(basispath) // "/bib"

        if (basis == 'diff') then
            do k = 1, ncenter
                basisname = trim(basispath) // '/' // trim(atoms(k)%ba)
                if (MASTER) then
                    inquire(file = basisname, exist = fileExists)
                    if (.not.fileExists)&
                            call abortp('(aoinputex,diff): basis file not found')
                end if
                call readFileParallel(mytid, basisname, lines, nl)
                idx = 1
                do
                    if (lines(idx)(1:4) == "####") then
                        read(lines(idx + 1), '(I2)', iostat = io) cha
                        call assert(io==0, '(aoinputex,diff): expecting charge')
                        idx = idx + 2
                        if (cha == atoms(k)%elemIdx) then
                            call readGaussianBasis(lines, idx, k, &
                                    socounter, iflag)
                            exit
                        endif
                    else
                        if (nl == idx + 1) then
                            call abortp('(readbasis,diff): atom basis not avaliable')
                        end if
                        idx = idx + 1
                    end if
                end do
            enddo
            somax = socounter
        else
            write(basislist, '(A20)') basis
            basisname = trim(basispath) // '/' // trim(basislist)
            if (MASTER) then
                inquire(file = basisname, exist = fileExists)
                if (.not.fileExists)&
                        call abortp('(aoinputex): basis file not found')
            end if
            call readFileParallel(mytid, basisname, lines, nl)
            do k = 1, ncenter
                idx = 1
                do
                    if (lines(idx)(1:4) == "####") then
                        read(lines(idx + 1), '(I2)', iostat = io) cha
                        call assert(io==0, '(aoinputex) expecting charge')
                        idx = idx + 2
                        if (cha == atoms(k)%elemIdx) then
                            call readGaussianBasis(lines, idx, k, &
                                    socounter, iflag)
                            exit
                        else if (cha==0) then
                            call abortp('(readbasis): atom basis not avaliable')
                        endif
                    else
                        idx = idx + 1
                    endif
                enddo
            enddo
            somax = socounter
        end if

        call calcNormFactors()

        !     ! initialize spline tables if spline option set
        if (spline) then
            call aospline(splinpnts)
        endif
        if (aosopt) then
           iu = iul
           call aooutputg(iu)
           call abortp('Normal termination after AO spline optimization')
        end if   

        if (cutao) then
            call ao_cut(aocutoff)
        endif

        call aos_initialize()

    end subroutine aoinputex


    !     -----------------------
    subroutine aoinputabs()
        !     -----------------------

        ! reading basis set in amolqc abs form
        integer, parameter :: MAXLINES = 10000
        integer, parameter :: MAXLEN = 120
        character(len = MAXLEN) :: lines(MAXLINES) ! lines array
        integer :: nl      ! actual # of lines
        integer k, socounter, cha, iflag, idx, iu
        character(len=180) basispath
        character(len = 150) basisFileName
        logical fileExists

        nbas = 0                        ! counts AO's, individual basis functions
        nbasf = 0                       ! counts distinct AO types
        socounter = 0                   ! counts the contracted basis functions
        call getAmolqcPath(basispath)
        call assert(len(trim(basispath))<176, &
                "aoinputabs: amolqc path length exceeds definition")
        basispath = trim(basispath) // "/bib"

        if (basis=='diff') then
            do k = 1, ncenter
                basisFileName = trim(basispath) // '/' // trim(atoms(k)%ba) // '.abs'
                if (MASTER) then
                    inquire(file = basisFileName, exist = fileExists)
                    if (.not.fileExists) then
                        call abortp('(aoinputabs,diff): basis file not found')
                    endif
                end if
                call readFileParallel(mytid, basisFileName, lines, nl)
                idx = 1
                do
                    if (lines(idx)(1:4) == "****") then
                        cha = getPSEIdx(lines(idx + 1)(1:2))
                        idx = idx + 2
                        if (cha == atoms(k)%elemIdx) then
                            call readAmolqcBasis(lines, idx, k, socounter, iflag)
                            exit
                        else if (cha==0) then
                            call abortp('(aoinputabs,diff): atom basis not avaliable')
                        end if
                    else
                        idx = idx + 1
                    end if
                end do
            enddo
            somax = socounter
        else
            basisFileName = trim(basispath) // '/' // trim(basis) // '.abs'
            if (MASTER) then
                inquire(file = basisFileName, exist = fileExists)
                if (.not.fileExists)&
                        call abortp('(aoinputabs): basis file not found')
            end if
            call readFileParallel(mytid, basisFileName, lines, nl)

            do k = 1, ncenter
                idx = 1
                do
                    if (lines(idx)(1:4) == "****") then
                        cha = getPSEIdx(lines(idx + 1)(1:2))
                        idx = idx + 2
                        if (cha == atoms(k)%elemIdx) then
                            call readAmolqcBasis(lines, idx, k, socounter, iflag)
                            exit
                        else if (cha==0) then
                            call abortp('(aoinputabs): atom basis not avaliable')
                        end if
                    else
                        idx = idx + 1
                    end if
                end do
            end do
            somax = socounter
        endif

        if (somax==0) spline = .false.

        call calcNormFactors()

        !     ! initialize spline tables if spline option set
        if (spline) then
            call aospline(splinpnts)
        endif
        if (aosopt) then
           iu = iul
           call aooutput(iu)
           call abortp('Normal termination after AO spline optimization')
        end if   

        if (cutao) then
            call ao_cut(aocutoff)
        endif

        call aos_initialize()

    end subroutine aoinputabs


    !     -----------------------
    subroutine aoinputgbs()
        !     -----------------------

        ! reading basis set in gaussian gbs form
        integer, parameter :: MAXLINES = 5000
        integer, parameter :: MAXLEN = 120
        character(len = MAXLEN) :: lines(MAXLINES) ! lines array
        integer :: nl      ! actual # of lines
        integer k, socounter, cha, io, iflag, idx, iu
        character(len=180) basispath, basisname
        logical fileExists

        nbas = 0                        ! counts AO's, individual basis functions
        nbasf = 0                       ! counts distinct AO types
        socounter = 0                   ! counts the contracted basis functions
        call getAmolqcPath(basispath)
        call assert(len(trim(basispath))<176, &
                "aoinputabs: amolqc path length exceeds definition")
        basispath = trim(basispath) // "/bib"
        basisname = trim(basispath) // '/' // basis
        if (MASTER) then
            inquire(file = basisname, exist = fileExists)
            if (.not.fileExists)&
                    call abortp('(aoinputgbs): basis file not found')
        end if
        call readFileParallel(mytid, basisname, lines, nl)

        do k = 1, ncenter
            idx = 1
            do
                if (lines(idx)(1:1) == "-") then
                    read(lines(idx + 1), '(I2)', iostat = io) cha
                    call assert(io==0, '(aoinputgbs): expecting charge')
                    idx = idx + 2
                    if (cha == atoms(k)%elemIdx) then
                        call readGaussianBasis(lines, idx, k, &
                                socounter, iflag)
                        exit
                    else if (cha==0) then
                        call abortp('(aoinputgbs): atom basis not avaliable')
                    endif
                else
                    idx = idx + 1
                endif
            enddo
        enddo
        somax = socounter
        if (somax == 0) spline = .false.

        call calcNormFactors()

        !     ! initialize spline tables if spline option set
        if (spline) then
            call aospline(splinpnts)
        endif
        if (aosopt) then
           iu = iul
           call aooutputg(iu)
           call abortp('Normal termination after AO spline optimization')
        end if   

        if (cutao) then
            call ao_cut(aocutoff)
        endif

        call aos_initialize()

    end subroutine aoinputgbs


    !     ----------------------------
    subroutine aoinput(lines)
        !     ----------------------------

        ! read general basis set from lines array
        ! then the normalization factors are calculated.
        character(len = *), intent(in) :: lines(:)! lines array
        integer idx, k, socounter, io, iflag, iu

        nbas = 0                         ! counts AO's, individual basis functions
        nbasf = 0                        ! counts distinct AO types
        socounter = 0

        idx = 2
        do k = 1, ncenter
            read(lines(idx), '(A)', iostat = io) atoms(k)%ba           ! name for basis of atom
            call assert(io==0, '(aoinput) expecting basis set name')
            idx = idx + 1
            call readAmolqcBasis(lines, idx, k, socounter, iflag)
            idx = idx + 1
        enddo
        somax = socounter

        call calcNormFactors()

        if (somax==0) spline = .false.   ! no splines necessary

        !     ! initialize spline tables if spline option set
        if (spline) then
            call aospline(splinpnts)
        endif
        if (aosopt) then
           iu = iul
           call aooutput(iu)
           call abortp('Normal termination after AO spline optimization')
        end if   

        if (cutao) then
            call ao_cut(aocutoff)
        endif

        call aos_initialize()

    end subroutine aoinput


    !     ------------------------
    subroutine aooutputg(iu)
        !     ------------------------
        !
        ! AO basis in gaussian forma like on input
        integer iu    ! file unit for writing
        integer i, j, k
        real(r8) scal

        scal = 1.0d0    ! scaling factor currently not used!
        i = 1
        do k = 1, ncenter
            write(iu, '(A)') getPSESymbol(atoms(k)%elemIdx)
            do while (bc(i)==k .and. i<=nbasf)
                if (cuspcor) then
                    write(iu, '(a2,i5,g10.3,4(G13.6))')  bl(i), ngto(i), scal, &
                            stoc(1, i), stoc(2, i), stoc(3, i), stoc(4, i)
                else
                    write(iu, '(a2,i5,g10.3)')  bl(i), ngto(i), scal
                endif
                call assert(ngto(i)<=cnmax, &
                        '(aooutputg): GTO contraction length too large')
                do j = 1, ngto(i)
                    write(iu, *) cntrctn(1, j, i), cntrctn(2, j, i)
                enddo
                i = i + 1
            enddo
            write(iu, '(a)') '****'
        enddo

    end subroutine aooutputg


    !     -----------------------
    subroutine aooutput(iu)
        !     -----------------------

        ! aoinput writes basis functions to log file unit 'iu'
        integer iu    ! file unit for writing
        integer i, j, k

        i = 1
        do k = 1, ncenter
            write(iu, '(A)') getPSESymbol(atoms(k)%elemIdx)
            do while (bc(i)==k .and. i<=nbasf)
                if (typ(i) == 'STO') then
                    write(iu, '(I1,A1,1X,A3,1X,G14.5)')&
                            bn(i), bl(i), typ(i), bzet(i)
                else if (typ(i) == 'GTO') then
                    if (cuspcor) then
                        write(iu, '(I1,A1,1X,A3,1X,I4,4(1X,G10.4))')&
                                bn(i), bl(i), typ(i), ngto(i), &
                                stoc(1, i), stoc(2, i), stoc(3, i), stoc(4, i)
                    else
                        write(iu, '(I1,A1,1X,A3,1X,I4)')&
                                bn(i), bl(i), typ(i), ngto(i)
                    endif
                    do j = 1, ngto(i)
                        write(iu, '(I5,2G14.7)') j, cntrctn(1, j, i), &
                                cntrctn(2, j, i)
                    enddo
                endif
                i = i + 1
            enddo
            write(iu, '(a)') '****'
        enddo

    end subroutine aooutput


    !     -----------------------
    subroutine aospline(np)
        !     -----------------------

        ! Version 1.0 (22.11.97): spline radial parts of AO's
        !                         modified Version 2.0 of aocalc (1.1 of getaos)
        !                         splines are used only for contracted GTO's
        !
        !     09.09.1999 SM     : cusp correction option added
        !                         analytic derivatives at r = 0

        ! parameter:
        integer np           ! number of spline points

        ! constants:
        real(r8), parameter :: sqr3 = 1.73205080756887729d0
        real(r8), parameter :: sqr5 = 2.236067977499789696d0
        real(r8), parameter :: alpha = 0.1d0    ! spline mapping factor

        ! variables
        integer i, bf, bf1, ic, k0, ispl
        real(r8) r1, r2, alp, u, ux, a, c
        real(r8) y0(np), y1(np), y2(np)
        real(r8) y00(6)           ! 1,2: y'(0)

        ! bf refers to the degenerate set of cartesian
        ! basis function (S:1,P:3,D:6,F:10) as input, which may be of type STO
        ! or contracted GTO.
        ! al refers to the individual basis function, as used in LCAO-MO's.
        ! (composed in subroutine mdetwf)
        ! i refers to the current electron.

        ! note: the so(bf) array maps basis functions to spline entries (objects)
        ! same basis function for same atoms can share spline object (same radial function)

        ! TODO: construct so(bf) table here using basis set string (like "c_avtz")

        if (somax == 0) then
            call abortp('(aospline): somax must not be 0')
        endif
        call csplinit(3 * somax, np, alpha)   ! 3 splines per functions (f,f',f'')
        call csplxarray(np)

        if (logmode >= 1) then
            write(iul, *) ' splining GTOs with ', np, ' points'
            if (cuspcor.and.logmode==2) write(iul, *)&
                    ' correcting cusp of the following basis functions:'
        endif

        !-----Calculation of radial parts of the AO's and their derivatives at spline points
        do ispl = 1, somax                 ! loop over spline-functions
            bf = 0
            do bf1 = 1, nbasf               ! search-loop over basis functions
                if (ispl==so(bf1)) then
                    bf = bf1
                    exit
                endif
            enddo

            if (bf==0) cycle             ! not found

            if (typ(bf) /= 'GTO') cycle   ! only GTO's are splined
            !        // only primitive cartesian gaussians: 1s,2p,3d,4f
            !        // i.e. no r factor. Thus nn is not used here.

            y00(1) = 0d0
            y00(2) = 0d0
            y00(3) = 0d0
            y00(4) = 0d0
            y00(5) = 0d0
            y00(6) = 0d0

            if (bl(bf) == 'S') then               ! 1s GTO
                do ic = 1, ngto(bf)                     ! loop over contraction
                    y00(1) = y00(1) + 0d0
                    y00(2) = y00(2) - 2d0 * cntrctn(1, ic, bf) * cntrctn(2, ic, bf)
                    y00(3) = y00(3) - 2d0 * cntrctn(1, ic, bf) * cntrctn(2, ic, bf)
                    y00(4) = y00(4) + 0d0
                    y00(5) = y00(5) + 0d0
                    y00(6) = y00(6) + 12d0 * cntrctn(1, ic, bf)**2 * cntrctn(2, ic, bf)
                enddo
                do i = 1, np - 1
                    r1 = csplx(i)
                    r2 = csplx(i)**2
                    y0(i) = 0.d0
                    y1(i) = 0.d0
                    y2(i) = 0.d0
                    do ic = 1, ngto(bf)                     ! loop over contraction
                        alp = cntrctn(1, ic, bf)
                        u = cntrctn(2, ic, bf) * exp(-alp * r2)
                        y0(i) = y0(i) + u
                        y1(i) = y1(i) - u * 2d0 * alp * r1
                        y2(i) = y2(i) + u * (4d0 * alp**2 * r2 - 2d0 * alp)
                    enddo
                enddo
                if (cuspcor) then                       ! Cusp-Korrektur
                    if (stoc(1, bf)/=(0d0)) then
                        call cuspcorrect (ngto, cntrctn, bf, REAL(atoms(bc(bf))%za, r8), &
                                stoc(1, bf), &
                                stoc(2, bf), stoc(3, bf), stoc(4, bf), &
                                y0, y1, y2, np, k0)
                        a = stoc(1, bf)
                        c = stoc(2, bf)
                        y00(1) = -a * c
                        y00(2) = a * a * c
                        y00(3) = a * a * c
                        y00(4) = -a * a * a * c
                        y00(5) = -a * a * a * c
                        y00(6) = a**4 * c
                    endif
                endif

            else if (bl(bf) == 'P') then          ! 2p GTO's
                do ic = 1, ngto(bf)                     ! loop over contraction
                    y00(1) = y00(1) + 0d0
                    y00(2) = y00(2) - 2d0 * cntrctn(1, ic, bf) * cntrctn(2, ic, bf)
                    y00(3) = y00(3) + 0d0
                    y00(4) = y00(4) + 4d0 * cntrctn(1, ic, bf)**2 * cntrctn(2, ic, bf)
                    y00(5) = y00(5) + 0d0
                    y00(6) = y00(6) + 28d0 * cntrctn(1, ic, bf)**2 * cntrctn(2, ic, bf)
                enddo
                do i = 1, np - 1
                    r2 = csplx(i)**2
                    y0(i) = 0.d0
                    y1(i) = 0.d0
                    y2(i) = 0.d0
                    do ic = 1, ngto(bf)                      ! loop over contraction
                        alp = cntrctn(1, ic, bf)
                        u = cntrctn(2, ic, bf) * exp(-alp * r2)
                        ux = -2d0 * alp * u
                        y0(i) = y0(i) + u
                        y1(i) = y1(i) + ux
                        y2(i) = y2(i) + ux * (5d0 - 2d0 * alp * r2)
                    enddo
                enddo
            else if (bl(bf) == 'D') then         ! 3d GTO
                do ic = 1, ngto(bf)                     ! loop over contraction
                    y00(1) = y00(1) + 0d0
                    y00(2) = y00(2) - 2d0 * cntrctn(1, ic, bf) * cntrctn(2, ic, bf)
                    y00(3) = y00(3) + 0d0
                    y00(4) = y00(4) + 4d0 * cntrctn(1, ic, bf)**2 * cntrctn(2, ic, bf)
                    y00(5) = y00(5) + 0d0
                    y00(6) = y00(6) + 36d0 * cntrctn(1, ic, bf)**2 * cntrctn(2, ic, bf)
                enddo
                do i = 1, np - 1
                    r2 = csplx(i)**2
                    y0(i) = 0.d0
                    y1(i) = 0.d0
                    y2(i) = 0.d0
                    do ic = 1, ngto(bf)                      ! loop over contraction
                        alp = cntrctn(1, ic, bf)
                        u = cntrctn(2, ic, bf) * exp(-alp * r2)
                        ux = -2d0 * alp * u
                        y0(i) = y0(i) + u
                        y1(i) = y1(i) + ux
                        y2(i) = y2(i) + ux * (7d0 - 2d0 * alp * r2)
                    enddo
                enddo
            else if (bl(bf) == 'F') then         ! 4f GTO
                do ic = 1, ngto(bf)                     ! loop over contraction
                    y00(1) = y00(1) + 0d0
                    y00(2) = y00(2) - 2d0 * cntrctn(1, ic, bf) * cntrctn(2, ic, bf)
                    y00(3) = y00(3) + 0d0
                    y00(4) = y00(4) + 4d0 * cntrctn(1, ic, bf)**2 * cntrctn(2, ic, bf)
                    y00(5) = y00(5) + 0d0
                    y00(6) = y00(6) + 44d0 * cntrctn(1, ic, bf)**2 * cntrctn(2, ic, bf)
                enddo
                do i = 1, np - 1
                    r2 = csplx(i)**2
                    y0(i) = 0.d0
                    y1(i) = 0.d0
                    y2(i) = 0.d0
                    do ic = 1, ngto(bf)                      ! loop over contraction
                        alp = cntrctn(1, ic, bf)
                        u = cntrctn(2, ic, bf) * exp(-alp * r2)
                        ux = -2d0 * alp * u
                        y0(i) = y0(i) + u
                        y1(i) = y1(i) + ux
                        y2(i) = y2(i) + ux * (9d0 - 2d0 * alp * r2)
                    enddo
                enddo
            else
                call abortp('(aospline): spline interpolation implemented only for s,p,d,f')
            endif  ! bl

            !        ! call cubic spline routine. evaluate cubic splines for all
            !        ! three functions. Faster would be spline evaluation only of y0
            !        ! and using derivatives of spline polynomial.
            y0(np) = 0
            y1(np) = 0
            y2(np) = 0
            call cspline(3 * ispl - 2, y0, np, y00(1), y00(2))
            call cspline(3 * ispl - 1, y1, np, y00(3), y00(4))
            call cspline(3 * ispl, y2, np, y00(5), y00(6))

        enddo  ! ispl-loop over spline-functions

        if (logmode >=2) write(iul, *)

    end subroutine aospline


    !     ------------------------------------------------
    logical function inquireBasisFile(basisFileName)
        !     ------------------------------------------------

        character(len = *), intent(in) :: basisFileName
        character(len = 180) :: basispath
        logical fexist
#ifdef MPI
        integer ierr
#endif

      call getAmolqcPath(basispath)
      basispath = trim(basispath)//'/bib'
      if (MASTER) then
         inquire(file=trim(basispath)//'/'//basisFileName,exist=fexist)
      end if
#ifdef MPI
      call mpi_bcast(fexist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif
      inquireBasisFile = fexist

    end function inquireBasisFile


    !     ----------------------------------------------------------
    subroutine readAmolqcBasis(lines, idx, k, socounter, iflag)
        !     ----------------------------------------------------------

        character(len = *), intent(in) :: lines(:)! lines array
        integer, intent(inout) :: idx     ! idx in line array
        !                                             ! on return: next line
        integer, intent(in) :: k           ! index of nucleus
        integer, intent(inout) :: socounter   ! counter for functions to spline
        integer, intent(out) :: iflag       ! indicate if properly finished
        integer i, j, io, nToken
        character(len = 40) :: token(10)

        do i = nbasf + 1, basmax              ! 'infinite' loop (basmax=max no. of orbs)
            call tokenize(lines(idx), token, nToken)
            if (token(1)(1:1) == '*') exit
            if (.not.(token(2)=='GTO' .or. token(2)=='STO')) then
                call abortp("readAmolqcBasis: illegal format reading basis")
            end if
            typ(i) = token(2)(1:3)
            if (typ(i) == 'STO') then
                read(token(1), '(i1)') bn(i)   ! convert to integer!
                bl(i) = token(1)(2:2)
                if (.not. (bn(i)>=1 .and. bn(i)<=5)) then
                    call abortp("readAmolqcBasis: n must be 1 ... 5")
                end if
                if (.not. (bl(i)=='S' .or. bl(i)=='P' .or. bl(i)=='D' .or. bl(i)=='F' .or. bl(i)=='G')) then
                    call abortp("readAmolqcBasis: l must be S,P,D,F,G")
                end if
                read(token(3), *, iostat = io) bzet(i)       ! convert to double
                if (io /= 0) call abortp("readAmolqcBasis: zeta not read")
                idx = idx + 1
            else if (typ(i) == 'GTO') then
                if (len(trim(token(1))) /= 1) then
                    call abortp("readAmolqcBasis: GTO requires S,P,D,F,G only")
                end if
                bl(i) = token(1)(1:1)
                if (bl(i) == 'S') then
                    bn(i) = 1
                else if (bl(i) == 'P') then
                    bn(i) = 2
                else if (bl(i) == 'D') then
                    bn(i) = 3
                else if (bl(i) == 'F') then
                    bn(i) = 4
                else if (bl(i) == 'G') then
                    bn(i) = 5
                else
                    call abortp("readAmolqcBasis: GTO only S,P,D,F,G")
                end if
                read(token(3), *) ngto(i)
                if (cuspcor .and. ngto(i) > 1) then     ! cusp correction parameters
                    if (nToken < 7) then
                        call abortp("readAmolqcBasis: expecting cusp" // &
                                " correction parameters")
                    endif
                    read(token(4), *) stoc(1, i)
                    read(token(5), *) stoc(2, i)
                    read(token(6), *) stoc(3, i)
                    read(token(7), *) stoc(4, i)
                    call getso(socounter, i, k)
                else if (spline) then
                    call getso(socounter, i, k)
                end if
                call assert(ngto(i)<=cnmax, &
                        "readAmolqcBasis: GTO contraction length too large")
                do j = 1, ngto(i)
                    read(lines(idx + j), *) cntrctn(1, j, i), cntrctn(2, j, i)
                end do
                idx = idx + ngto(i) + 1
            end if
            nbasf = nbasf + 1

            if (bl(i) == 'S') then
                nbas = nbas + 1
                AOIdx%nums(k) = AOIdx%nums(k) + 1
                AOIdx%s(k, AOIdx%nums(k)) = nbas
            else if (bl(i) == 'P') then
                nbas = nbas + 3
                AOIdx%nump(k) = AOIdx%nump(k) + 1
                AOIdx%p(k, AOIdx%nump(k)) = nbas - 2   ! first idx of px,py,pz
            else if (bl(i) == 'D') then
                nbas = nbas + 6
                AOIdx%numd(k) = AOIdx%numd(k) + 1
                AOIdx%d(k, AOIdx%numd(k)) = nbas - 5
            else if (bl(i) == 'F') then
                nbas = nbas + 10
                AOIdx%numf(k) = AOIdx%numf(k) + 1
                AOIdx%f(k, AOIdx%numf(k)) = nbas - 9
            else if (bl(i) == 'G') then
                nbas = nbas + 15
                AOIdx%numg(k) = AOIdx%numf(k) + 1
                AOIdx%f(k, AOIdx%numg(k)) = nbas - 14
            end if
            bc(i) = k
        end do
        if (i>basmax) iflag = 1

    end subroutine readAmolqcBasis


    !     ------------------------------------------------------------
    subroutine readGaussianBasis(lines, idx, k, socounter, iflag)
        !     ------------------------------------------------------------

        character(len = *), intent(in) :: lines(:)! lines array
        integer, intent(inout) :: idx     ! idx in line array
        !                                             ! on return: next line
        integer, intent(in) :: k           ! index of nucleus
        integer, intent(inout) :: socounter   ! counter for functions to spline
        integer, intent(out) :: iflag       ! indicate if properly finished
        integer i, j, nToken
        !     !!!character(len=120)     :: line
        character(len = 40) :: token(10)

        iflag = 0
        do i = nbasf + 1, basmax              ! 'infinite' loop (basmax=max no. of orbs)
            !        !!!read(iu,'(A)') line
            call tokenize(lines(idx), token, nToken)
            if (token(1)(1:1) == '*') exit
            typ(i) = 'GTO'
            if (len(trim(token(1))) /= 1) then
                call abortp("readGaussianBasis: GTO requires S,P,D,F,G only")
            end if
            bl(i) = strToUpper(token(1)(1:1))
            if (.not.(bl(i)=='S'.or.bl(i)=='P'.or.bl(i)=='D'&
                    .or.bl(i)=='F' .or.bl(i)=='G')) then
                call abortp("readGaussianBasis: GTO only S,P,D,F,G")
            end if
            read(token(2), *) ngto(i)
            if (cuspcor .and. ngto(i) > 1) then     ! cusp correction parameters
                if (nToken < 7) then
                    call abortp("readGaussianBasis: expecting cusp" // &
                            " correction parameters")
                end if
                read(token(4), *) stoc(1, i)
                read(token(5), *) stoc(2, i)
                read(token(6), *) stoc(3, i)
                read(token(7), *) stoc(4, i)
                call getso(socounter, i, k)
            else if (spline) then
                call getso(socounter, i, k)
            end if
            call assert(ngto(i)<=cnmax, &
                    "readGaussianBasis: GTO contraction length too large")
            do j = 1, ngto(i)
                !           !!read(iu,*) cntrctn(1,j,i),cntrctn(2,j,i)
                read(lines(idx + j), *) cntrctn(1, j, i), cntrctn(2, j, i)
            end do
            idx = idx + ngto(i) + 1
            nbasf = nbasf + 1

            if (bl(i) == 'S') then
                nbas = nbas + 1
                AOIdx%nums(k) = AOIdx%nums(k) + 1
                AOIdx%s(k, AOIdx%nums(k)) = nbas
                bn(i) = 1
            else if (bl(i) == 'P') then
                nbas = nbas + 3
                AOIdx%nump(k) = AOIdx%nump(k) + 1
                AOIdx%p(k, AOIdx%nump(k)) = nbas - 2   ! first idx of px,py,pz
                bn(i) = 2
            else if (bl(i) == 'D') then
                nbas = nbas + 6
                AOIdx%numd(k) = AOIdx%numd(k) + 1
                AOIdx%d(k, AOIdx%numd(k)) = nbas - 5
                bn(i) = 3
            else if (bl(i) == 'F') then
                nbas = nbas + 10
                AOIdx%numf(k) = AOIdx%numf(k) + 1
                AOIdx%f(k, AOIdx%numf(k)) = nbas - 9
                bn(i) = 4
            else if (bl(i) == 'G') then
                nbas = nbas + 15
                AOIdx%numg(k) = AOIdx%numg(k) + 1
                AOIdx%f(k, AOIdx%numg(k)) = nbas - 14
                bn(i) = 5
            end if
            bc(i) = k
        end do
        if (i>basmax) iflag = 1
        idx = idx + 1

    end subroutine readGaussianBasis


    !     ----------------------------
    subroutine calcNormFactors()
        !     ----------------------------

        !     * Calculate Normalizations factors *

        !     // The MO coefficients in Gaussian or GAMESS pertain to
        !     // normalized basis functions (=contracted primitive Gaussians)

        !     // Calculating Normalization Constants For Cartesian AO's
        !     // both GTO and STO
        !     //   norm(i) = radial norm(i) * angular norm(i)
        !     // for GTO, the norm is directly multiplied with the contraction coeff
        !     // for d and f, the angular norm for d_xx, f_xxx, resp. is used
        !     // in eloc and mdet, the norm for the other d,f AO's are corrected
        integer i, j, k
        real(r8) tmp, fak, alp, n

        if (.not. normalize) return
        normalize = .false.

        if (logmode >= 3) write(iul, *) "Normalization factors of basis functions:"
        do i = 1, nbasf
            !        // Angular norm factor
            if (bl(i) == 'S') then
                tmp = 1d0 / (2d0 * sqrt(pi))
            else if (bl(i) == 'P') then
                tmp = sqrt(3d0 / pi) / 2d0
            else if (bl(i) == 'D') then
                tmp = sqrt(5d0 / pi) / 2d0       ! for d_xx,d_yy,d_zz only!
            else if (bl(i) == 'F') then
                tmp = sqrt(7d0 / pi) / 2d0       ! for f_xxx only!
            else if (bl(i) == 'G') then
                tmp = sqrt(9d0 / pi) / 2d0       ! for g_xxxx only!
            else
                call abortp('(calcNormFactors): l>3 not yet implemented')
            endif
            !        // Radial norm factor
            if (typ(i) == 'STO') then
                fak = 1d0
                do j = 2, 2 * bn(i)
                    fak = fak * j
                enddo
                norm(i) = tmp * (2d0 * bzet(i))**(bn(i) + 0.5d0) / sqrt(fak)
            else
                !           ! Check if it's already normalized
                n = calcGTONorm(i)
                if (logmode >= 3) then
                    write(iul, '(A7,I3,A,G15.8)') "Before ", i, ": ", n
                endif
                !           ! Check if the GTO is already normalized
                if(abs(n - 1d0) < 1d-5) cycle

                !           ! GTO
                fak = 1d0                                   ! get (2n-1)!!
                do k = 3, 2 * bn(i) - 1, 2
                    fak = fak * k
                enddo
                do j = 1, ngto(i)
                    alp = cntrctn(1, j, i)
                    cntrctn(2, j, i) = cntrctn(2, j, i) * tmp * &
                            sqrt(2**(bn(i) + 1) / fak * sqrt((2d0 * alp)**(2 * bn(i) + 1) / pi))
                enddo

                !           ! Calculate norm factor again
                n = calcGTONorm(i)
                if (logmode >= 2 .and. abs(n - 1d0) >= 1d-5) then
                    write(iul, '(A,I3,A,G15.8)') "Warning: CGTO #", i, " multiplied by norm factor ", n
                endif
                if (logmode >= 3) then
                    write(iul, '(A7,I3,A,G15.8)') "After ", i, ": ", n
                endif
                cntrctn(2, :, i) = cntrctn(2, :, i) * n

                n = calcGTONorm(i)
                call assert(abs(n - 1d0) < 1d-5, "GTO normalization failed")
            endif
        enddo

    end subroutine calcNormFactors

    !     -----------------------------------------------------------
    real(r8) function calcGTONorm(i)
        !     -----------------------------------------------------------
        integer, intent(in) :: i ! GTO index

        integer :: j, k
        real(r8) :: fak, sum, n

        fak = 1d0
        do k = 3, 2 * (bn(i) - 1) - 1, 2
            fak = fak * k
        enddo

        sum = 0d0
        do j = 1, ngto(i)
            do k = 1, ngto(i)
                sum = sum + cntrctn(2, j, i) * cntrctn(2, k, i) / &
                        (cntrctn(1, j, i) + cntrctn(1, k, i))**(bn(i) + 0.5d0)
            enddo
        enddo
        n = (pi**1.5d0 * fak * sum / 2**(bn(i) - 1)) ** (-0.5d0)

        calcGTONorm = n
    end function calcGTONorm


    !     -----------------------------------------------------------
    subroutine getso(socounter, i, k)
        !     -----------------------------------------------------------

        integer, intent(in) :: i, k
        integer, intent(inout) :: socounter
        integer, save :: k_old, sonew, k_now

        if (ngto(i) == 1) then
            so(i) = 0
        else if (ngto(i) > 1) then
            if (k == 1) then
                socounter = socounter + 1
                sonew = sonew + 1
                so(i) = socounter
                k_now = k
            else if (k > 1) then
                if (basis == 'diff' .or. basis == 'default') then
                    if (k /= k_now .and. atoms(k)%elemIdx==atoms(k_now)%elemIdx&
                            .and. atoms(k)%ba == atoms(k_now)%ba) then
                        socounter = socounter - sonew
                        sonew = 1
                        socounter = socounter + 1
                        so(i) = socounter
                        k_now = k
                        k_old = k - 1
                    else if (k == k_now .and. k /= k_old) then
                        socounter = socounter + 1
                        sonew = sonew + 1
                        so(i) = socounter
                    else if (k /= k_now&
                            .and. atoms(k)%elemIdx/=atoms(k_now)%elemIdx) then
                        sonew = sonew - sonew + 1
                        socounter = socounter + 1
                        so(i) = socounter
                        k_now = k
                        k_old = k - 1
                    else if (k /= k_now&
                            .and. atoms(k)%ba /= atoms(k_now)%ba) then
                        sonew = sonew - sonew + 1
                        socounter = socounter + 1
                        so(i) = socounter
                        k_now = k
                        k_old = k - 1
                    end if
                else if (basis /= 'diff' .and. basis /= 'default') then
                    if (k /= k_now&
                            .and. atoms(k)%elemIdx==atoms(k_now)%elemIdx) then
                        socounter = socounter - sonew
                        sonew = 1
                        socounter = socounter + 1
                        so(i) = socounter
                        k_now = k
                        k_old = k - 1
                    else if (k == k_now .and. k /= k_old) then
                        socounter = socounter + 1
                        sonew = sonew + 1
                        so(i) = socounter
                    else if (k /= k_now&
                            .and. atoms(k)%elemIdx/=atoms(k_now)%elemIdx) then
                        sonew = sonew - sonew + 1
                        socounter = socounter + 1
                        so(i) = socounter
                        k_now = k
                        k_old = k - 1
                    end if

                endif
            endif
        endif

    end subroutine getso

END MODULE aosData_m
