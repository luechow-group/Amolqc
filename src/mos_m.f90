! Copyright (C) 1996, 2005-2007, 2012 Arne Luechow
! Copyright (C) 2000 Sebastian Manten
!
! SPDX-License-Identifier: GPL-3.0-or-later

!     ----------
MODULE mos_m
    !     ----------

    !  Module for MO calculation

    ! this version for array implementation of dense matrices
    ! (needs reorganisation)
    !
    ! note alternative direct AO/MO calculation !
    ! note sparse matrix list implementation in molist
    !
    ! $Id: mos_m.f,v 1.1.1.1 2007/04/25 13:42:20 luechow Exp $
    !
    ! $Log: mos_m.f,v $
    ! Revision 1.1.1.1  2007/04/25 13:42:20  luechow
    ! QMC program amolqc. rewritten in 2006. AL
    !
    ! Revision 1.9  2006/04/28 13:05:30  annika
    ! needed for regional analysis
    !
    ! Revision 1.8  2006/04/20 19:23:49  luechow
    ! gaussFOrder now in all ao routines. walkerStatistics added. trajectory option
    ! added. Note: new $ecp block in .in and .wf.
    !
    ! Revision 1.7  2005/03/11 13:33:12  luechow
    ! adding 'linscal' routines using sparse matrix algorithms with umfpack lib
    !
    !
    use kinds_m, only : r8
    use error_m
    use wfData_m
    use aosData_m, only : typ, bl, uao, uxao, uyao, uzao, u2ao, mAOElecConfigs
    implicit none

    real(r8), allocatable :: cmo(:, :)            ! MO coefficients
    real(r8), allocatable :: mat(:, :, :), mat1x(:, :, :), mat1y(:, :, :), &
            mat1z(:, :, :), mat2(:, :, :)

    !cc needed for aomo_calc; remains unallocated if AO/MO evaluation is done
    !cc by aocalc and mocalc
    real(r8), allocatable :: cmoa(:)   !cmo array used for simultanious AO/MO evaluation
    integer, allocatable :: nmos(:), mo_o(:)

    integer :: mMOElecConfigs = 0

CONTAINS

    !     ------------------------------
    subroutine mos_initialize(nec)
        !     ------------------------------

        integer, intent(in), optional :: nec
        integer nn, alstat
        if (present(nec)) then
            mMOElecConfigs = nec
        else
            mMOElecConfigs = 1
        endif
        nn = mMOElecConfigs

        call assert(nbas>0 .and. ne>0, &
                'aos_initialize: nbas or ne not yet set')

        if (allocated(mat)) then
            if (size(mat, 1)/=norb .or. size(mat, 2)/=ne&
                    .or. size(mat, 3) /= nn) then
                call mos_deallocate()
            else
                return
            endif
        endif

        call assert(norb>0 .and. ne>0, &
                'mos_initialize failed: sizes not set')
        allocate(mat(norb, ne, nn), mat1x(norb, ne, nn), mat1y(norb, ne, nn), &
                mat1z(norb, ne, nn), mat2(norb, ne, nn), stat = alstat)
        call assert(alstat==0, 'mos_initialize: allocation failed')

    end subroutine mos_initialize

    !     ---------------------------
    subroutine mos_deallocate()
        !     ---------------------------

        integer alstat
        if (allocated(mat)) then
            deallocate(mat, mat1x, mat1y, mat1z, mat2, stat = alstat)
            call assert(alstat==0, 'aos_deallocate: failed')
        endif
    end subroutine mos_deallocate


    !     ----------------------------
    subroutine moinput(lines)
        !     ----------------------------

        ! moinput reads MO related input (the MO coeffs in cmo) from file unit iu.
        !
        character(len = *), intent(in) :: lines(:)! input lines
        integer i, ii, j, jj, alstat, idx
        character fmt*19

        !     // Read in number of MO's
        read(lines(2), *) norb

        call assert(.not. allocated(cmo), &
                'moinput: mo coeffs already allocated')
        call assert(nbas>0 .and. norb>0, &
                'mos_initialize failed: sizes not set')
        allocate(cmo(nbas, norb), stat = alstat)
        call assert(alstat==0, 'moinput: allocation failed')

        call mos_initialize()

        !     // The orbital coefficients may be given directly in the formats
        !     // of Gaussian (punch=MO), or GAMESS (from the $VEC section
        !     // of the .dat file), or in Free Format
        !     // 'tmx' refers to the tm2xx output (turbomole)
        !     // note the different order of f functions in Gaussian and Turbomole/Gamess
        idx = 4
        if (evfmt=='tmx' .or. evfmt=='gau') then
            fmt = '(5D15.8)'
            do i = 1, norb
                read(lines(idx), *) ii
                jj = 0; idx = idx + 1
                do
                    read(lines(idx), fmt, pad='no') (cmo(j, ii), j = jj + 1, min(jj + 5, nbas))
                    jj = jj + 5; idx = idx + 1
                    if (jj>=nbas) exit
                end do
            end do
        else if (evfmt .eq. 'gms') then
            fmt = '(5X,5E15.8)'
            do i = 1, norb
                jj = 0
                do
                    read(lines(idx), fmt, pad='no') (cmo(j, i), j = jj + 1, min(jj + 5, nbas))
                    jj = jj + 5; idx = idx + 1
                    if (jj>=nbas) exit
                end do
            end do
        else if (evfmt .eq. 'fre' .or. evfmt .eq. 'mol') then
            do i = 1, norb
                read(lines(idx), *) ii
                jj = 0; idx = idx + 1
                do
                    read(lines(idx), *, pad='no') (cmo(j, ii), j = jj + 1, min(jj + 5, nbas))
                    jj = jj + 5; idx = idx + 1
                    if (jj>=nbas) exit
                end do
            enddo
        else
            call abortp('(moinput): unrecognized evfmt')
        endif

        !cccccccc

        if (aomocomb) then
            do j = 1, norb
                if (typ(j)=='STO') then
                    call abortp('aomo_calc does not support STOs yet')
                endif
            enddo
            !cc
            if (logmode>=2) then
                write(iul, *)
                write(iul, *) 'AO/MO evaluation carried', &
                        ' out using aomo_calc.'
                write(iul, *)
            endif
            if (cutmo) then
                call mocut(mocutoff)
            else
                call reorder_cmo()
            endif
        else
            if (logmode>=2) then
                write(iul, *)
                write(iul, *) 'AO/MO evaluation carried', &
                        ' out using aocalc and mocalc.'
                write(iul, *)
            endif
        endif

    end subroutine moinput

    !=================================================

    !     ---------------------
    subroutine mocut(thr)
        !     ---------------------

        ! 1) reorder MO coefficients --- afterwords they are seemingly completely messed up...
        ! but in the correct order for aomo_calc
        ! 2) apply MO cutoff. (Simplest way to handle localized orbitals; possibly also efficient
        ! for canonical orbitals (screening for near zeros)).

        real(r8) :: thr
        real(r8) :: tmp

        integer :: bf, al, d, nd, moc
        integer :: alstat
        integer :: cnt, j

        allocate(cmoa(nbas * norb), stat = alstat)
        allocate(nmos(nbas), mo_o(nbas * norb), stat = alstat)
        if (alstat.ne.0) call abortp('(moinput):allocation error')

        al = 1
        moc = 0

        do bf = 1, nbasf
            if (bl(bf).eq.'S') then
                nd = 0
            elseif (bl(bf).eq.'P') then
                nd = 2
            elseif (bl(bf).eq.'D') then
                nd = 5
            elseif (bl(bf).eq.'F') then
                nd = 9
            else
                call abortp('(getaos): wrong GTO')
            endif
            do d = 0, nd
                cnt = 0
                do j = 1, norb
                    tmp = cmo(al + d, j)
                    if (abs(tmp).gt.thr) then
                        moc = moc + 1
                        cnt = cnt + 1
                        mo_o(moc) = j
                        cmoa(moc) = tmp
                    endif
                enddo
                nmos(al + d) = cnt
            enddo
            al = al + 1 + nd
        enddo

        !cc some output

        write(iul, *)
        write(iul, *) 'overall number of MO-coefficients (nbas*norb): ', &
                nbas * norb
        write(iul, '('' MO-cutoff: '',G8.2)') thr
        write(iul, *) 'number of MO-coefficients after applying cutoff: ', &
                moc
        write(iul, *)

    end subroutine mocut


    !     ------------------------
    subroutine reorder_cmo()
        !     ------------------------

        integer bf, al, d, nd, moc
        integer alstat, j

        allocate(cmoa(nbas * norb), stat = alstat)
        if (alstat.ne.0) call abortp('(moinput):allocation error')

        !cc reorder MO coefficients --- afterwords they are seemingly completely messed up...
        !cc but in the correct order for aomo_calc

        al = 1
        moc = 0

        do bf = 1, nbasf
            if (bl(bf).eq.'S') then
                nd = 0
            elseif (bl(bf).eq.'P') then
                nd = 2
            elseif (bl(bf).eq.'D') then
                nd = 5
            elseif (bl(bf).eq.'F') then
                nd = 9
            else
                call abortp('(getaos): wrong GTO')
            endif
            do j = 1, norb
                do d = 0, nd
                    moc = moc + 1
                    cmoa(moc) = cmo(al + d, j)
                enddo
            enddo
            al = al + 1 + nd
        enddo

    end subroutine reorder_cmo

    !=================================================
    !
    !     -----------------------
    subroutine mooutput(iu)
        !     -----------------------
        !
        ! mooutput writes MOs to file unit iu.
        !
        integer iu
        integer i, j, bidx, bidxmax, max, min
        character fmt*19

        write(iu, '(i4)') norb
        write(iu, *)
        if (evfmt=='tmx' .or. evfmt=='gau') then
            fmt = '(5D15.8)'
            do i = 1, norb
                write(iu, '(i4)') i
                write(iu, fmt) (cmo(j, i), j = 1, nbas)
            enddo
        else if (evfmt .eq. 'gms') then
            fmt = '(I2,I3,1P,5E15.8)'
            bidxmax = CEILING(nbas / real(5))
            bidx = 1
            max = 0
            i = 1
            do
                if (bidx>bidxmax) bidx = 1
                min = max + 1
                max = max + 5
                if (max.GE.nbas) then
                    max = nbas
                endif
                write(iu, fmt) i, bidx, (cmo(j, i), j = min, max)
                bidx = bidx + 1
                if(max==nbas) then
                    i = i + 1
                    bidx = 1
                    max = 0
                endif
                if(i>norb) exit
            enddo
        else if (evfmt .eq. 'fre' .or. evfmt=='mol') then
            do i = 1, norb
                write(iu, *) i
                write(iu, *) (cmo(j, i), j = 1, nbas)
            enddo
        else
            call abortp('(mooutput): unrecognized evfmt')
        endif

    end subroutine mooutput

    !================================================

    !     ---------------------
    subroutine mocalc(ie)
        !     ---------------------

        ! mocalc calculates the MO's (+ derivatives) for one (ie) or
        ! all (ie=0) electrons. mocalc accesses the AO arrays uao etc.
        ! assumed to be up to date (e.g. by previous call to aocalc)

        ! input parameter:
        integer ie                      ! electron to calculate MO's for (0=all)
        ! variables:
        integer i, j, n
        real(r8)  zero, one
        character*1 uplo1, uplo2

        !     ! Calculate all MO's and derivatives for all or one electron(s)
        !     ! for all electron configurations

        call assert(size(mat, 3)>=mAOElecConfigs, &
                '(mocalc): sizes do not match')
        mMOElecConfigs = mAOElecConfigs

        do n = 1, mMOElecConfigs

            if (ie == 0) then                 ! all electrons new
                if (.not. useLAlib_mos) then
                    do i = 1, ne
                        do j = 1, norb
                            mat(j, i, n) = dot_product(cmo(:, j), uao(:, i, n))
                            mat1x(j, i, n) = dot_product(cmo(:, j), uxao(:, i, n))
                            mat1y(j, i, n) = dot_product(cmo(:, j), uyao(:, i, n))
                            mat1z(j, i, n) = dot_product(cmo(:, j), uzao(:, i, n))
                            mat2(j, i, n) = dot_product(cmo(:, j), u2ao(:, i, n))
                        enddo
                    enddo
                else
                    !
                    !       Perform Matrix multiplication by BLAS routine
                    !
                    one = 1._r8
                    zero = 0._r8
                    uplo1 = 'T'
                    uplo2 = 'N'
                    call dgemm(uplo1, uplo2, norb, ne, nbas, one, cmo, nbas, &
                            uao(:, :, n), nbas, zero, mat(:, :, n), norb)
                    call dgemm(uplo1, uplo2, norb, ne, nbas, one, cmo, nbas, &
                            uxao(:, :, n), nbas, zero, mat1x(:, :, n), norb)
                    call dgemm(uplo1, uplo2, norb, ne, nbas, one, cmo, nbas, &
                            uyao(:, :, n), nbas, zero, mat1y(:, :, n), norb)
                    call dgemm(uplo1, uplo2, norb, ne, nbas, one, cmo, nbas, &
                            uzao(:, :, n), nbas, zero, mat1z(:, :, n), norb)
                    call dgemm(uplo1, uplo2, norb, ne, nbas, one, cmo, nbas, &
                            u2ao(:, :, n), nbas, zero, mat2(:, :, n), norb)
                endif

            else
                i = ie
                if (.not. useLAlib_mos) then
                    do j = 1, norb
                        mat(j, i, n) = dot_product(cmo(:, j), uao(:, i, n))
                        mat1x(j, i, n) = dot_product(cmo(:, j), uxao(:, i, n))
                        mat1y(j, i, n) = dot_product(cmo(:, j), uyao(:, i, n))
                        mat1z(j, i, n) = dot_product(cmo(:, j), uzao(:, i, n))
                        mat2(j, i, n) = dot_product(cmo(:, j), u2ao(:, i, n))
                    enddo
                else
                    one = 1._r8
                    zero = 0._r8
                    uplo1 = 'T'
                    call DGEMV(uplo1,nbas,norb,one,cmo(:,:),nbas,uao(:, i, n),1,zero,mat(:, i, n),1)
                    call DGEMV(uplo1,nbas,norb,one,cmo(:,:),nbas,uxao(:, i, n),1,zero,mat1x(:, i, n),1)
                    call DGEMV(uplo1,nbas,norb,one,cmo(:,:),nbas,uyao(:, i, n),1,zero,mat1y(:, i, n),1)
                    call DGEMV(uplo1,nbas,norb,one,cmo(:,:),nbas,uzao(:, i, n),1,zero,mat1z(:, i, n),1)
                    call DGEMV(uplo1,nbas,norb,one,cmo(:,:),nbas,u2ao(:, i, n),1,zero,mat2(:, i, n),1)
                end if
            endif

        enddo

    end subroutine mocalc


    !================================================

    !     ----------------------
    subroutine mo1calc(ie)
        !     ----------------------

        ! mocalc calculates the MO's for one (ie) or all (ie=0) electrons
        ! mocalc accesses the AO arrays uao etc. assumed to be up to date
        ! (e.g. by previous call to aocalc)

        ! input parameter:
        integer ie            ! electron to calculate MO's for (0=all)
        ! variables:
        integer i, j
        real(r8)  zero, one
        character*1 uplo1, uplo2

        mMOElecConfigs = 1

        !     ! Calculate all MO's and derivatives for all or one electron(s)

        if (ie == 0) then                 ! all electrons new
            if (.not. useLAlib_mos) then
                do i = 1, ne
                    do j = 1, norb
                        mat(j, i, 1) = dot_product(cmo(:, j), uao(:, i, 1))
                    enddo
                enddo
            else
                !
                !       Perform Matrix multiplication by BLAS routine
                !
                one = 1._r8
                zero = 0._r8
                uplo1 = 'T'
                uplo2 = 'N'
                call dgemm(uplo1, uplo2, norb, ne, nbas, one, cmo(:,:), nbas, &
                        uao(:,:,1), nbas, zero, mat(:,:,1), norb)
            endif

        else
            i = ie
            if (.not. useLAlib_mos) then
                do j = 1, norb
                    mat(j, i, 1) = dot_product(cmo(:, j), uao(:, i, 1))
                enddo
            else
                one = 1._r8
                zero = 0._r8
                uplo1 = 'T'
                call DGEMV(uplo1,nbas,norb,one,cmo(:,:),nbas,uao(:, i, 1),1,zero,mat(:, i, 1),1)
            end if
        endif

    end subroutine mo1calc

END MODULE mos_m

