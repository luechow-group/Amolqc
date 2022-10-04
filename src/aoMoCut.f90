! Copyright (C) 2004 Christian Diedrich
!
! SPDX-License-Identifier: GPL-3.0-or-later

!cc aomocut_m.f90 contains the routines which are used if AO- and MO-cutoff are desired.
!cc In principle this routine could also be used for calculations without any cutoff
!cc or with AO-cutoff only but the routines in aomo_m.f90 are faster for that purpose.
!cc !!! WARNING !!!
!cc The ordering of the MO-coefficients in 'cmoa' is different depending on
!cc which routines are used. 'aomocut_calc' requires the cmoa array to be initialized by
!cc 'mocut' (in mos_m.f) !!!
!cc !!! WARNING !!!

subroutine aomocut_calc(ie, x, y, z, rai)

    use kinds_m, only : r8
    use wfData_m
    use aoMo_m

    implicit none

    ! input parameters:
    integer ie                 ! if >0 only AO's for electron ie recalculated

    real(r8), dimension(nmax) :: x, y, z    ! x,y,z coordinates of position vector
    real(r8), dimension(amax, nmax) :: rai ! r_ai electron-nucleus distances

    ! constants:
    real(r8) :: sqr3, sqr5
    parameter (sqr3 = 1.73205080756887729d0, sqr5 = 2.236067977499789696d0)
    ! local variables
    integer :: bf, a, i, i1, i2, ii, ic, al
    integer :: j, d, moc, n
    real(r8) :: xx, yy, zz, rr, r2, alp, nrm, u, ux, dx, dy, dz, tmp, &
            dx2, dy2, dz2, dxyz, dxdy, dxdz, dydz
    real(r8) :: dy2dx, dx2dy, dx2dz, dz2dx, dy2dz, dz2dy, dxdydz
    !cc
    real(r8), dimension(5) :: tmps !second index means:
    real(r8), dimension(0:2, 5) :: tmpp !1:   AO itsself
    real(r8), dimension(0:5, 5) :: tmpd !2-4: derivative with respect to x,y,z
    real(r8), dimension(0:9, 5) :: tmpf !5:   laplacian
    logical gaussFOrder       ! .t.: Gaussian order for f function
    !                               ! .f.: Gamess==Turbomole order used



    ! bf refers to the degenerate set of cartesian
    ! basis function (S:1,P:3,D:6,F:10) as input, which may be of type STO
    ! or contracted GTO.
    ! al refers to the individual basis function, as used in LCAO-MO's.
    ! (composed in subroutine mdetwf)
    ! i refers to the current electron.

    !-----Calculation of the AO's and their derivatives

    if (evfmt=='gau' .or. evfmt=='mol') then
        gaussFOrder = .true.
    else
        gaussFOrder = .false.
    endif

    if (ie == 0) then                     ! AO's for all electrons
        i1 = 1
        i2 = ne
    else
        i1 = ie                              ! only AO for electron ie
        i2 = ie
    endif

    do i = i1, i2                              ! loop over electrons
        xx = x(i)
        yy = y(i)
        zz = z(i)

        al = 1                               !Pointer for individual basis function
        moc = 0                              !Pointer for cmoa array

        ! initialisation

        mat(1:norb, i, 1) = 0d0
        mat1x(1:norb, i, 1) = 0d0
        mat1y(1:norb, i, 1) = 0d0
        mat1z(1:norb, i, 1) = 0d0
        mat2(1:norb, i, 1) = 0d0

        do bf = 1, nbasf                        ! loop over basis functions
            a = bc(bf)                         ! center of AO
            rr = rai(a, i)                      ! r_ai

            if (cutao) then                    !AO - Cutoff
                if (rr>aocuts(bf)) then       ! --> do nothing but adjusting the counters
                    if (bl(bf) == 'S') then
                        moc = moc + nmos(al)
                        al = al + 1
                    elseif (bl(bf) == 'P') then
                        do n = 0, 2
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 3
                    elseif (bl(bf) == 'D') then
                        do n = 0, 5
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 6
                    elseif (bl(bf) == 'F') then
                        do n = 0, 9
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 10
                    else
                        call abortp('(getaos): wrong GTO')
                    endif
                    cycle  !continue with next basis function
                endif
            endif

            r2 = rr * rr

            if (bl(bf) == 'S') then                 ! 1s GTO

                dx = xx - atoms(a)%cx
                dy = yy - atoms(a)%cy
                dz = zz - atoms(a)%cz

                tmps = 0d0

                do ic = 1, ngto(bf)                     ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)
                    ux = -2d0 * alp * u

                    tmps(1) = tmps(1) + u
                    tmps(2) = tmps(2) + ux * dx
                    tmps(3) = tmps(3) + ux * dy
                    tmps(4) = tmps(4) + ux * dz
                    tmps(5) = tmps(5) + ux * (3d0 - 2d0 * alp * r2)
                enddo
                !cc MO calculation
                do n = 1, nmos(al)
                    moc = moc + 1
                    tmp = cmoa(moc)
                    j = mo_o(moc)
                    mat(j, i, 1) = mat(j, i, 1) + tmp * tmps(1)
                    mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmps(2)
                    mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmps(3)
                    mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmps(4)
                    mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmps(5)
                enddo

                al = al + 1

            else if (bl(bf) == 'P') then             ! 2p GTO's
                !              // do all 3 P simultaneously (same exponent is required)
                !              // order p_x,p_y,p_z

                dx = xx - atoms(a)%cx
                dx2 = dx * dx
                dy = yy - atoms(a)%cy
                dy2 = dy * dy
                dz = zz - atoms(a)%cz
                dz2 = dz * dz

                dxdy = dx * dy
                dxdz = dx * dz
                dydz = dy * dz

                tmpp = 0d0

                do ic = 1, ngto(bf)                      ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)
                    ux = -2d0 * alp * u

                    tmpp(0, 1) = tmpp(0, 1) + dx * u
                    tmpp(1, 1) = tmpp(1, 1) + dy * u
                    tmpp(2, 1) = tmpp(2, 1) + dz * u

                    tmpp(0, 2) = tmpp(0, 2) + u + ux * dx2
                    tmpp(1, 2) = tmpp(1, 2) + ux * dxdy
                    tmpp(2, 2) = tmpp(2, 2) + ux * dxdz
                    tmpp(0, 3) = tmpp(0, 3) + ux * dxdy
                    tmpp(1, 3) = tmpp(1, 3) + u + ux * dy2
                    tmpp(2, 3) = tmpp(2, 3) + ux * dydz
                    tmpp(0, 4) = tmpp(0, 4) + ux * dxdz
                    tmpp(1, 4) = tmpp(1, 4) + ux * dydz
                    tmpp(2, 4) = tmpp(2, 4) + u + ux * dz2

                    tmp = (5d0 - 2d0 * alp * r2) * ux
                    tmpp(0, 5) = tmpp(0, 5) + tmp * dx
                    tmpp(1, 5) = tmpp(1, 5) + tmp * dy
                    tmpp(2, 5) = tmpp(2, 5) + tmp * dz
                enddo

                !cc MO calculation
                do d = 0, 2
                    do n = 1, nmos(al + d)
                        moc = moc + 1
                        tmp = cmoa(moc)
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + tmp * tmpp(d, 1)
                        mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpp(d, 2)
                        mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpp(d, 3)
                        mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpp(d, 4)
                        mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpp(d, 5)
                    enddo
                enddo

                al = al + 3

            else if (bl(bf) == 'D') then         ! 3d GTO
                !              // do all 6 D simultaneously (same exponent is required)
                !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)

                dx = xx - atoms(a)%cx
                dx2 = dx * dx
                dy = yy - atoms(a)%cy
                dy2 = dy * dy
                dz = zz - atoms(a)%cz
                dz2 = dz * dz

                dxdy = dx * dy
                dxdydz = dxdy * dz
                dy2dx = dxdy * dy
                dx2dy = dxdy * dx
                dxdz = dx * dz
                dx2dz = dxdz * dx
                dz2dx = dxdz * dz
                dydz = dy * dz
                dy2dz = dydz * dy
                dz2dy = dydz * dz

                tmpd = 0d0

                do ic = 1, ngto(bf)                      ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)
                    ux = -2d0 * alp * u

                    tmpd(0, 1) = tmpd(0, 1) + dx2 * u
                    tmpd(1, 1) = tmpd(1, 1) + dy2 * u
                    tmpd(2, 1) = tmpd(2, 1) + dz2 * u

                    tmpd(0, 2) = tmpd(0, 2) + (2d0 * u + ux * dx2) * dx
                    tmpd(1, 2) = tmpd(1, 2) + dy2dx * ux
                    tmpd(2, 2) = tmpd(2, 2) + dz2dx * ux
                    tmpd(0, 3) = tmpd(0, 3) + dx2dy * ux
                    tmpd(1, 3) = tmpd(1, 3) + (2d0 * u + ux * dy2) * dy
                    tmpd(2, 3) = tmpd(2, 3) + dz2dy * ux
                    tmpd(0, 4) = tmpd(0, 4) + dx2dz * ux
                    tmpd(1, 4) = tmpd(1, 4) + dy2dz * ux
                    tmpd(2, 4) = tmpd(2, 4) + (2d0 * u + ux * dz2) * dz
                    tmp = (7d0 - 2d0 * alp * r2) * ux
                    tmpd(0, 5) = tmpd(0, 5) + 2d0 * u + dx2 * tmp
                    tmpd(1, 5) = tmpd(1, 5) + 2d0 * u + dy2 * tmp
                    tmpd(2, 5) = tmpd(2, 5) + 2d0 * u + dz2 * tmp

                    u = sqr3 * u                   ! correction of norm
                    ux = sqr3 * ux                 ! N(dxx)*sqr3 = N(dxy)

                    tmpd(3, 1) = tmpd(3, 1) + dxdy * u
                    tmpd(4, 1) = tmpd(4, 1) + dxdz * u
                    tmpd(5, 1) = tmpd(5, 1) + dydz * u
                    tmp = ux * dxdydz
                    tmpd(3, 2) = tmpd(3, 2) + (u + ux * dx2) * dy
                    tmpd(4, 2) = tmpd(4, 2) + (u + ux * dx2) * dz
                    tmpd(5, 2) = tmpd(5, 2) + tmp
                    tmpd(3, 3) = tmpd(3, 3) + (u + ux * dy2) * dx
                    tmpd(4, 3) = tmpd(4, 3) + tmp
                    tmpd(5, 3) = tmpd(5, 3) + (u + ux * dy2) * dz
                    tmpd(3, 4) = tmpd(3, 4) + tmp
                    tmpd(4, 4) = tmpd(4, 4) + (u + ux * dz2) * dx
                    tmpd(5, 4) = tmpd(5, 4) + (u + ux * dz2) * dy
                    tmp = (7d0 - 2d0 * alp * r2) * ux
                    tmpd(3, 5) = tmpd(3, 5) + tmp * dxdy
                    tmpd(4, 5) = tmpd(4, 5) + tmp * dxdz
                    tmpd(5, 5) = tmpd(5, 5) + tmp * dydz

                enddo

                !cc MO calculation

                do d = 0, 5
                    do n = 1, nmos(al + d)
                        moc = moc + 1
                        tmp = cmoa(moc)
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + tmp * tmpd(d, 1)
                        mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpd(d, 2)
                        mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpd(d, 3)
                        mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpd(d, 4)
                        mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpd(d, 5)
                    enddo
                enddo

                al = al + 6

            else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! f GTO
                !              // do all 10 F simultaneously (same exponent is required)
                !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

                tmpf = 0d0

                do ic = 1, ngto(bf)                      ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)
                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz
                    dxyz = dx * dy * dz
                    ux = -2d0 * alp * u

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0, 1) = tmpf(0, 1) + dx2 * dx * u
                    tmpf(1, 1) = tmpf(1, 1) + dy2 * dy * u
                    tmpf(2, 1) = tmpf(2, 1) + dz2 * dz * u

                    tmpf(0, 2) = tmpf(0, 2) + (3d0 * u + ux * dx2) * dx2
                    tmpf(1, 2) = tmpf(1, 2) + dy2 * dy * ux * dx
                    tmpf(2, 2) = tmpf(2, 2) + dz2 * dz * ux * dx
                    tmpf(0, 3) = tmpf(0, 3) + dx2 * dx * ux * dy
                    tmpf(1, 3) = tmpf(1, 3) + (3d0 * u + ux * dy2) * dy2
                    tmpf(2, 3) = tmpf(2, 3) + dz2 * dz * ux * dy
                    tmpf(0, 4) = tmpf(0, 4) + dx2 * dx * ux * dz
                    tmpf(1, 4) = tmpf(1, 4) + dy2 * dy * ux * dz
                    tmpf(2, 4) = tmpf(2, 4) + (3d0 * u + ux * dz2) * dz2
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(0, 5) = tmpf(0, 5) + (6d0 * u + dx2 * tmp) * dx
                    tmpf(1, 5) = tmpf(1, 5) + (6d0 * u + dy2 * tmp) * dy
                    tmpf(2, 5) = tmpf(2, 5) + (6d0 * u + dz2 * tmp) * dz

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    ux = sqr5 * ux                 ! N(fxxx)*sqrt(5) = N(fxxy)

                    tmpf(3, 1) = tmpf(3, 1) + dx2 * dy * u
                    tmpf(4, 1) = tmpf(4, 1) + dx2 * dz * u
                    tmpf(5, 1) = tmpf(5, 1) + dy2 * dx * u
                    tmpf(6, 1) = tmpf(6, 1) + dy2 * dz * u
                    tmpf(7, 1) = tmpf(7, 1) + dz2 * dx * u
                    tmpf(8, 1) = tmpf(8, 1) + dz2 * dy * u

                    ! derivatives
                    tmp = ux * dxyz
                    tmpf(3, 2) = tmpf(3, 2) + (2d0 * u + ux * dx2) * dx * dy
                    tmpf(4, 2) = tmpf(4, 2) + (2d0 * u + ux * dx2) * dx * dz
                    tmpf(5, 2) = tmpf(5, 2) + (u + ux * dx2) * dy2
                    tmpf(6, 2) = tmpf(6, 2) + tmp * dy
                    tmpf(7, 2) = tmpf(7, 2) + (u + ux * dx2) * dz2
                    tmpf(8, 2) = tmpf(8, 2) + tmp * dz
                    tmpf(3, 3) = tmpf(3, 3) + (u + ux * dy2) * dx2
                    tmpf(4, 3) = tmpf(4, 3) + tmp * dx
                    tmpf(5, 3) = tmpf(5, 3) + (2d0 * u + ux * dy2) * dx * dy
                    tmpf(6, 3) = tmpf(6, 3) + (2d0 * u + ux * dy2) * dy * dz
                    tmpf(7, 3) = tmpf(7, 3) + tmp * dz
                    tmpf(8, 3) = tmpf(8, 3) + (u + ux * dy2) * dz2
                    tmpf(3, 4) = tmpf(3, 4) + tmp * dx
                    tmpf(4, 4) = tmpf(4, 4) + (u + ux * dz2) * dx2
                    tmpf(5, 4) = tmpf(5, 4) + tmp * dy
                    tmpf(6, 4) = tmpf(6, 4) + (u + ux * dz2) * dy2
                    tmpf(7, 4) = tmpf(7, 4) + (2d0 * u + ux * dz2) * dx * dz
                    tmpf(8, 4) = tmpf(8, 4) + (2d0 * u + ux * dz2) * dy * dz
                    ! laplacians
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(3, 5) = tmpf(3, 5) + (2d0 * u + dx2 * tmp) * dy
                    tmpf(4, 5) = tmpf(4, 5) + (2d0 * u + dx2 * tmp) * dz
                    tmpf(5, 5) = tmpf(5, 5) + (2d0 * u + dy2 * tmp) * dx
                    tmpf(6, 5) = tmpf(6, 5) + (2d0 * u + dy2 * tmp) * dz
                    tmpf(7, 5) = tmpf(7, 5) + (2d0 * u + dz2 * tmp) * dx
                    tmpf(8, 5) = tmpf(8, 5) + (2d0 * u + dz2 * tmp) * dy

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    ux = sqr3 * ux                ! N(fxxx)*sqrt(15)=
                    !                                             ! N(fxxy)*sqrt(3)=N(fxyz)

                    tmpf(9, 1) = tmpf(9, 1) + dxyz * u

                    tmpf(9, 2) = tmpf(9, 2) + (u + ux * dx2) * dy * dz
                    tmpf(9, 3) = tmpf(9, 3) + (u + ux * dy2) * dx * dz
                    tmpf(9, 4) = tmpf(9, 4) + (u + ux * dz2) * dx * dy
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(9, 5) = tmpf(9, 5) + dxyz * tmp

                enddo

                !cc MO calculation
                do d = 0, 9
                    do n = 1, nmos(al + d)
                        moc = moc + 1
                        tmp = cmoa(moc)
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + tmp * tmpf(d, 1)
                        mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpf(d, 2)
                        mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpf(d, 3)
                        mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpf(d, 4)
                        mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpf(d, 5)
                    enddo
                enddo

                al = al + 10

            else if (bl(bf)=='F'.and.gaussFOrder) then     ! f GTO
                !              // do all 10 F simultaneously (same exponent is required)
                !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

                tmpf = 0d0

                do ic = 1, ngto(bf)                      ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)
                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz
                    dxyz = dx * dy * dz
                    ux = -2d0 * alp * u

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0, 1) = tmpf(0, 1) + dx2 * dx * u
                    tmpf(1, 1) = tmpf(1, 1) + dy2 * dy * u
                    tmpf(2, 1) = tmpf(2, 1) + dz2 * dz * u

                    tmpf(0, 2) = tmpf(0, 2) + (3d0 * u + ux * dx2) * dx2
                    tmpf(1, 2) = tmpf(1, 2) + dy2 * dy * ux * dx
                    tmpf(2, 2) = tmpf(2, 2) + dz2 * dz * ux * dx
                    tmpf(0, 3) = tmpf(0, 3) + dx2 * dx * ux * dy
                    tmpf(1, 3) = tmpf(1, 3) + (3d0 * u + ux * dy2) * dy2
                    tmpf(2, 3) = tmpf(2, 3) + dz2 * dz * ux * dy
                    tmpf(0, 4) = tmpf(0, 4) + dx2 * dx * ux * dz
                    tmpf(1, 4) = tmpf(1, 4) + dy2 * dy * ux * dz
                    tmpf(2, 4) = tmpf(2, 4) + (3d0 * u + ux * dz2) * dz2
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(0, 5) = tmpf(0, 5) + (6d0 * u + dx2 * tmp) * dx
                    tmpf(1, 5) = tmpf(1, 5) + (6d0 * u + dy2 * tmp) * dy
                    tmpf(2, 5) = tmpf(2, 5) + (6d0 * u + dz2 * tmp) * dz

                    !                 // f_xyy, f_xxy, f_xxz, f_xzz, f_yzz, f_yyz
                    u = sqr5 * u                   ! correction of norm
                    ux = sqr5 * ux                 ! N(fxxx)*sqrt(5) = N(fxxy)

                    tmpf(4, 1) = tmpf(4, 1) + dx2 * dy * u
                    tmpf(5, 1) = tmpf(5, 1) + dx2 * dz * u
                    tmpf(3, 1) = tmpf(3, 1) + dy2 * dx * u
                    tmpf(8, 1) = tmpf(8, 1) + dy2 * dz * u
                    tmpf(6, 1) = tmpf(6, 1) + dz2 * dx * u
                    tmpf(7, 1) = tmpf(7, 1) + dz2 * dy * u

                    ! derivatives
                    tmp = ux * dxyz
                    tmpf(4, 2) = tmpf(4, 2) + (2d0 * u + ux * dx2) * dx * dy
                    tmpf(5, 2) = tmpf(5, 2) + (2d0 * u + ux * dx2) * dx * dz
                    tmpf(3, 2) = tmpf(3, 2) + (u + ux * dx2) * dy2
                    tmpf(8, 2) = tmpf(8, 2) + tmp * dy
                    tmpf(6, 2) = tmpf(6, 2) + (u + ux * dx2) * dz2
                    tmpf(7, 2) = tmpf(7, 2) + tmp * dz
                    tmpf(4, 3) = tmpf(4, 3) + (u + ux * dy2) * dx2
                    tmpf(5, 3) = tmpf(5, 3) + tmp * dx
                    tmpf(3, 3) = tmpf(3, 3) + (2d0 * u + ux * dy2) * dx * dy
                    tmpf(8, 3) = tmpf(8, 3) + (2d0 * u + ux * dy2) * dy * dz
                    tmpf(6, 3) = tmpf(6, 3) + tmp * dz
                    tmpf(7, 3) = tmpf(7, 3) + (u + ux * dy2) * dz2
                    tmpf(4, 4) = tmpf(4, 4) + tmp * dx
                    tmpf(5, 4) = tmpf(5, 4) + (u + ux * dz2) * dx2
                    tmpf(3, 4) = tmpf(3, 4) + tmp * dy
                    tmpf(8, 4) = tmpf(8, 4) + (u + ux * dz2) * dy2
                    tmpf(6, 4) = tmpf(6, 4) + (2d0 * u + ux * dz2) * dx * dz
                    tmpf(7, 4) = tmpf(7, 4) + (2d0 * u + ux * dz2) * dy * dz
                    ! laplacians
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(4, 5) = tmpf(4, 5) + (2d0 * u + dx2 * tmp) * dy
                    tmpf(5, 5) = tmpf(5, 5) + (2d0 * u + dx2 * tmp) * dz
                    tmpf(3, 5) = tmpf(3, 5) + (2d0 * u + dy2 * tmp) * dx
                    tmpf(8, 5) = tmpf(8, 5) + (2d0 * u + dy2 * tmp) * dz
                    tmpf(6, 5) = tmpf(6, 5) + (2d0 * u + dz2 * tmp) * dx
                    tmpf(7, 5) = tmpf(7, 5) + (2d0 * u + dz2 * tmp) * dy

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    ux = sqr3 * ux                ! N(fxxx)*sqrt(15)=
                    !                                             ! N(fxxy)*sqrt(3)=N(fxyz)

                    tmpf(9, 1) = tmpf(9, 1) + dxyz * u

                    tmpf(9, 2) = tmpf(9, 2) + (u + ux * dx2) * dy * dz
                    tmpf(9, 3) = tmpf(9, 3) + (u + ux * dy2) * dx * dz
                    tmpf(9, 4) = tmpf(9, 4) + (u + ux * dz2) * dx * dy
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(9, 5) = tmpf(9, 5) + dxyz * tmp

                enddo

                !cc MO calculation
                do d = 0, 9
                    do n = 1, nmos(al + d)
                        moc = moc + 1
                        tmp = cmoa(moc)
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + tmp * tmpf(d, 1)
                        mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpf(d, 2)
                        mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpf(d, 3)
                        mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpf(d, 4)
                        mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpf(d, 5)
                    enddo
                enddo

                al = al + 10

            else
                call abortp('(getaos): wrong GTO')
            endif  ! bl
        enddo    ! bf-loop over basis functions
    enddo       ! i-loop over electrons

end subroutine aomocut_calc


subroutine aomocut1_calc(ie, x, y, z, rai)

    use wfData_m
    use aoMo_m

    implicit none

    integer ie                 ! if >0 only AO's for electron ie recalculated

    real(r8), dimension(nmax) :: x, y, z  ! x,y,z coordinates of position vector
    real(r8), dimension(amax, nmax) :: rai    ! r_ai electron-nucleus distances


    ! constants:
    real(r8) :: sqr3, sqr5
    parameter (sqr3 = 1.73205080756887729d0, sqr5 = 2.236067977499789696d0)
    ! variables
    integer :: bf, a, i, ii, i1, i2, ic, al
    integer :: j, d, moc, n
    real(r8) :: xx, yy, zz, rr, alp, nrm, u, dx, dy, dz, r2, dx2, dy2, dz2
    real(r8) :: dxdy, dxdz, dydz
    real(r8) :: tmp

    !cc
    real(r8) :: tmps
    real(r8), dimension(0:2) :: tmpp
    real(r8), dimension(0:5) :: tmpd
    real(r8), dimension(0:9) :: tmpf
    logical gaussFOrder       ! .t.: Gaussian order for f function
    !                               ! .f.: Gamess==Turbomole order used


    !-----Calculation of the AO's

    if (evfmt=='gau' .or. evfmt=='mol') then
        gaussFOrder = .true.
    else
        gaussFOrder = .false.
    endif

    if (ie == 0) then                     ! AO's for all electrons
        i1 = 1
        i2 = ne
    else
        i1 = ie                              ! only AO for electron ie
        i2 = ie
    endif

    do i = i1, i2                              ! loop over electrons
        xx = x(i)
        yy = y(i)
        zz = z(i)

        al = 1
        moc = 0                               !Pointer for cmoa array

        ! initialisation
        mat(1:norb, i, 1) = 0d0

        do bf = 1, nbasf                        ! loop over basis functions
            a = bc(bf)                         ! center of AO
            rr = rai(a, i)                      ! r_ai

            if (cutao) then                    !AO - Cutoff
                if (rr>aocuts(bf)) then       ! --> do nothing but adjusting the counters
                    if (bl(bf) == 'S') then
                        moc = moc + nmos(al)
                        al = al + 1
                    elseif (bl(bf) == 'P') then
                        do n = 0, 2
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 3
                    elseif (bl(bf) == 'D') then
                        do n = 0, 5
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 6
                    elseif (bl(bf) == 'F') then
                        do n = 0, 9
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 10
                    else
                        call abortp('(getaos): wrong GTO')
                    endif
                    cycle  !continue with next basis function
                endif
            endif

            r2 = rr * rr

            if (bl(bf) == 'S') then                 ! 1s GTO

                tmps = 0d0

                do ic = 1, ngto(bf) ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)

                    tmps = tmps + u
                enddo

                !cc MO calculation
                do n = 1, nmos(al)
                    moc = moc + 1
                    j = mo_o(moc)
                    mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmps
                enddo

                al = al + 1

            else if (bl(bf) == 'P') then             ! 2p GTO's
                !              // do all 3 P simultaneously (same exponent is required)
                !              // order p_x,p_y,p_z

                dx = xx - atoms(a)%cx
                dy = yy - atoms(a)%cy
                dz = zz - atoms(a)%cz

                tmpp = 0d0

                do ic = 1, ngto(bf)                      ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)

                    tmpp(0) = tmpp(0) + dx * u
                    tmpp(1) = tmpp(1) + dy * u
                    tmpp(2) = tmpp(2) + dz * u
                enddo

                !cc MO calculation

                do d = 0, 2
                    do n = 1, nmos(al + d)
                        moc = moc + 1
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpp(d)
                    enddo
                enddo

                al = al + 3

            else if (bl(bf) == 'D') then         ! 3d GTO
                !              // do all 6 D simultaneously (same exponent is required)
                !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)

                dx = xx - atoms(a)%cx
                dx2 = dx * dx
                dy = yy - atoms(a)%cy
                dy2 = dy * dy
                dz = zz - atoms(a)%cz
                dz2 = dz * dz

                dxdy = dx * dy
                dxdz = dx * dz
                dydz = dy * dz

                tmpd = 0d0

                do ic = 1, ngto(bf)                      ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)

                    tmpd(0) = tmpd(0) + dx2 * u
                    tmpd(1) = tmpd(1) + dy2 * u
                    tmpd(2) = tmpd(2) + dz2 * u

                    u = sqr3 * u                   ! correction of norm for last 3

                    tmpd(3) = tmpd(3) + dxdy * u
                    tmpd(4) = tmpd(4) + dxdz * u
                    tmpd(5) = tmpd(5) + dydz * u

                enddo

                !cc MO calculation
                do d = 0, 5
                    do n = 1, nmos(al + d)
                        moc = moc + 1
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpd(d)
                    enddo
                enddo

                al = al + 6

            else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! f GTO
                !              // do all 10 F simultaneously (same exponent is required)
                !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

                tmpf = 0d0

                do ic = 1, ngto(bf)                      ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)
                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0) = tmpf(0) + dx2 * dx * u
                    tmpf(1) = tmpf(1) + dy2 * dy * u
                    tmpf(2) = tmpf(2) + dz2 * dz * u

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm

                    tmpf(3) = tmpf(3) + dx2 * dy * u
                    tmpf(4) = tmpf(4) + dx2 * dz * u
                    tmpf(5) = tmpf(5) + dy2 * dx * u
                    tmpf(6) = tmpf(6) + dy2 * dz * u
                    tmpf(7) = tmpf(7) + dz2 * dx * u
                    tmpf(8) = tmpf(8) + dz2 * dy * u

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm

                    tmpf(9) = tmpf(9) + dx * dy * dz * u

                enddo
                do d = 0, 9
                    do n = 1, nmos(al + d)
                        moc = moc + 1
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpf(d)
                    enddo
                enddo

                al = al + 10

            else if (bl(bf)=='F'.and.gaussFOrder) then     ! f GTO
                !              // do all 10 F simultaneously (same exponent is required)
                !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

                tmpf = 0d0

                do ic = 1, ngto(bf)                      ! loop over contraction
                    alp = cntrctn(1, ic, bf)
                    u = cntrctn(2, ic, bf) * exp(-alp * r2)
                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0) = tmpf(0) + dx2 * dx * u
                    tmpf(1) = tmpf(1) + dy2 * dy * u
                    tmpf(2) = tmpf(2) + dz2 * dz * u

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm

                    tmpf(3) = tmpf(3) + dy2 * dx * u
                    tmpf(4) = tmpf(4) + dx2 * dy * u
                    tmpf(5) = tmpf(5) + dx2 * dz * u
                    tmpf(6) = tmpf(6) + dz2 * dx * u
                    tmpf(7) = tmpf(7) + dz2 * dy * u
                    tmpf(8) = tmpf(8) + dy2 * dz * u

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm

                    tmpf(9) = tmpf(9) + dx * dy * dz * u

                enddo
                do d = 0, 9
                    do n = 1, nmos(al + d)
                        moc = moc + 1
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpf(d)
                    enddo
                enddo

                al = al + 10

            else
                call abortp('(getaos): wrong GTO')
            endif  ! bl
        enddo     ! bf-loop over basis functions
    enddo        ! i-loop over electrons

end subroutine aomocut1_calc


!     ----------------------------------

subroutine aomocutspl_calc(ie, x, y, z, rai)

    use wfData_m
    use aoMo_m
    use cspline_m

    implicit none

    ! input parameters:
    integer ie                 ! if >0 only AO's for electron ie recalculated

    real(r8), dimension(nmax) :: x, y, z    ! x,y,z coordinates of position vector
    real(r8), dimension(amax, nmax) :: rai ! r_ai electron-nucleus distances

    ! constants:
    real(r8) :: sqr3, sqr5
    parameter (sqr3 = 1.73205080756887729d0, sqr5 = 2.236067977499789696d0)
    ! local variables
    integer :: bf, a, i, i1, i2, ii, ic, al
    integer :: j, d, moc, n
    integer :: spl, ispl
    real(r8) :: xx, yy, zz, rr, r2, alp, nrm, u, ux, uxx, u2, dx, dy, dz, tmp, &
            dx2, dy2, dz2, dxyz, dxdy, dxdz, dydz
    real(r8) :: df
    !cc
    real(r8), dimension(5) :: tmps !second index means:
    real(r8), dimension(0:2, 5) :: tmpp !1:   AO itsself
    real(r8), dimension(0:5, 5) :: tmpd !2-4: derivative with respect to x,y,z
    real(r8), dimension(0:9, 5) :: tmpf !5:   laplacian
    logical gaussFOrder       ! .t.: Gaussian order for f function
    !                               ! .f.: Gamess==Turbomole order used



    ! bf refers to the degenerate set of cartesian
    ! basis function (S:1,P:3,D:6,F:10) as input, which may be of type STO
    ! or contracted GTO.
    ! al refers to the individual basis function, as used in LCAO-MO's.
    ! (composed in subroutine mdetwf)
    ! i refers to the current electron.

    !-----Calculation of the AO's and their derivatives

    if (evfmt=='gau' .or. evfmt=='mol') then
        gaussFOrder = .true.
    else
        gaussFOrder = .false.
    endif

    if (ie == 0) then                     ! AO's for all electrons
        i1 = 1
        i2 = ne
    else
        i1 = ie                              ! only AO for electron ie
        i2 = ie
    endif

    do i = i1, i2                              ! loop over electrons
        xx = x(i)
        yy = y(i)
        zz = z(i)

        al = 1                               !Pointer for individual basis function
        moc = 0                              !Pointer for cmoa array

        ! initialisation

        mat(1:norb, i, 1) = 0d0
        mat1x(1:norb, i, 1) = 0d0
        mat1y(1:norb, i, 1) = 0d0
        mat1z(1:norb, i, 1) = 0d0
        mat2(1:norb, i, 1) = 0d0

        do bf = 1, nbasf                        ! loop over basis functions
            a = bc(bf)                         ! center of AO
            rr = rai(a, i)                      ! r_ai

            if (cutao) then                    !AO - Cutoff
                if (rr>aocuts(bf)) then       ! --> do nothing but adjusting the counters
                    if (bl(bf) == 'S') then
                        moc = moc + nmos(al)
                        al = al + 1
                    elseif (bl(bf) == 'P') then
                        do n = 0, 2
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 3
                    elseif (bl(bf) == 'D') then
                        do n = 0, 5
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 6
                    elseif (bl(bf) == 'F') then
                        do n = 0, 9
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 10
                    else
                        call abortp('(getaos): wrong GTO')
                    endif
                    cycle  !continue with next basis function
                endif
            endif

            r2 = rr * rr

            if (so(bf)==0) then !only 1 GTO in contraction, no splines used !

                if (bl(bf) == 'S') then                 ! 1s GTO

                    dx = xx - atoms(a)%cx
                    dy = yy - atoms(a)%cy
                    dz = zz - atoms(a)%cz

                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)
                    ux = -2d0 * alp * u

                    !cc tmps(1) not used here because u directly in MO loop
                    tmps(2) = ux * dx
                    tmps(3) = ux * dy
                    tmps(4) = ux * dz
                    tmps(5) = ux * (3d0 - 2d0 * alp * r2)

                    !cc MO calculation
                    do n = 1, nmos(al)
                        moc = moc + 1
                        tmp = cmoa(moc)
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + tmp * u
                        mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmps(2)
                        mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmps(3)
                        mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmps(4)
                        mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmps(5)
                    enddo

                    al = al + 1

                else if (bl(bf) == 'P') then             ! 2p GTO's
                    !              // do all 3 P simultaneously (same exponent is required)
                    !              // order p_x,p_y,p_z

                    dx = xx - atoms(a)%cx
                    dy = yy - atoms(a)%cy
                    dz = zz - atoms(a)%cz

                    dxdy = dx * dy
                    dxdz = dx * dz
                    dydz = dy * dz

                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)
                    ux = -2d0 * alp * u

                    tmpp(0, 1) = dx * u
                    tmpp(1, 1) = dy * u
                    tmpp(2, 1) = dz * u

                    tmpp(0, 2) = u + ux * dx * dx
                    tmpp(1, 2) = ux * dxdy
                    tmpp(2, 2) = ux * dxdz
                    tmpp(0, 3) = ux * dxdy
                    tmpp(1, 3) = u + ux * dy * dy
                    tmpp(2, 3) = ux * dydz
                    tmpp(0, 4) = ux * dxdz
                    tmpp(1, 4) = ux * dydz
                    tmpp(2, 4) = u + ux * dz * dz

                    tmp = (5d0 - 2d0 * alp * r2) * ux
                    tmpp(0, 5) = tmp * dx
                    tmpp(1, 5) = tmp * dy
                    tmpp(2, 5) = tmp * dz

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)

                    do d = 0, 2
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            tmp = cmoa(moc)
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + tmp * tmpp(d, 1)
                            mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpp(d, 2)
                            mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpp(d, 3)
                            mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpp(d, 4)
                            mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpp(d, 5)
                        enddo
                    enddo

                    al = al + 3

                else if (bl(bf) == 'D') then         ! 3d GTO
                    !              // do all 6 D simultaneously (same exponent is required)
                    !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)

                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz

                    dxdy = dx * dy
                    dxdz = dx * dz
                    dydz = dy * dz
                    !cc
                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)
                    ux = -2d0 * alp * u

                    tmpd(0, 1) = dx2 * u
                    tmpd(1, 1) = dy2 * u
                    tmpd(2, 1) = dz2 * u

                    tmpd(0, 2) = (2d0 * u + ux * dx2) * dx
                    tmpd(1, 2) = dxdy * dy * ux
                    tmpd(2, 2) = dxdz * dz * ux
                    tmpd(0, 3) = dxdy * dx * ux
                    tmpd(1, 3) = (2d0 * u + ux * dy2) * dy
                    tmpd(2, 3) = dydz * dz * ux
                    tmpd(0, 4) = dxdz * dx * ux
                    tmpd(1, 4) = dydz * dy * ux
                    tmpd(2, 4) = (2d0 * u + ux * dz2) * dz
                    tmp = (7d0 - 2d0 * alp * r2) * ux
                    tmpd(0, 5) = 2d0 * u + dx2 * tmp
                    tmpd(1, 5) = 2d0 * u + dy2 * tmp
                    tmpd(2, 5) = 2d0 * u + dz2 * tmp

                    u = sqr3 * u                   ! correction of norm
                    ux = sqr3 * ux                 ! N(dxx)*sqr3 = N(dxy)

                    tmpd(3, 1) = dxdy * u
                    tmpd(4, 1) = dxdz * u
                    tmpd(5, 1) = dydz * u
                    tmp = ux * dx * dy * dz
                    tmpd(3, 2) = (u + ux * dx2) * dy
                    tmpd(4, 2) = (u + ux * dx2) * dz
                    tmpd(5, 2) = tmp
                    tmpd(3, 3) = (u + ux * dy2) * dx
                    tmpd(4, 3) = tmp
                    tmpd(5, 3) = (u + ux * dy2) * dz
                    tmpd(3, 4) = tmp
                    tmpd(4, 4) = (u + ux * dz2) * dx
                    tmpd(5, 4) = (u + ux * dz2) * dy
                    tmp = (7d0 - 2d0 * alp * r2) * ux
                    tmpd(3, 5) = tmp * dxdy
                    tmpd(4, 5) = tmp * dxdz
                    tmpd(5, 5) = tmp * dydz

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)

                    do d = 0, 5
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            tmp = cmoa(moc)
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + tmp * tmpd(d, 1)
                            mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpd(d, 2)
                            mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpd(d, 3)
                            mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpd(d, 4)
                            mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpd(d, 5)
                        enddo
                    enddo

                    al = al + 6

                else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! f GTO
                    !              // do all 10 F simultaneously (same exponent is required)
                    !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                    !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)
                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz
                    dxyz = dx * dy * dz
                    ux = -2d0 * alp * u

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0, 1) = dx2 * dx * u
                    tmpf(1, 1) = dy2 * dy * u
                    tmpf(2, 1) = dz2 * dz * u

                    tmpf(0, 2) = (3d0 * u + ux * dx2) * dx2
                    tmpf(1, 2) = dy2 * dy * ux * dx
                    tmpf(2, 2) = dz2 * dz * ux * dx
                    tmpf(0, 3) = dx2 * dx * ux * dy
                    tmpf(1, 3) = (3d0 * u + ux * dy2) * dy2
                    tmpf(2, 3) = dz2 * dz * ux * dy
                    tmpf(0, 4) = dx2 * dx * ux * dz
                    tmpf(1, 4) = dy2 * dy * ux * dz
                    tmpf(2, 4) = (3d0 * u + ux * dz2) * dz2
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(0, 5) = (6d0 * u + dx2 * tmp) * dx
                    tmpf(1, 5) = (6d0 * u + dy2 * tmp) * dy
                    tmpf(2, 5) = (6d0 * u + dz2 * tmp) * dz

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    ux = sqr5 * ux                 ! N(fxxx)*sqrt(5) = N(fxxy)

                    tmpf(3, 1) = dx2 * dy * u
                    tmpf(4, 1) = dx2 * dz * u
                    tmpf(5, 1) = dy2 * dx * u
                    tmpf(6, 1) = dy2 * dz * u
                    tmpf(7, 1) = dz2 * dx * u
                    tmpf(8, 1) = dz2 * dy * u

                    ! derivatives
                    tmp = ux * dxyz
                    tmpf(3, 2) = (2d0 * u + ux * dx2) * dx * dy
                    tmpf(4, 2) = (2d0 * u + ux * dx2) * dx * dz
                    tmpf(5, 2) = (u + ux * dx2) * dy2
                    tmpf(6, 2) = tmp * dy
                    tmpf(7, 2) = (u + ux * dx2) * dz2
                    tmpf(8, 2) = tmp * dz
                    tmpf(3, 3) = (u + ux * dy2) * dx2
                    tmpf(4, 3) = tmp * dx
                    tmpf(5, 3) = (2d0 * u + ux * dy2) * dx * dy
                    tmpf(6, 3) = (2d0 * u + ux * dy2) * dy * dz
                    tmpf(7, 3) = tmp * dz
                    tmpf(8, 3) = (u + ux * dy2) * dz2
                    tmpf(3, 4) = tmp * dx
                    tmpf(4, 4) = (u + ux * dz2) * dx2
                    tmpf(5, 4) = tmp * dy
                    tmpf(6, 4) = (u + ux * dz2) * dy2
                    tmpf(7, 4) = (2d0 * u + ux * dz2) * dx * dz
                    tmpf(8, 4) = (2d0 * u + ux * dz2) * dy * dz
                    ! laplacians
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(3, 5) = (2d0 * u + dx2 * tmp) * dy
                    tmpf(4, 5) = (2d0 * u + dx2 * tmp) * dz
                    tmpf(5, 5) = (2d0 * u + dy2 * tmp) * dx
                    tmpf(6, 5) = (2d0 * u + dy2 * tmp) * dz
                    tmpf(7, 5) = (2d0 * u + dz2 * tmp) * dx
                    tmpf(8, 5) = (2d0 * u + dz2 * tmp) * dy

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    ux = sqr3 * ux                ! N(fxxx)*sqrt(15)=
                    !                                          ! N(fxxy)*sqrt(3)=N(fxyz)
                    tmpf(9, 1) = dxyz * u

                    tmpf(9, 2) = (u + ux * dx2) * dy * dz
                    tmpf(9, 3) = (u + ux * dy2) * dx * dz
                    tmpf(9, 4) = (u + ux * dz2) * dx * dy
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(9, 5) = dxyz * tmp


                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)

                    do d = 0, 9
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            tmp = cmoa(moc)
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + tmp * tmpf(d, 1)
                            mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpf(d, 2)
                            mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpf(d, 3)
                            mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpf(d, 4)
                            mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpf(d, 5)
                        enddo
                    enddo

                    al = al + 10

                else if (bl(bf)=='F'.and.gaussFOrder) then     ! f GTO
                    !              // do all 10 F simultaneously (same exponent is required)
                    !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                    !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)
                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz
                    dxyz = dx * dy * dz
                    ux = -2d0 * alp * u

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0, 1) = dx2 * dx * u
                    tmpf(1, 1) = dy2 * dy * u
                    tmpf(2, 1) = dz2 * dz * u

                    tmpf(0, 2) = (3d0 * u + ux * dx2) * dx2
                    tmpf(1, 2) = dy2 * dy * ux * dx
                    tmpf(2, 2) = dz2 * dz * ux * dx
                    tmpf(0, 3) = dx2 * dx * ux * dy
                    tmpf(1, 3) = (3d0 * u + ux * dy2) * dy2
                    tmpf(2, 3) = dz2 * dz * ux * dy
                    tmpf(0, 4) = dx2 * dx * ux * dz
                    tmpf(1, 4) = dy2 * dy * ux * dz
                    tmpf(2, 4) = (3d0 * u + ux * dz2) * dz2
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(0, 5) = (6d0 * u + dx2 * tmp) * dx
                    tmpf(1, 5) = (6d0 * u + dy2 * tmp) * dy
                    tmpf(2, 5) = (6d0 * u + dz2 * tmp) * dz

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    ux = sqr5 * ux                 ! N(fxxx)*sqrt(5) = N(fxxy)

                    tmpf(4, 1) = dx2 * dy * u
                    tmpf(5, 1) = dx2 * dz * u
                    tmpf(3, 1) = dy2 * dx * u
                    tmpf(8, 1) = dy2 * dz * u
                    tmpf(6, 1) = dz2 * dx * u
                    tmpf(7, 1) = dz2 * dy * u

                    ! derivatives
                    tmp = ux * dxyz
                    tmpf(4, 2) = (2d0 * u + ux * dx2) * dx * dy
                    tmpf(5, 2) = (2d0 * u + ux * dx2) * dx * dz
                    tmpf(3, 2) = (u + ux * dx2) * dy2
                    tmpf(8, 2) = tmp * dy
                    tmpf(6, 2) = (u + ux * dx2) * dz2
                    tmpf(7, 2) = tmp * dz
                    tmpf(4, 3) = (u + ux * dy2) * dx2
                    tmpf(5, 3) = tmp * dx
                    tmpf(3, 3) = (2d0 * u + ux * dy2) * dx * dy
                    tmpf(8, 3) = (2d0 * u + ux * dy2) * dy * dz
                    tmpf(6, 3) = tmp * dz
                    tmpf(7, 3) = (u + ux * dy2) * dz2
                    tmpf(4, 4) = tmp * dx
                    tmpf(5, 4) = (u + ux * dz2) * dx2
                    tmpf(3, 4) = tmp * dy
                    tmpf(8, 4) = (u + ux * dz2) * dy2
                    tmpf(6, 4) = (2d0 * u + ux * dz2) * dx * dz
                    tmpf(7, 4) = (2d0 * u + ux * dz2) * dy * dz
                    ! laplacians
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(4, 5) = (2d0 * u + dx2 * tmp) * dy
                    tmpf(5, 5) = (2d0 * u + dx2 * tmp) * dz
                    tmpf(3, 5) = (2d0 * u + dy2 * tmp) * dx
                    tmpf(8, 5) = (2d0 * u + dy2 * tmp) * dz
                    tmpf(6, 5) = (2d0 * u + dz2 * tmp) * dx
                    tmpf(7, 5) = (2d0 * u + dz2 * tmp) * dy

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    ux = sqr3 * ux                ! N(fxxx)*sqrt(15)=
                    !                                          ! N(fxxy)*sqrt(3)=N(fxyz)
                    tmpf(9, 1) = dxyz * u

                    tmpf(9, 2) = (u + ux * dx2) * dy * dz
                    tmpf(9, 3) = (u + ux * dy2) * dx * dz
                    tmpf(9, 4) = (u + ux * dz2) * dx * dy
                    tmp = (9d0 - 2d0 * alp * r2) * ux
                    tmpf(9, 5) = dxyz * tmp


                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)

                    do d = 0, 9
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            tmp = cmoa(moc)
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + tmp * tmpf(d, 1)
                            mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpf(d, 2)
                            mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpf(d, 3)
                            mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpf(d, 4)
                            mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpf(d, 5)
                        enddo
                    enddo

                    al = al + 10

                else
                    call abortp('(getaos): wrong GTO')
                endif  ! bl

            else ! CGTO (more than one primitive gaussian) --> use splines

                !cc            r2 = rr*rr
                spl = (csplnpnt - 1) * rr / (csalpha + rr) + 1
                df = rr - csplx(spl)

                if (bl(bf) == 'S') then                 ! 1s GTO

                    ispl = 3 * so(bf) - 2
                    tmps(1) = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    ux = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    u2 = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dy = yy - atoms(a)%cy
                    dz = zz - atoms(a)%cz

                    tmps(2) = ux * dx / rr
                    tmps(3) = ux * dy / rr
                    tmps(4) = ux * dz / rr
                    tmps(5) = u2 + 2 * ux / rr

                    !cc MO calculation
                    do n = 1, nmos(al)
                        moc = moc + 1
                        tmp = cmoa(moc)
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + tmp * tmps(1)
                        mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmps(2)
                        mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmps(3)
                        mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmps(4)
                        mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmps(5)
                    enddo

                    al = al + 1

                else if (bl(bf) == 'P') then             ! 2p GTO's
                    !              // do all 3 P simultaneously (same exponent is required)
                    !              // order p_x,p_y,p_z

                    ispl = 3 * so(bf) - 2
                    u = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    ux = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    uxx = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dy = yy - atoms(a)%cy
                    dz = zz - atoms(a)%cz

                    dxdy = dx * dy
                    dxdz = dx * dz
                    dydz = dy * dz

                    !cc
                    tmpp(0, 1) = dx * u
                    tmpp(1, 1) = dy * u
                    tmpp(2, 1) = dz * u

                    tmpp(0, 2) = u + ux * dx * dx
                    tmpp(1, 2) = ux * dxdy
                    tmpp(2, 2) = ux * dxdz
                    tmpp(0, 3) = ux * dxdy
                    tmpp(1, 3) = u + ux * dy * dy
                    tmpp(2, 3) = ux * dydz
                    tmpp(0, 4) = ux * dxdz
                    tmpp(1, 4) = ux * dydz
                    tmpp(2, 4) = u + ux * dz * dz

                    tmpp(0, 5) = uxx * dx
                    tmpp(1, 5) = uxx * dy
                    tmpp(2, 5) = uxx * dz

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 2
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            tmp = cmoa(moc)
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + tmp * tmpp(d, 1)
                            mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpp(d, 2)
                            mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpp(d, 3)
                            mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpp(d, 4)
                            mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpp(d, 5)
                        enddo
                    enddo

                    al = al + 3

                else if (bl(bf) == 'D') then         ! 3d GTO
                    !              // do all 6 D simultaneously (same exponent is required)
                    !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)

                    ispl = 3 * so(bf) - 2
                    u = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    ux = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    uxx = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz
                    dxdy = dx * dy
                    dxdz = dx * dz
                    dydz = dy * dz

                    tmpd(0, 1) = dx2 * u
                    tmpd(1, 1) = dy2 * u
                    tmpd(2, 1) = dz2 * u

                    tmpd(0, 2) = (2d0 * u + ux * dx2) * dx
                    tmpd(1, 2) = dy2 * ux * dx
                    tmpd(2, 2) = dz2 * ux * dx
                    tmpd(0, 3) = dx2 * ux * dy
                    tmpd(1, 3) = (2d0 * u + ux * dy2) * dy
                    tmpd(2, 3) = dz2 * ux * dy
                    tmpd(0, 4) = dx2 * ux * dz
                    tmpd(1, 4) = dy2 * ux * dz
                    tmpd(2, 4) = (2d0 * u + ux * dz2) * dz

                    tmpd(0, 5) = 2d0 * u + dx2 * uxx
                    tmpd(1, 5) = 2d0 * u + dy2 * uxx
                    tmpd(2, 5) = 2d0 * u + dz2 * uxx

                    u = sqr3 * u                   ! correction of norm
                    ux = sqr3 * ux                 ! N(dxx)*sqr3 = N(dxy)
                    uxx = sqr3 * uxx

                    tmpd(3, 1) = dxdy * u
                    tmpd(4, 1) = dxdz * u
                    tmpd(5, 1) = dydz * u
                    tmp = ux * dx * dy * dz
                    tmpd(3, 2) = (u + ux * dx2) * dy
                    tmpd(4, 2) = (u + ux * dx2) * dz
                    tmpd(5, 2) = tmp
                    tmpd(3, 3) = (u + ux * dy2) * dx
                    tmpd(4, 3) = tmp
                    tmpd(5, 3) = (u + ux * dy2) * dz
                    tmpd(3, 4) = tmp
                    tmpd(4, 4) = (u + ux * dz2) * dx
                    tmpd(5, 4) = (u + ux * dz2) * dy

                    tmpd(3, 5) = uxx * dxdy
                    tmpd(4, 5) = uxx * dxdz
                    tmpd(5, 5) = uxx * dydz

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 5
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            tmp = cmoa(moc)
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + tmp * tmpd(d, 1)
                            mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpd(d, 2)
                            mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpd(d, 3)
                            mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpd(d, 4)
                            mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpd(d, 5)
                        enddo
                    enddo

                    al = al + 6

                else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! f GTO
                    !              // do all 10 F simultaneously (same exponent is required)
                    !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                    !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

                    ispl = 3 * so(bf) - 2
                    u = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    ux = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    uxx = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz
                    dxyz = dx * dy * dz

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0, 1) = dx2 * dx * u
                    tmpf(1, 1) = dy2 * dy * u
                    tmpf(2, 1) = dz2 * dz * u

                    tmpf(0, 2) = (3d0 * u + ux * dx2) * dx2
                    tmpf(1, 2) = dy2 * dy * ux * dx
                    tmpf(2, 2) = dz2 * dz * ux * dx
                    tmpf(0, 3) = dx2 * dx * ux * dy
                    tmpf(1, 3) = (3d0 * u + ux * dy2) * dy2
                    tmpf(2, 3) = dz2 * dz * ux * dy
                    tmpf(0, 4) = dx2 * dx * ux * dz
                    tmpf(1, 4) = dy2 * dy * ux * dz
                    tmpf(2, 4) = (3d0 * u + ux * dz2) * dz2

                    tmpf(0, 5) = (6d0 * u + dx2 * uxx) * dx
                    tmpf(1, 5) = (6d0 * u + dy2 * uxx) * dy
                    tmpf(2, 5) = (6d0 * u + dz2 * uxx) * dz

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    ux = sqr5 * ux                 ! N(fxxx)*sqrt(5) = N(fxxy)
                    uxx = sqr5 * uxx

                    tmpf(3, 1) = dx2 * dy * u
                    tmpf(4, 1) = dx2 * dz * u
                    tmpf(5, 1) = dy2 * dx * u
                    tmpf(6, 1) = dy2 * dz * u
                    tmpf(7, 1) = dz2 * dx * u
                    tmpf(8, 1) = dz2 * dy * u

                    ! derivatives
                    tmp = ux * dxyz
                    tmpf(3, 2) = (2d0 * u + ux * dx2) * dx * dy
                    tmpf(4, 2) = (2d0 * u + ux * dx2) * dx * dz
                    tmpf(5, 2) = (u + ux * dx2) * dy2
                    tmpf(6, 2) = tmp * dy
                    tmpf(7, 2) = (u + ux * dx2) * dz2
                    tmpf(8, 2) = tmp * dz
                    tmpf(3, 3) = (u + ux * dy2) * dx2
                    tmpf(4, 3) = tmp * dx
                    tmpf(5, 3) = (2d0 * u + ux * dy2) * dx * dy
                    tmpf(6, 3) = (2d0 * u + ux * dy2) * dy * dz
                    tmpf(7, 3) = tmp * dz
                    tmpf(8, 3) = (u + ux * dy2) * dz2
                    tmpf(3, 4) = tmp * dx
                    tmpf(4, 4) = (u + ux * dz2) * dx2
                    tmpf(5, 4) = tmp * dy
                    tmpf(6, 4) = (u + ux * dz2) * dy2
                    tmpf(7, 4) = (2d0 * u + ux * dz2) * dx * dz
                    tmpf(8, 4) = (2d0 * u + ux * dz2) * dy * dz
                    ! laplacians
                    tmpf(3, 5) = (2d0 * u + dx2 * uxx) * dy
                    tmpf(4, 5) = (2d0 * u + dx2 * uxx) * dz
                    tmpf(5, 5) = (2d0 * u + dy2 * uxx) * dx
                    tmpf(6, 5) = (2d0 * u + dy2 * uxx) * dz
                    tmpf(7, 5) = (2d0 * u + dz2 * uxx) * dx
                    tmpf(8, 5) = (2d0 * u + dz2 * uxx) * dy

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    ux = sqr3 * ux                ! N(fxxx)*sqrt(15)=
                    uxx = sqr3 * uxx              ! N(fxxy)*sqrt(3)=N(fxyz)

                    tmpf(9, 1) = dxyz * u

                    tmpf(9, 2) = (u + ux * dx2) * dy * dz
                    tmpf(9, 3) = (u + ux * dy2) * dx * dz
                    tmpf(9, 4) = (u + ux * dz2) * dx * dy
                    tmpf(9, 5) = dxyz * uxx

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 9
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            tmp = cmoa(moc)
                            j = mo_o(j)
                            mat(j, i, 1) = mat(j, i, 1) + tmp * tmpf(d, 1)
                            mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpf(d, 2)
                            mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpf(d, 3)
                            mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpf(d, 4)
                            mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpf(d, 5)
                        enddo
                    enddo

                    al = al + 10

                else if (bl(bf)=='F'.and.gaussFOrder) then     ! f GTO
                    !              // do all 10 F simultaneously (same exponent is required)
                    !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                    !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

                    ispl = 3 * so(bf) - 2
                    u = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    ux = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))
                    ispl = ispl + 1
                    uxx = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz
                    dxyz = dx * dy * dz

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0, 1) = dx2 * dx * u
                    tmpf(1, 1) = dy2 * dy * u
                    tmpf(2, 1) = dz2 * dz * u

                    tmpf(0, 2) = (3d0 * u + ux * dx2) * dx2
                    tmpf(1, 2) = dy2 * dy * ux * dx
                    tmpf(2, 2) = dz2 * dz * ux * dx
                    tmpf(0, 3) = dx2 * dx * ux * dy
                    tmpf(1, 3) = (3d0 * u + ux * dy2) * dy2
                    tmpf(2, 3) = dz2 * dz * ux * dy
                    tmpf(0, 4) = dx2 * dx * ux * dz
                    tmpf(1, 4) = dy2 * dy * ux * dz
                    tmpf(2, 4) = (3d0 * u + ux * dz2) * dz2

                    tmpf(0, 5) = (6d0 * u + dx2 * uxx) * dx
                    tmpf(1, 5) = (6d0 * u + dy2 * uxx) * dy
                    tmpf(2, 5) = (6d0 * u + dz2 * uxx) * dz

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    ux = sqr5 * ux                 ! N(fxxx)*sqrt(5) = N(fxxy)
                    uxx = sqr5 * uxx

                    tmpf(3, 1) = dx2 * dy * u
                    tmpf(4, 1) = dx2 * dz * u
                    tmpf(5, 1) = dy2 * dx * u
                    tmpf(6, 1) = dy2 * dz * u
                    tmpf(7, 1) = dz2 * dx * u
                    tmpf(8, 1) = dz2 * dy * u

                    ! derivatives
                    tmp = ux * dxyz
                    tmpf(4, 2) = (2d0 * u + ux * dx2) * dx * dy
                    tmpf(5, 2) = (2d0 * u + ux * dx2) * dx * dz
                    tmpf(3, 2) = (u + ux * dx2) * dy2
                    tmpf(8, 2) = tmp * dy
                    tmpf(6, 2) = (u + ux * dx2) * dz2
                    tmpf(7, 2) = tmp * dz
                    tmpf(4, 3) = (u + ux * dy2) * dx2
                    tmpf(5, 3) = tmp * dx
                    tmpf(3, 3) = (2d0 * u + ux * dy2) * dx * dy
                    tmpf(8, 3) = (2d0 * u + ux * dy2) * dy * dz
                    tmpf(6, 3) = tmp * dz
                    tmpf(7, 3) = (u + ux * dy2) * dz2
                    tmpf(4, 4) = tmp * dx
                    tmpf(5, 4) = (u + ux * dz2) * dx2
                    tmpf(3, 4) = tmp * dy
                    tmpf(8, 4) = (u + ux * dz2) * dy2
                    tmpf(6, 4) = (2d0 * u + ux * dz2) * dx * dz
                    tmpf(7, 4) = (2d0 * u + ux * dz2) * dy * dz
                    ! laplacians
                    tmpf(4, 5) = (2d0 * u + dx2 * uxx) * dy
                    tmpf(5, 5) = (2d0 * u + dx2 * uxx) * dz
                    tmpf(3, 5) = (2d0 * u + dy2 * uxx) * dx
                    tmpf(8, 5) = (2d0 * u + dy2 * uxx) * dz
                    tmpf(6, 5) = (2d0 * u + dz2 * uxx) * dx
                    tmpf(7, 5) = (2d0 * u + dz2 * uxx) * dy

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    ux = sqr3 * ux                ! N(fxxx)*sqrt(15)=
                    uxx = sqr3 * uxx              ! N(fxxy)*sqrt(3)=N(fxyz)

                    tmpf(9, 1) = dxyz * u

                    tmpf(9, 2) = (u + ux * dx2) * dy * dz
                    tmpf(9, 3) = (u + ux * dy2) * dx * dz
                    tmpf(9, 4) = (u + ux * dz2) * dx * dy
                    tmpf(9, 5) = dxyz * uxx

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 9
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            tmp = cmoa(moc)
                            j = mo_o(j)
                            mat(j, i, 1) = mat(j, i, 1) + tmp * tmpf(d, 1)
                            mat1x(j, i, 1) = mat1x(j, i, 1) + tmp * tmpf(d, 2)
                            mat1y(j, i, 1) = mat1y(j, i, 1) + tmp * tmpf(d, 3)
                            mat1z(j, i, 1) = mat1z(j, i, 1) + tmp * tmpf(d, 4)
                            mat2(j, i, 1) = mat2(j, i, 1) + tmp * tmpf(d, 5)
                        enddo
                    enddo

                    al = al + 10

                else
                    call abortp('(getaos): wrong GTO')
                endif ! bl
            endif  ! CGTO or primitive gaussian function
        enddo    ! bf-loop over basis functions
    enddo       ! i-loop over electrons

end subroutine aomocutspl_calc


!     ----------------------------------

subroutine aomocut1spl_calc(ie, x, y, z, rai)

    use wfData_m
    use aoMo_m
    use cspline_m

    implicit none

    ! input parameters:
    integer ie                 ! if >0 only AO's for electron ie recalculated

    real(r8), dimension(nmax) :: x, y, z    ! x,y,z coordinates of position vector
    real(r8), dimension(amax, nmax) :: rai ! r_ai electron-nucleus distances

    ! constants:
    real(r8) :: sqr3, sqr5
    parameter (sqr3 = 1.73205080756887729d0, sqr5 = 2.236067977499789696d0)
    ! local variables
    integer :: bf, a, i, i1, i2, ii, ic, al
    integer :: j, d, moc, n
    integer :: spl, ispl
    real(r8) :: xx, yy, zz, rr, r2, alp, nrm, u, ux, uxx, u2, dx, dy, dz, tmp, &
            dx2, dy2, dz2, dxdy, dxdz, dydz
    real(r8) :: df
    !cc
    real(r8) :: tmps
    real(r8), dimension(0:2) :: tmpp
    real(r8), dimension(0:5) :: tmpd
    real(r8), dimension(0:9) :: tmpf
    logical gaussFOrder       ! .t.: Gaussian order for f function
    !                               ! .f.: Gamess==Turbomole order used



    ! bf refers to the degenerate set of cartesian
    ! basis function (S:1,P:3,D:6,F:10) as input, which may be of type STO
    ! or contracted GTO.
    ! al refers to the individual basis function, as used in LCAO-MO's.
    ! (composed in subroutine mdetwf)
    ! i refers to the current electron.

    !-----Calculation of the AO's and their derivatives

    if (evfmt=='gau' .or. evfmt=='mol') then
        gaussFOrder = .true.
    else
        gaussFOrder = .false.
    endif

    if (ie == 0) then                     ! AO's for all electrons
        i1 = 1
        i2 = ne
    else
        i1 = ie                              ! only AO for electron ie
        i2 = ie
    endif

    do i = i1, i2                              ! loop over electrons
        xx = x(i)
        yy = y(i)
        zz = z(i)

        al = 1                               !Pointer for indivudial basis function
        moc = 0                              !Pointer for cmoa array

        ! initialisation

        mat(1:norb, i, 1) = 0d0

        do bf = 1, nbasf                        ! loop over basis functions
            a = bc(bf)                         ! center of AO
            rr = rai(a, i)                      ! r_ai

            if (cutao) then                    !AO - Cutoff
                if (rr>aocuts(bf)) then       ! --> do nothing but adjusting the counters
                    if (bl(bf) == 'S') then
                        moc = moc + nmos(al)
                        al = al + 1
                    elseif (bl(bf) == 'P') then
                        do n = 0, 2
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 3
                    elseif (bl(bf) == 'D') then
                        do n = 0, 5
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 6
                    elseif (bl(bf) == 'F') then
                        do n = 0, 9
                            moc = moc + nmos(al + n)
                        enddo
                        al = al + 10
                    else
                        call abortp('(getaos): wrong GTO')
                    endif
                    cycle  !continue with next basis function
                endif
            endif

            r2 = rr * rr

            if (so(bf)==0) then !only 1 GTO in contraction, no splines used !

                if (bl(bf) == 'S') then                 ! 1s GTO

                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)

                    !cc MO calculation
                    do n = 1, nmos(al)
                        moc = moc + 1
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * u
                    enddo

                    al = al + 1

                else if (bl(bf) == 'P') then             ! 2p GTO's
                    !              // do all 3 P simultaneously (same exponent is required)
                    !              // order p_x,p_y,p_z

                    dx = xx - atoms(a)%cx
                    dy = yy - atoms(a)%cy
                    dz = zz - atoms(a)%cz

                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)

                    tmpp(0) = dx * u
                    tmpp(1) = dy * u
                    tmpp(2) = dz * u

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 2
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpp(d)
                        enddo
                    enddo

                    al = al + 3

                else if (bl(bf) == 'D') then         ! 3d GTO
                    !              // do all 6 D simultaneously (same exponent is required)
                    !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)

                    dx = xx - atoms(a)%cx
                    dy = yy - atoms(a)%cy
                    dz = zz - atoms(a)%cz
                    !cc
                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)

                    tmpd(0) = dx * dx * u
                    tmpd(1) = dy * dy * u
                    tmpd(2) = dz * dz * u

                    u = sqr3 * u                   ! correction of norm
                    !                                           ! N(dxx)*sqr3 = N(dxy)

                    tmpd(3) = dx * dy * u
                    tmpd(4) = dx * dz * u
                    tmpd(5) = dy * dz * u

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 5
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpd(d)
                        enddo
                    enddo

                    al = al + 6

                else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! f GTO
                    !              // do all 10 F simultaneously (same exponent is required)
                    !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                    !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)

                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0) = dx2 * dx * u
                    tmpf(1) = dy2 * dy * u
                    tmpf(2) = dz2 * dz * u

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    !                                           ! N(fxxx)*sqrt(5) = N(fxxy)

                    tmpf(3) = dx2 * dy * u
                    tmpf(4) = dx2 * dz * u
                    tmpf(5) = dy2 * dx * u
                    tmpf(6) = dy2 * dz * u
                    tmpf(7) = dz2 * dx * u
                    tmpf(8) = dz2 * dy * u

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    !                                          ! N(fxxx)*sqrt(15)=
                    !                                          ! N(fxxy)*sqrt(3)=N(fxyz)
                    tmpf(9) = dx * dy * dz * u

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 9
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpf(d)
                        enddo
                    enddo

                    al = al + 10

                else if (bl(bf)=='F'.and.gaussFOrder) then     ! f GTO
                    !              // do all 10 F simultaneously (same exponent is required)
                    !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                    !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

                    alp = cntrctn(1, 1, bf)
                    u = cntrctn(2, 1, bf) * exp(-alp * r2)

                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0) = dx2 * dx * u
                    tmpf(1) = dy2 * dy * u
                    tmpf(2) = dz2 * dz * u

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    !                                           ! N(fxxx)*sqrt(5) = N(fxxy)

                    tmpf(3) = dy2 * dx * u
                    tmpf(4) = dx2 * dy * u
                    tmpf(5) = dx2 * dz * u
                    tmpf(6) = dz2 * dx * u
                    tmpf(7) = dz2 * dy * u
                    tmpf(8) = dy2 * dz * u

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    !                                          ! N(fxxx)*sqrt(15)=
                    !                                          ! N(fxxy)*sqrt(3)=N(fxyz)
                    tmpf(9) = dx * dy * dz * u

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 9
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpf(d)
                        enddo
                    enddo

                    al = al + 10

                else
                    call abortp('(getaos): wrong GTO')
                endif  ! bl

            else ! CGTO (more than one primitive gaussian) --> use splines

                !cc            r2 = rr*rr
                spl = (csplnpnt - 1) * rr / (csalpha + rr) + 1
                df = rr - csplx(spl)

                if (bl(bf) == 'S') then                 ! 1s GTO

                    ispl = 3 * so(bf) - 2
                    tmps = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    !cc MO calculation
                    do n = 1, nmos(al)
                        moc = moc + 1
                        j = mo_o(moc)
                        mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmps
                    enddo

                    al = al + 1

                else if (bl(bf) == 'P') then             ! 2p GTO's
                    !              // do all 3 P simultaneously (same exponent is required)
                    !              // order p_x,p_y,p_z

                    ispl = 3 * so(bf) - 2
                    u = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dy = yy - atoms(a)%cy
                    dz = zz - atoms(a)%cz

                    tmpp(0) = dx * u
                    tmpp(1) = dy * u
                    tmpp(2) = dz * u

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 2
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpp(d)
                        enddo
                    enddo

                    al = al + 3

                else if (bl(bf) == 'D') then         ! 3d GTO
                    !              // do all 6 D simultaneously (same exponent is required)
                    !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)

                    ispl = 3 * so(bf) - 2
                    u = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dy = yy - atoms(a)%cy
                    dz = zz - atoms(a)%cz

                    tmpd(0) = dx * dx * u
                    tmpd(1) = dy * dy * u
                    tmpd(2) = dz * dz * u

                    u = sqr3 * u                   ! correction of norm
                    !                                           ! N(dxx)*sqr3 = N(dxy)

                    tmpd(3) = dx * dy * u
                    tmpd(4) = dx * dz * u
                    tmpd(5) = dy * dz * u

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 5
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpd(d)
                        enddo
                    enddo

                    al = al + 6

                else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! f GTO
                    !              // do all 10 F simultaneously (same exponent is required)
                    !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                    !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

                    ispl = 3 * so(bf) - 2
                    u = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0) = dx2 * dx * u
                    tmpf(1) = dy2 * dy * u
                    tmpf(2) = dz2 * dz * u

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    !                                           ! N(fxxx)*sqrt(5) = N(fxxy)

                    tmpf(3) = dx2 * dy * u
                    tmpf(4) = dx2 * dz * u
                    tmpf(5) = dy2 * dx * u
                    tmpf(6) = dy2 * dz * u
                    tmpf(7) = dz2 * dx * u
                    tmpf(8) = dz2 * dy * u

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    !                                          ! N(fxxx)*sqrt(15)=
                    !                                          ! N(fxxy)*sqrt(3)=N(fxyz)

                    tmpf(9) = dx * dy * dz * u

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 9
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpf(d)
                        enddo
                    enddo

                    al = al + 10

                else if (bl(bf)=='F'.and.gaussFOrder) then     ! f GTO
                    !              // do all 10 F simultaneously (same exponent is required)
                    !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                    !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

                    ispl = 3 * so(bf) - 2
                    u = cspla(ispl, spl) + df * (csplb(ispl, spl)&
                            + df * (csplc(ispl, spl) + df * cspld(ispl, spl)))

                    dx = xx - atoms(a)%cx
                    dx2 = dx * dx
                    dy = yy - atoms(a)%cy
                    dy2 = dy * dy
                    dz = zz - atoms(a)%cz
                    dz2 = dz * dz

                    !                 // f_xxx, f_yyy, f_zzz
                    tmpf(0) = dx2 * dx * u
                    tmpf(1) = dy2 * dy * u
                    tmpf(2) = dz2 * dz * u

                    !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                    u = sqr5 * u                   ! correction of norm
                    !                                           ! N(fxxx)*sqrt(5) = N(fxxy)

                    tmpf(3) = dy2 * dx * u
                    tmpf(4) = dx2 * dy * u
                    tmpf(5) = dx2 * dz * u
                    tmpf(6) = dz2 * dx * u
                    tmpf(7) = dz2 * dy * u
                    tmpf(8) = dy2 * dz * u

                    !                 // f_xyz
                    u = sqr3 * u                  ! correction of norm
                    !                                          ! N(fxxx)*sqrt(15)=
                    !                                          ! N(fxxy)*sqrt(3)=N(fxyz)

                    tmpf(9) = dx * dy * dz * u

                    !cc MO calculation (It's up to the compiler to unroll this loop if that's faster)
                    do d = 0, 9
                        do n = 1, nmos(al + d)
                            moc = moc + 1
                            j = mo_o(moc)
                            mat(j, i, 1) = mat(j, i, 1) + cmoa(moc) * tmpf(d)
                        enddo
                    enddo

                    al = al + 10

                else
                    call abortp('(getaos): wrong GTO')
                endif ! bl
            endif  ! CGTO or primitive gaussian function
        enddo    ! bf-loop over basis functions
    enddo       ! i-loop over electrons

end subroutine aomocut1spl_calc