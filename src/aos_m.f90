! Copyright (C) 1998-1999, 2005, 2007-2008 2012-2014 Arne Luechow
! Copyright (C) 1999, 2000 Sebastian Manten
! Copyright (C) 2003-2004 Christian Diedrich
! Copyright (C) 2004-2005 Annika Bande
! Copyright (C) 2005 Tony C. Scott
! Copyright (C) 2015 Kaveh Haghighi Mood
! Copyright (C) 2020 Vladimir Terzi
!
! SPDX-License-Identifier: GPL-3.0-or-later

! F90 module for handling the atomic basis functions (AO's)

! use allocatable arrays for uao etc.
! factored out basis input and norm factors
! allocation, input and output moved to aosData_m
! here only calculations based on aosData_m

! Modifications: //TODO Translate
!     SM
!     Gaussianausgaben als AO-Eingabe moeglich
!     Basisfunktionen aus einem GTO werden nicht mehr ueber splines
!     berechnet sondern analytisch, das ist effektiver weil hier keine
!     Schleife ueber die Kontraktion benoetigt wird
!     Nur noch eine Spline-Tabelle fuer identische Basisfunktionen bei
!     gleichen Kernen, Zugriff erfolgt ueber Zuordnungstabelle so(basmax)

module aos_m
    use kinds_m, only : r8
    use wfData_m
    use aosData_m, only: typ, bl, norm, bzet, ngto, cntrctn, &
        uao, uxao, uyao, uzao, u2ao, &
        so, mAOElecConfigs
    use cspline_m
    use eConfigs_m

    implicit none

    ! constants:
    character(len = *), parameter :: AO_TYPES = 'SPDFGHIK'
    real(r8), parameter :: SQRT_3 = SQRT(3._r8), SQRT_5 = SQRT(5._r8), &
        SQRT_7 = SQRT(7._r8), SQRT_5_3 = SQRT_5 / SQRT_3

contains
    subroutine aocalc(elecIdx, elecConfs, distTensor)
        ! aocalc calculates all atomic orbitals for position vectors x,y,z
        ! (in configurations elecConfs) including all required derivatives.
        ! This is the version that calculates several electron configurations
        ! Currently sequential calculation of electron configurations.

        ! on entry: requires rai distance matrix properly dimensioned and set
        !           for all input position vector x,y,z.
        !           for elecIdx>0 position vector must be modified only compared to
        !           last call to aocalc only at electron elecIdx
        ! on exit:  updates AO data structure in aosData_m (i.e. AOs and derivatives)

        ! Version 3.0 (03.04.2020)  extended STO and CGTO calculation for all AOs
        !                           lower than H
        ! Version 2.1 (28.04.1998)  references to walker eliminated
        ! Version 2.0 (20.03.1998)  module version
        ! Version 1.1 (25.06.1996)  deleted 5D code; one orb coeff for
        !                           all p,d,f orbs; simultaneous evaluation
        !                           for p,d,f orbs. contracted GTO's added.
        ! Version 1.0 (26.01.1996)

        ! input parameters:
        integer, intent(in) :: elecIdx  ! if > 0 only AO's for electron elecIdx recalculated
        type(eConfigArray), intent(inout) :: elecConfs  ! electron configurations (intent is in)
        real(r8), intent(in) :: distTensor(:, :, :)  ! electron-nucleus distances
        ! variables:
        integer :: w, i, initial, final, j, bf, atomIdx, k, n, l
        integer :: fxxx, fyyy, fzzz, &
            fxyy, fxxy, fxxz, fxzz, fyzz, fyyz, &
            fxyz
        integer :: gxxxx, gyyyy, gzzzz, &
            gxxxy, gxxxz, gyyyz, gxyyy, gxzzz, gyzzz, &
            gxxyy, gxxzz, gyyzz, &
            gxxyz, gxyyz, gxyzz
        real(r8), pointer :: xs(:), ys(:), zs(:)
        real(r8) :: xi, yi, zi, r, r2, ir2
        real(r8) :: u, ur, u1r, ur1r, u2r, ur2r, t, tu, tur
        real(r8) :: urx, ury, urz, &
            urx2, ury2, urz2, urxy, urxz, uryz, &
            urx3, ury3, urz3, urxyz
        real(r8) :: x, y, z, &
            x2, y2, z2, xy, xz, yz, &
            x3, y3, z3, x2y, x2z, y2z, xy2, xz2, yz2, xyz, &
            x4, y4, z4, x3y, x3z, y3z, xy3, xz3, yz3, x2y2, x2z2, y2z2, x2yz, xy2z, xyz2

        ! bf refers to the degenerate set of cartesian
        ! basis function (S:1,P:3,D:6,F:10) as input, which may be of type STO
        ! or contracted GTO.
        ! j refers to the individual basis function, as used in LCAO-MO's.
        ! (composed in subroutine mdetwf)
        ! i refers to the current electron.

        !-----Calculation of the AO's and their derivatives

        mAOElecConfigs = eConfigArray_size(elecConfs)
        do w = 1, mAOElecConfigs
            call eConfigArray_getPtr(elecConfs, w, xs, ys, zs)

            ! check input data for NaN or Inf
            call assert(ALL(distTensor(:, :, w) < HUGE(1._r8)), &
                "aocalc: illegal distTensor values")
            call assert(ALL(ABS(xs) < HUGE(1._r8)) .and. &
                ALL(ABS(ys) < HUGE(1._r8)) .and. &
                ALL(ABS(zs) < HUGE(1._r8)), &
                "aocalc: illegal x, y, z coords in elecConfs")

            if (evfmt == 'gau' .or. evfmt == 'mol') then  ! Gaussian order
                fxxx = 0
                fyyy = 1
                fzzz = 2
                fxyy = 3
                fxxy = 4
                fxxz = 5
                fxzz = 6
                fyzz = 7
                fyyz = 8
                fxyz = 9
                gzzzz = 0
                gyzzz = 1
                gyyzz = 2
                gyyyz = 3
                gyyyy = 4
                gxzzz = 5
                gxyzz = 6
                gxyyz = 7
                gxyyy = 8
                gxxzz = 9
                gxxyz = 10
                gxxyy = 11
                gxxxz = 12
                gxxxy = 13
                gxxxx = 14
            else  ! GAMESS (US) and TURBOMOLE order
                fxxx = 0
                fyyy = 1
                fzzz = 2
                fxxy = 3
                fxxz = 4
                fxyy = 5
                fyyz = 6
                fxzz = 7
                fyzz = 8
                fxyz = 9
                gxxxx = 0
                gyyyy = 1
                gzzzz = 2
                gxxxy = 3
                gxxxz = 4
                gxyyy = 5
                gyyyz = 6
                gxzzz = 7
                gyzzz = 8
                gxxyy = 9
                gxxzz = 10
                gyyzz = 11
                gxxyz = 12
                gxyyz = 13
                gxyzz = 14
            end if

            if (elecIdx == 0) then  ! AO's for all electrons
                initial = 1
                final = ne
            else  ! only AO for electron elecIdx
                initial = elecIdx
                final = elecIdx
            end if

            do i = initial, final  ! loop over electrons
                xi = xs(i)
                yi = ys(i)
                zi = zs(i)
                j = 1
                do bf = 1, nbasf
                    atomIdx = bc(bf)
                    x = xi - atoms(atomIdx)%cx
                    y = yi - atoms(atomIdx)%cy
                    z = zi - atoms(atomIdx)%cz
                    r = distTensor(atomIdx, i, w)
                    r2 = r * r
                    ir2 = 1 / r2
                    l = INDEX(AO_TYPES, bl(bf))
                    if (typ(bf) == 'STO') then
                        n = bn(bf) - l
                        t = - bzet(bf) * r
                        ur = norm(bf) * r ** n * EXP(t)
                        u = ur * ir2
                        t = t + n
                        u2r = u * (t * (t + 2 * l) - n)
                        u = u * t
                    else
                        u = 0
                        ur = 0
                        u2r = 0
                        do k = 1, ngto(bf)
                            t = - cntrctn(1, k, bf) * r2
                            tur = cntrctn(2, k, bf) * EXP(t)
                            t = t * 2
                            tu = tur * t
                            u = u + tu
                            ur = ur + tur
                            u2r = u2r + tu * (t + 2 * l + 1)
                        end do
                        u = u * ir2
                        u2r = u2r * ir2
                    end if
                    select case (l)
                    case (1)
                        n = j
                        uao(n, i, w) = ur
                        uxao(n, i, w) = u * x
                        uyao(n, i, w) = u * y
                        uzao(n, i, w) = u * z
                        u2ao(n, i, w) = u2r
                    case (2)
                        n = j
                        u1r = u * x
                        uao(n, i, w) = ur * x
                        uxao(n, i, w) = u1r * x + ur
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * x
                        n = j + 1
                        u1r = u * y
                        uao(n, i, w) = ur * y
                        uxao(n, i, w) = uyao(j, i, w)
                        uyao(n, i, w) = u1r * y + ur
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * y
                        n = j + 2
                        u1r = u * z
                        uao(n, i, w) = ur * z
                        uxao(n, i, w) = uzao(j, i, w)
                        uyao(n, i, w) = uzao(j + 1, i, w)
                        uzao(n, i, w) = u1r * z + ur
                        u2ao(n, i, w) = u2r * z
                    case (3)
                        x2 = x * x
                        y2 = y * y
                        z2 = z * z
                        xy = x * y
                        xz = x * z
                        yz = y * z

                        ur1r = 2 * ur
                        ur2r = ur1r
                        n = j
                        u1r = u * x2
                        uao(n, i, w) = ur * x2
                        uxao(n, i, w) = (u1r + ur1r) * x
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * x2 + ur2r
                        n = j + 1
                        u1r = u * y2
                        uao(n, i, w) = ur * y2
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = (u1r + ur1r) * y
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * y2 + ur2r
                        n = j + 2
                        u1r = u * z2
                        uao(n, i, w) = ur * z2
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = (u1r + ur1r) * z
                        u2ao(n, i, w) = u2r * z2 + ur2r

                        u = u * SQRT_3
                        ur = ur * SQRT_3
                        u2r = u2r * SQRT_3
                        urx = ur * x
                        ury = ur * y
                        urz = ur * z
                        n = j + 3
                        u1r = u * xy
                        uao(n, i, w) = ur * xy
                        uxao(n, i, w) = u1r * x + ury
                        uyao(n, i, w) = u1r * y + urx
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * xy
                        n = j + 4
                        u1r = u * xz
                        uao(n, i, w) = ur * xz
                        uxao(n, i, w) = u1r * x + urz
                        uyao(n, i, w) = uzao(j + 3, i, w)
                        uzao(n, i, w) = u1r * z + urx
                        u2ao(n, i, w) = u2r * xz
                        n = j + 5
                        u1r = u * yz
                        uao(n, i, w) = ur * yz
                        uxao(n, i, w) = uzao(j + 3, i, w)
                        uyao(n, i, w) = u1r * y + urz
                        uzao(n, i, w) = u1r * z + ury
                        u2ao(n, i, w) = u2r * yz
                    case (4)
                        x2 = x * x
                        y2 = y * y
                        z2 = z * z
                        xy = x * y
                        xz = x * z
                        yz = y * z
                        x3 = x2 * x
                        y3 = y2 * y
                        z3 = z2 * z
                        x2y = x2 * y
                        x2z = x2 * z
                        y2z = y2 * z
                        xy2 = xy * y
                        xz2 = xz * z
                        yz2 = yz * z
                        xyz = xy * z

                        ur1r = 3 * ur
                        ur2r = 6 * ur
                        n = j + fxxx
                        u1r = u * x3
                        uao(n, i, w) = ur * x3
                        uxao(n, i, w) = u1r * x + ur1r * x2
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * x3 + ur2r * x
                        n = j + fyyy
                        u1r = u * y3
                        uao(n, i, w) = ur * y3
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = u1r * y + ur1r * y2
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * y3 + ur2r * y
                        n = j + fzzz
                        u1r = u * z3
                        uao(n, i, w) = ur * z3
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z + ur1r * z2
                        u2ao(n, i, w) = u2r * z3 + ur2r * z

                        u = u * SQRT_5
                        ur = ur * SQRT_5
                        u2r = u2r * SQRT_5
                        ur1r = 2 * ur
                        ur2r = ur1r
                        urx2 = ur * x2
                        ury2 = ur * y2
                        urz2 = ur * z2
                        urxy = ur1r * xy
                        urxz = ur1r * xz
                        uryz = ur1r * yz
                        n = j + fxxy
                        u1r = u * x2y
                        uao(n, i, w) = ur * x2y
                        uxao(n, i, w) = u1r * x + urxy
                        uyao(n, i, w) = u1r * y + urx2
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * x2y + ur2r * y
                        n = j + fxxz
                        u1r = u * x2z
                        uao(n, i, w) = ur * x2z
                        uxao(n, i, w) = u1r * x + urxz
                        uyao(n, i, w) = uzao(j + fxxy, i, w)
                        uzao(n, i, w) = u1r * z + urx2
                        u2ao(n, i, w) = u2r * x2z + ur2r * z
                        n = j + fyyz
                        u1r = u * y2z
                        uao(n, i, w) = ur * y2z
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = u1r * y + uryz
                        uzao(n, i, w) = u1r * z + ury2
                        u2ao(n, i, w) = u2r * y2z + ur2r * z
                        n = j + fxyy
                        u1r = u * xy2
                        uao(n, i, w) = ur * xy2
                        uxao(n, i, w) = u1r * x + ury2
                        uyao(n, i, w) = u1r * y + urxy
                        uzao(n, i, w) = uxao(j + fyyz, i, w)
                        u2ao(n, i, w) = u2r * xy2 + ur2r * x
                        n = j + fxzz
                        u1r = u * xz2
                        uao(n, i, w) = ur * xz2
                        uxao(n, i, w) = u1r * x + urz2
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z + urxz
                        u2ao(n, i, w) = u2r * xz2 + ur2r * x
                        n = j + fyzz
                        u1r = u * yz2
                        uao(n, i, w) = ur * yz2
                        uxao(n, i, w) = uyao(j + fxzz, i, w)
                        uyao(n, i, w) = u1r * y + urz2
                        uzao(n, i, w) = u1r * z + uryz
                        u2ao(n, i, w) = u2r * yz2 + ur2r * y

                        u = u * SQRT_3
                        ur = ur * SQRT_3
                        u2r = u2r * SQRT_3
                        n = j + fxyz
                        u1r = u * xyz
                        uao(n, i, w) = ur * xyz
                        uxao(n, i, w) = u1r * x + ur * yz
                        uyao(n, i, w) = u1r * y + ur * xz
                        uzao(n, i, w) = u1r * z + ur * xy
                        u2ao(n, i, w) = u2r * xyz
                    case (5)
                        x2 = x * x
                        y2 = y * y
                        z2 = z * z
                        xy = x * y
                        xz = x * z
                        yz = y * z
                        x3 = x2 * x
                        y3 = y2 * y
                        z3 = z2 * z
                        x2y = x2 * y
                        x2z = x2 * z
                        y2z = y2 * z
                        xy2 = xy * y
                        xz2 = xz * z
                        yz2 = yz * z
                        xyz = xy * z
                        x4 = x3 * x
                        y4 = y3 * y
                        z4 = z3 * z
                        x3y = x3 * y
                        x3z = x3 * z
                        y3z = y3 * z
                        xy3 = xy2 * y
                        xz3 = xz2 * z
                        yz3 = yz2 * z
                        x2y2 = x2y * y
                        x2z2 = x2z * z
                        y2z2 = y2z * z
                        x2yz = x2y * z
                        xy2z = xy2 * z
                        xyz2 = xyz * z

                        ur1r = 4 * ur
                        ur2r = 12 * ur
                        n = j + gxxxx
                        u1r = u * x4
                        uao(n, i, w) = ur * x4
                        uxao(n, i, w) = u1r * x + ur1r * x3
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * x4 + ur2r * x2
                        n = j + gyyyy
                        u1r = u * y4
                        uao(n, i, w) = ur * y4
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = u1r * y + ur1r * y3
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * y4 + ur2r * y2
                        n = j + gzzzz
                        u1r = u * z4
                        uao(n, i, w) = ur * z4
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z + ur1r * z3
                        u2ao(n, i, w) = u2r * z4 + ur2r * z2

                        u = u * SQRT_7
                        ur = ur * SQRT_7
                        u2r = u2r * SQRT_7
                        ur1r = 3 * ur
                        ur2r = 6 * ur
                        urx3 = ur * x3
                        ury3 = ur * y3
                        urz3 = ur * z3
                        n = j + gxxxy
                        u1r = u * x3y
                        uao(n, i, w) = ur * x3y
                        uxao(n, i, w) = u1r * x + ur1r * x2y
                        uyao(n, i, w) = u1r * y + urx3
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * x3y + ur2r * xy
                        n = j + gxxxz
                        u1r = u * x3z
                        uao(n, i, w) = ur * x3z
                        uxao(n, i, w) = u1r * x + ur1r * x2z
                        uyao(n, i, w) = uzao(j + gxxxy, i, w)
                        uzao(n, i, w) = u1r * z + urx3
                        u2ao(n, i, w) = u2r * x3z + ur2r * xz
                        n = j + gyyyz
                        u1r = u * y3z
                        uao(n, i, w) = ur * y3z
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = u1r * y + ur1r * y2z
                        uzao(n, i, w) = u1r * z + ury3
                        u2ao(n, i, w) = u2r * y3z + ur2r * yz
                        n = j + gxyyy
                        u1r = u * xy3
                        uao(n, i, w) = ur * xy3
                        uxao(n, i, w) = u1r * x + ury3
                        uyao(n, i, w) = u1r * y + ur1r * xy2
                        uzao(n, i, w) = uxao(j + gyyyz, i, w)
                        u2ao(n, i, w) = u2r * xy3 + ur2r * xy
                        n = j + gxzzz
                        u1r = u * xz3
                        uao(n, i, w) = ur * xz3
                        uxao(n, i, w) = u1r * x + urz3
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z + ur1r * xz2
                        u2ao(n, i, w) = u2r * xz3 + ur2r * xz
                        n = j + gyzzz
                        u1r = u * yz3
                        uao(n, i, w) = ur * yz3
                        uxao(n, i, w) = uyao(j + gxzzz, i, w)
                        uyao(n, i, w) = u1r * y + urz3
                        uzao(n, i, w) = u1r * z + ur1r * yz2
                        u2ao(n, i, w) = u2r * yz3 + ur2r * yz

                        u = u * SQRT_5_3
                        ur = ur * SQRT_5_3
                        u2r = u2r * SQRT_5_3
                        ur1r = 2 * ur
                        ur2r = ur1r
                        n = j + gxxyy
                        u1r = u * x2y2
                        uao(n, i, w) = ur * x2y2
                        uxao(n, i, w) = u1r * x + ur1r * xy2
                        uyao(n, i, w) = u1r * y + ur1r * x2y
                        uzao(n, i, w) = u1r * z
                        u2ao(n, i, w) = u2r * x2y2 + ur2r * (x2 + y2)
                        n = j + gxxzz
                        u1r = u * x2z2
                        uao(n, i, w) = ur * x2z2
                        uxao(n, i, w) = u1r * x + ur1r * xz2
                        uyao(n, i, w) = u1r * y
                        uzao(n, i, w) = u1r * z + ur1r * x2z
                        u2ao(n, i, w) = u2r * x2z2 + ur2r * (x2 + z2)
                        n = j + gyyzz
                        u1r = u * y2z2
                        uao(n, i, w) = ur * y2z2
                        uxao(n, i, w) = u1r * x
                        uyao(n, i, w) = u1r * y + ur1r * yz2
                        uzao(n, i, w) = u1r * z + ur1r * y2z
                        u2ao(n, i, w) = u2r * y2z2 + ur2r * (y2 + z2)

                        u = u * SQRT_3
                        ur = ur * SQRT_3
                        u2r = u2r * SQRT_3
                        ur1r = 2 * ur
                        ur2r = ur1r
                        urxyz = ur1r * xyz
                        n = j + gxxyz
                        u1r = u * x2yz
                        uao(n, i, w) = ur * x2yz
                        uxao(n, i, w) = u1r * x + urxyz
                        uyao(n, i, w) = u1r * y + ur * x2z
                        uzao(n, i, w) = u1r * z + ur * x2y
                        u2ao(n, i, w) = u2r * x2yz + ur2r * yz
                        n = j + gxyyz
                        u1r = u * xy2z
                        uao(n, i, w) = ur * xy2z
                        uxao(n, i, w) = u1r * x + ur * y2z
                        uyao(n, i, w) = u1r * y + urxyz
                        uzao(n, i, w) = u1r * z + ur * xy2
                        u2ao(n, i, w) = u2r * xy2z + ur2r * xz
                        n = j + gxyzz
                        u1r = u * xyz2
                        uao(n, i, w) = ur * xyz2
                        uxao(n, i, w) = u1r * x + ur * yz2
                        uyao(n, i, w) = u1r * y + ur * xz2
                        uzao(n, i, w) = u1r * z + urxyz
                        u2ao(n, i, w) = u2r * xyz2 + ur2r * xy
                    case default
                        call abortp("aocalc: AOs higher than G are not implemented")
                    end select
                    j = j + l * (l + 1) / 2
                end do
            end do
        end do
    end subroutine aocalc

    ! ==================================================
    !BH

    !     --------------------------------
    subroutine aosplcalc(ie, ec, rrai)
        !     --------------------------------

        ! aosplcalc calculates all atomic orbitals for position vector x,y,z
        ! including all required derivatives, using splines.

        ! on entry: requires rai distance matrix properly dimensioned and set
        !           for input position vector x,y,z.
        !           for ie>0 position vector must be modified only compared to
        !           last call to aocalc only at electron ie
        ! on exit:  updates AO data structure in aos.h (i.e. AOs and derivatives)

        ! Version 1.1 (28.4.98):  reference to walker eliminated.
        ! Version 1.0 (22.11.97): AO calculation using splined radial parts
        !                         modified Version 2.0 of aocalc (1.1 of getaos)
        !                         splines are used only for contracted GTO's
        integer, intent(in) :: ie                 ! if >0 only AO's for electron ie recalculated
        type(eConfigArray), intent(inout) :: ec                 ! electron configurations
        real(r8), intent(in) :: rrai(:, :, :)     ! r_ai electron-nucleus distances
        ! constants:
        real(r8) sqr3, sqr5, sqr7
        parameter (sqr3 = 1.73205080756887729d0, sqr5 = 2.236067977499789696d0, &
                sqr7 = 2.645751311064591d0)
        ! variables
        integer bf, a, i, j, i1, i2, nn, al, ispl, n, w
        real(r8) rai(ncenter, ne)
        real(r8), pointer :: x(:), y(:), z(:)
        real(r8) xx, yy, zz, rr, r2, alp, nrm, u, ux, uxx, u2, dx, dy, dz, tmp, &
                dx2, dy2, dz2, dxyz, df
        logical gaussFOrder       ! .t.: Gaussian order for f function
        !                               ! .f.: Gamess==Turbomole order used


        ! bf refers to the degenerate set of cartesian
        ! basis function (S:1,P:3,D:6,F:10) as input, which may be of typ STO
        ! or contracted GTO.
        ! al refers to the individual basis function, as used in LCAO-MO's.
        ! (composed in subroutine mdetwf)
        ! i refers to the current electron.
        !-----Calculation of the AO's and their derivatives

        mAOElecConfigs = eConfigArray_size(ec)
        do w = 1, mAOElecConfigs

            call eConfigArray_getPtr(ec, w, x, y, z)
            rai = rrai(:, :, w)

            !     ! check input data for NaN or Inf
            call assert(all(rai<huge(1.d0)), "aosplcalc: illegal rai values")
            call assert(all(abs(x(1:ne))<huge(1.d0)) .and.&
                    all(abs(y(1:ne))<huge(1.d0)) .and. all(abs(z(1:ne))<huge(1.d0)), &
                    "aosplcalc: illegal x,y,z coords in ec")

            if (evfmt=='gau' .or. evfmt=='mol') then
                gaussFOrder = .true.
            else
                gaussFOrder = .false.
            end if

            if (ie == 0) then                     ! AO's for all electrons
                i1 = 1
                i2 = ne
            else
                i1 = ie                              ! only AO for electron ie
                i2 = ie
            end if

            do i = i1, i2                              ! loop over electrons

                xx = x(i)
                yy = y(i)
                zz = z(i)

                al = 1

                !        ! loop over basis functions:
                do n = 1, nbasf
                    bf = n
                    a = bc(bf)                         ! center of AO
                    rr = rai(a, i)                      ! r_ai
                    nn = bn(bf)                        ! n quantum no. of AO

                    if (typ(bf) == 'STO') then      ! STO-type basis function

                        alp = bzet(bf)                    ! orbital exponent
                        nrm = norm(bf)                    ! normalization constant

                        if (nn == 1) then                    ! 1s orbital
                            u = nrm * exp(-alp * rr)
                            ux = -alp * u / rr
                            uao(al, i, w) = u
                            uxao(al, i, w) = ux * (xx - atoms(a)%cx)
                            uyao(al, i, w) = ux * (yy - atoms(a)%cy)
                            uzao(al, i, w) = ux * (zz - atoms(a)%cz)
                            u2ao(al, i, w) = alp * (alp - 2d0 / rr) * u
                            al = al + 1
                        else if (nn == 2) then
                            if (bl(bf) == 'S') then           ! 2s orbital
                                !                 // 2s orbital
                                u = nrm * exp(-alp * rr)
                                ux = (1d0 - alp * rr) * u / rr
                                uao(al, i, w) = rr * u
                                uxao(al, i, w) = ux * (xx - atoms(a)%cx)
                                uyao(al, i, w) = ux * (yy - atoms(a)%cy)
                                uzao(al, i, w) = ux * (zz - atoms(a)%cz)
                                u2ao(al, i, w) = alp * (alp * rr - 2d0) * u + 2d0 * ux
                                al = al + 1
                            else if (bl(bf) == 'P') then      ! 2p orbital
                                u = nrm * exp(-alp * rr)
                                dx = xx - atoms(a)%cx
                                dy = yy - atoms(a)%cy
                                dz = zz - atoms(a)%cz
                                ux = -alp * u / rr
                                uao(al, i, w) = u * dx
                                uao(al + 1, i, w) = u * dy
                                uao(al + 2, i, w) = u * dz
                                uxao(al, i, w) = ux * dx * dx + u
                                uxao(al + 1, i, w) = ux * dx * dy
                                uxao(al + 2, i, w) = ux * dx * dz
                                uyao(al, i, w) = uxao(al + 1, i, w)
                                uyao(al + 1, i, w) = ux * dy * dy + u
                                uyao(al + 2, i, w) = ux * dy * dz
                                uzao(al, i, w) = uxao(al + 2, i, w)
                                uzao(al + 1, i, w) = uyao(al + 2, i, w)
                                uzao(al + 2, i, w) = ux * dz * dz + u
                                tmp = alp * (alp - 4d0 / rr) * u
                                u2ao(al, i, w) = tmp * dx
                                u2ao(al + 1, i, w) = tmp * dy
                                u2ao(al + 2, i, w) = tmp * dz
                                al = al + 3
                            else
                                call abortp(' aossplcalc: n=2 and l > 2')
                            end if

                        else if (nn == 3) then
                            if (bl(bf) == 'S') then
                                !                 // 3s orbital
                                u = nrm * exp(-alp * rr)
                                ux = (2d0 - alp * rr) * u
                                uao(al, i, w) = rr * rr * u
                                uxao(al, i, w) = ux * (xx - atoms(a)%cx)
                                uyao(al, i, w) = ux * (yy - atoms(a)%cy)
                                uzao(al, i, w) = ux * (zz - atoms(a)%cz)
                                tmp = (2d0 + alp * rr * (-4d0 + alp * rr)) * u
                                u2ao(al, i, w) = tmp + 2d0 * ux
                                al = al + 1
                            else if (bl(bf) == 'P') then
                                call abortp(' aosplcalc: 3p orbitals not implemented')
                            else if (bl(bf) == 'D') then         ! 3d orbital (6D)
                                !                 // Norm is different for d_xx and d_xy !

                                u = nrm * exp(-alp * rr)
                                dx = xx - atoms(a)%cx
                                dx2 = dx * dx
                                dy = yy - atoms(a)%cy
                                dy2 = dy * dy
                                dz = zz - atoms(a)%cz
                                dz2 = dz * dz
                                ux = -alp * u / rr
                                tmp = alp * (alp - 6d0 / rr)

                                uao(al, i, w) = dx2 * u
                                uao(al + 1, i, w) = dy2 * u
                                uao(al + 2, i, w) = dz2 * u
                                uxao(al, i, w) = dx * (2d0 * u + ux * dx2)
                                uxao(al + 1, i, w) = ux * dy2 * dx
                                uxao(al + 2, i, w) = ux * dz2 * dx
                                uyao(al, i, w) = ux * dx2 * dy
                                uyao(al + 1, i, w) = dy * (2d0 * u + ux * dy2)
                                uyao(al + 2, i, w) = ux * dz2 * dy
                                uzao(al, i, w) = ux * dx2 * dz
                                uzao(al + 1, i, w) = ux * dy2 * dz
                                uzao(al + 2, i, w) = dz * (2d0 * u + ux * dz2)
                                u2ao(al, i, w) = u * (2d0 + tmp * dx2)
                                u2ao(al + 1, i, w) = u * (2d0 + tmp * dy2)
                                u2ao(al + 2, i, w) = u * (2d0 + tmp * dz2)

                                u = sqr3 * u                   ! correction of norm for last 3
                                ux = sqr3 * ux

                                uao(al + 3, i, w) = u * dx * dy
                                uao(al + 4, i, w) = u * dx * dz
                                uao(al + 5, i, w) = u * dy * dz

                                uxao(al + 3, i, w) = dy * (u + ux * dx2)
                                uxao(al + 4, i, w) = dz * (u + ux * dx2)
                                uxao(al + 5, i, w) = ux * dx * dy * dz
                                uyao(al + 3, i, w) = dx * (u + ux * dy2)
                                uyao(al + 4, i, w) = uxao(al + 5, i, w)
                                uyao(al + 5, i, w) = dz * (u + ux * dy2)
                                uzao(al + 3, i, w) = uxao(al + 5, i, w)
                                uzao(al + 4, i, w) = dx * (u + ux * dz2)
                                uzao(al + 5, i, w) = dy * (u + ux * dz2)

                                u2ao(al + 3, i, w) = u * dx * dy * tmp
                                u2ao(al + 4, i, w) = u * dx * dz * tmp
                                u2ao(al + 5, i, w) = u * dy * dz * tmp
                                al = al + 6
                            end if     ! bl
                        end if        ! nn


                        !          // Contracted GTO's as basis function (AO)
                        !          // Evaluate with splines
                    else
                        if (so(bf)==0) then !only 1 GTO in contraction, no splines used !

                            r2 = rr * rr

                            !           // only primitive cartesian gaussians: 1s,2p,3d,4f
                            !           // i.e. no r factor. Thus nn is not used here.

                            if (bl(bf) == 'S') then                 ! 1s GTO
                                alp = cntrctn(1, 1, bf)
                                u = cntrctn(2, 1, bf) * exp(-alp * r2)
                                ux = -2d0 * alp * u
                                uao(al, i, w) = u
                                uxao(al, i, w) = ux * (xx - atoms(a)%cx)
                                uyao(al, i, w) = ux * (yy - atoms(a)%cy)
                                uzao(al, i, w) = ux * (zz - atoms(a)%cz)
                                u2ao(al, i, w) = ux * (3d0 - 2d0 * alp * r2)
                                al = al + 1

                            else if (bl(bf) == 'P') then             ! 2p GTO's
                                !              // do all 3 P simultaneously (same exponent is required)
                                !              // order p_x,p_y,p_z
                                alp = cntrctn(1, 1, bf)
                                u = cntrctn(2, 1, bf) * exp(-alp * r2)
                                dx = xx - atoms(a)%cx
                                dy = yy - atoms(a)%cy
                                dz = zz - atoms(a)%cz
                                ux = -2d0 * alp * u
                                uao(al, i, w) = dx * u
                                uao(al + 1, i, w) = dy * u
                                uao(al + 2, i, w) = dz * u
                                uxao(al, i, w) = u + ux * dx * dx
                                uxao(al + 1, i, w) = ux * dx * dy
                                uxao(al + 2, i, w) = ux * dx * dz
                                uyao(al, i, w) = ux * dx * dy
                                uyao(al + 1, i, w) = u + ux * dy * dy
                                uyao(al + 2, i, w) = ux * dy * dz
                                uzao(al, i, w) = ux * dx * dz
                                uzao(al + 1, i, w) = ux * dy * dz
                                uzao(al + 2, i, w) = u + ux * dz * dz
                                tmp = (5d0 - 2d0 * alp * r2) * ux
                                u2ao(al, i, w) = tmp * dx
                                u2ao(al + 1, i, w) = tmp * dy
                                u2ao(al + 2, i, w) = tmp * dz
                                al = al + 3

                            else if (bl(bf) == 'D') then         ! 3d GTO
                                !              // do all 6 D simultaneously (same exponent is required)
                                !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
                                alp = cntrctn(1, 1, bf)
                                u = cntrctn(2, 1, bf) * exp(-alp * r2)
                                dx = xx - atoms(a)%cx
                                dx2 = dx * dx
                                dy = yy - atoms(a)%cy
                                dy2 = dy * dy
                                dz = zz - atoms(a)%cz
                                dz2 = dz * dz
                                ux = -2d0 * alp * u
                                uao(al, i, w) = dx2 * u
                                uao(al + 1, i, w) = dy2 * u
                                uao(al + 2, i, w) = dz2 * u
                                uxao(al, i, w) = (2d0 * u + ux * dx2) * dx
                                uxao(al + 1, i, w) = dy2 * ux * dx
                                uxao(al + 2, i, w) = dz2 * ux * dx
                                uyao(al, i, w) = dx2 * ux * dy
                                uyao(al + 1, i, w) = (2d0 * u + ux * dy2) * dy
                                uyao(al + 2, i, w) = dz2 * ux * dy
                                uzao(al, i, w) = dx2 * ux * dz
                                uzao(al + 1, i, w) = dy2 * ux * dz
                                uzao(al + 2, i, w) = (2d0 * u + ux * dz2) * dz
                                tmp = (7d0 - 2d0 * alp * r2) * ux
                                u2ao(al, i, w) = 2d0 * u + dx2 * tmp
                                u2ao(al + 1, i, w) = 2d0 * u + dy2 * tmp
                                u2ao(al + 2, i, w) = 2d0 * u + dz2 * tmp
                                u = sqr3 * u                   ! correction of norm for last 3
                                ux = sqr3 * ux
                                uao(al + 3, i, w) = dx * dy * u
                                uao(al + 4, i, w) = dx * dz * u
                                uao(al + 5, i, w) = dy * dz * u
                                tmp = ux * dx * dy * dz
                                uxao(al + 3, i, w) = (u + ux * dx2) * dy
                                uxao(al + 4, i, w) = (u + ux * dx2) * dz
                                uxao(al + 5, i, w) = tmp
                                uyao(al + 3, i, w) = (u + ux * dy2) * dx
                                uyao(al + 4, i, w) = tmp
                                uyao(al + 5, i, w) = (u + ux * dy2) * dz
                                uzao(al + 3, i, w) = tmp
                                uzao(al + 4, i, w) = (u + ux * dz2) * dx
                                uzao(al + 5, i, w) = (u + ux * dz2) * dy
                                tmp = (7d0 - 2d0 * alp * r2) * ux
                                u2ao(al + 3, i, w) = tmp * dx * dy
                                u2ao(al + 4, i, w) = tmp * dx * dz
                                u2ao(al + 5, i, w) = tmp * dy * dz
                                al = al + 6

                            else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! 3f GTO
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
                                !              // f_xxx, f_yyy, f_zzz
                                uao(al, i, w) = dx2 * dx * u
                                uao(al + 1, i, w) = dy2 * dy * u
                                uao(al + 2, i, w) = dz2 * dz * u
                                uxao(al, i, w) = (3d0 * u + ux * dx2) * dx2
                                uxao(al + 1, i, w) = dy2 * dy * ux * dx
                                uxao(al + 2, i, w) = dz2 * dz * ux * dx
                                uyao(al, i, w) = dx2 * dx * ux * dy
                                uyao(al + 1, i, w) = (3d0 * u + ux * dy2) * dy2
                                uyao(al + 2, i, w) = dz2 * dz * ux * dy
                                uzao(al, i, w) = dx2 * dx * ux * dz
                                uzao(al + 1, i, w) = dy2 * dy * ux * dz
                                uzao(al + 2, i, w) = (3d0 * u + ux * dz2) * dz2
                                tmp = (9d0 - 2d0 * alp * r2) * ux
                                u2ao(al, i, w) = (6d0 * u + dx2 * tmp) * dx
                                u2ao(al + 1, i, w) = (6d0 * u + dy2 * tmp) * dy
                                u2ao(al + 2, i, w) = (6d0 * u + dz2 * tmp) * dz
                                !              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                                u = sqr5 * u                   ! correction of norm
                                ux = sqr5 * ux
                                uao(al + 3, i, w) = dx2 * dy * u
                                uao(al + 4, i, w) = dx2 * dz * u
                                uao(al + 5, i, w) = dy2 * dx * u
                                uao(al + 6, i, w) = dy2 * dz * u
                                uao(al + 7, i, w) = dz2 * dx * u
                                uao(al + 8, i, w) = dz2 * dy * u

                                tmp = ux * dxyz
                                uxao(al + 3, i, w) = (2d0 * u + ux * dx2) * dx * dy
                                uxao(al + 4, i, w) = (2d0 * u + ux * dx2) * dx * dz
                                uxao(al + 5, i, w) = (u + ux * dx2) * dy2
                                uxao(al + 6, i, w) = tmp * dy
                                uxao(al + 7, i, w) = (u + ux * dx2) * dz2
                                uxao(al + 8, i, w) = tmp * dz
                                uyao(al + 3, i, w) = (u + ux * dy2) * dx2
                                uyao(al + 4, i, w) = tmp * dx
                                uyao(al + 5, i, w) = (2d0 * u + ux * dy2) * dx * dy
                                uyao(al + 6, i, w) = (2d0 * u + ux * dy2) * dy * dz
                                uyao(al + 7, i, w) = tmp * dz
                                uyao(al + 8, i, w) = (u + ux * dy2) * dz2
                                uzao(al + 3, i, w) = tmp * dx
                                uzao(al + 4, i, w) = (u + ux * dz2) * dx2
                                uzao(al + 5, i, w) = tmp * dy
                                uzao(al + 6, i, w) = (u + ux * dz2) * dy2
                                uzao(al + 7, i, w) = (2d0 * u + ux * dz2) * dx * dz
                                uzao(al + 8, i, w) = (2d0 * u + ux * dz2) * dy * dz

                                tmp = (9d0 - 2d0 * alp * r2) * ux
                                u2ao(al + 3, i, w) = (2d0 * u + dx2 * tmp) * dy
                                u2ao(al + 4, i, w) = (2d0 * u + dx2 * tmp) * dz
                                u2ao(al + 5, i, w) = (2d0 * u + dy2 * tmp) * dx
                                u2ao(al + 6, i, w) = (2d0 * u + dy2 * tmp) * dz
                                u2ao(al + 7, i, w) = (2d0 * u + dz2 * tmp) * dx
                                u2ao(al + 8, i, w) = (2d0 * u + dz2 * tmp) * dy
                                !              // f_xyz
                                u = sqr3 * u                  ! correction of norm
                                ux = sqr3 * ux
                                uao(al + 9, i, w) = dxyz * u
                                uxao(al + 9, i, w) = (u + ux * dx2) * dy * dz
                                uyao(al + 9, i, w) = (u + ux * dy2) * dx * dz
                                uzao(al + 9, i, w) = (u + ux * dz2) * dx * dy
                                tmp = (9d0 - 2d0 * alp * r2) * ux
                                u2ao(al + 9, i, w) = dxyz * tmp
                                al = al + 10

                            else if (bl(bf)=='F'.and.gaussFOrder) then     ! 3f GTO
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
                                !              // f_xxx, f_yyy, f_zzz
                                uao(al, i, w) = dx2 * dx * u
                                uao(al + 1, i, w) = dy2 * dy * u
                                uao(al + 2, i, w) = dz2 * dz * u
                                uxao(al, i, w) = (3d0 * u + ux * dx2) * dx2
                                uxao(al + 1, i, w) = dy2 * dy * ux * dx
                                uxao(al + 2, i, w) = dz2 * dz * ux * dx
                                uyao(al, i, w) = dx2 * dx * ux * dy
                                uyao(al + 1, i, w) = (3d0 * u + ux * dy2) * dy2
                                uyao(al + 2, i, w) = dz2 * dz * ux * dy
                                uzao(al, i, w) = dx2 * dx * ux * dz
                                uzao(al + 1, i, w) = dy2 * dy * ux * dz
                                uzao(al + 2, i, w) = (3d0 * u + ux * dz2) * dz2
                                tmp = (9d0 - 2d0 * alp * r2) * ux
                                u2ao(al, i, w) = (6d0 * u + dx2 * tmp) * dx
                                u2ao(al + 1, i, w) = (6d0 * u + dy2 * tmp) * dy
                                u2ao(al + 2, i, w) = (6d0 * u + dz2 * tmp) * dz
                                !              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                                u = sqr5 * u                   ! correction of norm
                                ux = sqr5 * ux

                                uao(al + 4, i, w) = dx2 * dy * u
                                uao(al + 5, i, w) = dx2 * dz * u
                                uao(al + 3, i, w) = dy2 * dx * u
                                uao(al + 8, i, w) = dy2 * dz * u
                                uao(al + 6, i, w) = dz2 * dx * u
                                uao(al + 7, i, w) = dz2 * dy * u

                                tmp = ux * dxyz
                                uxao(al + 4, i, w) = (2d0 * u + ux * dx2) * dx * dy
                                uxao(al + 5, i, w) = (2d0 * u + ux * dx2) * dx * dz
                                uxao(al + 3, i, w) = (u + ux * dx2) * dy2
                                uxao(al + 8, i, w) = tmp * dy
                                uxao(al + 6, i, w) = (u + ux * dx2) * dz2
                                uxao(al + 7, i, w) = tmp * dz
                                uyao(al + 4, i, w) = (u + ux * dy2) * dx2
                                uyao(al + 5, i, w) = tmp * dx
                                uyao(al + 3, i, w) = (2d0 * u + ux * dy2) * dx * dy
                                uyao(al + 8, i, w) = (2d0 * u + ux * dy2) * dy * dz
                                uyao(al + 6, i, w) = tmp * dz
                                uyao(al + 7, i, w) = (u + ux * dy2) * dz2
                                uzao(al + 4, i, w) = tmp * dx
                                uzao(al + 5, i, w) = (u + ux * dz2) * dx2
                                uzao(al + 3, i, w) = tmp * dy
                                uzao(al + 8, i, w) = (u + ux * dz2) * dy2
                                uzao(al + 6, i, w) = (2d0 * u + ux * dz2) * dx * dz
                                uzao(al + 7, i, w) = (2d0 * u + ux * dz2) * dy * dz

                                tmp = (9d0 - 2d0 * alp * r2) * ux
                                u2ao(al + 4, i, w) = (2d0 * u + dx2 * tmp) * dy
                                u2ao(al + 5, i, w) = (2d0 * u + dx2 * tmp) * dz
                                u2ao(al + 3, i, w) = (2d0 * u + dy2 * tmp) * dx
                                u2ao(al + 8, i, w) = (2d0 * u + dy2 * tmp) * dz
                                u2ao(al + 6, i, w) = (2d0 * u + dz2 * tmp) * dx
                                u2ao(al + 7, i, w) = (2d0 * u + dz2 * tmp) * dy

                                !              // f_xyz
                                u = sqr3 * u                  ! correction of norm
                                ux = sqr3 * ux
                                uao(al + 9, i, w) = dxyz * u
                                uxao(al + 9, i, w) = (u + ux * dx2) * dy * dz
                                uyao(al + 9, i, w) = (u + ux * dy2) * dx * dz
                                uzao(al + 9, i, w) = (u + ux * dz2) * dx * dy
                                tmp = (9d0 - 2d0 * alp * r2) * ux
                                u2ao(al + 9, i, w) = dxyz * tmp
                                al = al + 10

                            else if (bl(bf)=='G') then     ! 5g GTO
                                !              // do all 15 cartesian G simultaneously (same exponent is required)

                                uao(al:al + 14, i, w) = 0d0
                                uxao(al:al + 14, i, w) = 0d0
                                uyao(al:al + 14, i, w) = 0d0
                                uzao(al:al + 14, i, w) = 0d0
                                u2ao(al:al + 14, i, w) = 0d0

                                if (gaussFOrder) then
                                    call internal_GaussianOrderGFunctions()
                                else
                                    call internal_GamessOrderGFunctions()
                                end if

                                al = al + 15

                            else
                                call abortp('(aossplcalc): wrong GTO')
                            end if  ! bl

                        else     !more then 1 GTO in contraction, splines used !
                            r2 = rr * rr
                            j = (csplnpnt - 1) * rr / (csalpha + rr) + 1
                            df = rr - csplx(j)

                            !           // only primitive cartesian gaussians: 1s,2p,3d,4f
                            !           // i.e. no r factor. Thus nn is not used here.

                            if (bl(bf) == 'S') then                 ! 1s GTO

                                ispl = 3 * so(bf) - 2
                                uao(al, i, w) = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                ux = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                u2 = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                uxao(al, i, w) = ux * (xx - atoms(a)%cx) / rr
                                uyao(al, i, w) = ux * (yy - atoms(a)%cy) / rr
                                uzao(al, i, w) = ux * (zz - atoms(a)%cz) / rr
                                u2ao(al, i, w) = u2 + 2 * ux / rr

                                al = al + 1

                            else if (bl(bf) == 'P') then             ! 2p GTO's

                                !              // do all 3 P simultaneously (same exponent is required)
                                !              // order p_x,p_y,p_z
                                ispl = 3 * so(bf) - 2
                                u = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                ux = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                uxx = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                                dx = xx - atoms(a)%cx
                                dy = yy - atoms(a)%cy
                                dz = zz - atoms(a)%cz
                                uao(al, i, w) = dx * u
                                uao(al + 1, i, w) = dy * u
                                uao(al + 2, i, w) = dz * u
                                uxao(al, i, w) = u + ux * dx * dx
                                uxao(al + 1, i, w) = ux * dx * dy
                                uxao(al + 2, i, w) = ux * dx * dz
                                uyao(al, i, w) = ux * dx * dy
                                uyao(al + 1, i, w) = u + ux * dy * dy
                                uyao(al + 2, i, w) = ux * dy * dz
                                uzao(al, i, w) = ux * dx * dz
                                uzao(al + 1, i, w) = ux * dy * dz
                                uzao(al + 2, i, w) = u + ux * dz * dz
                                u2ao(al, i, w) = uxx * dx
                                u2ao(al + 1, i, w) = uxx * dy
                                u2ao(al + 2, i, w) = uxx * dz
                                al = al + 3

                            else if (bl(bf) == 'D') then         ! 3d GTO

                                !              // do all 6 D simultaneously (same exponent is required)
                                !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
                                ispl = 3 * so(bf) - 2
                                u = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                ux = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                uxx = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                                dx = xx - atoms(a)%cx
                                dx2 = dx * dx
                                dy = yy - atoms(a)%cy
                                dy2 = dy * dy
                                dz = zz - atoms(a)%cz
                                dz2 = dz * dz
                                uao(al, i, w) = dx2 * u
                                uao(al + 1, i, w) = dy2 * u
                                uao(al + 2, i, w) = dz2 * u
                                uxao(al, i, w) = (2d0 * u + ux * dx2) * dx
                                uxao(al + 1, i, w) = dy2 * ux * dx
                                uxao(al + 2, i, w) = dz2 * ux * dx
                                uyao(al, i, w) = dx2 * ux * dy
                                uyao(al + 1, i, w) = (2d0 * u + ux * dy2) * dy
                                uyao(al + 2, i, w) = dz2 * ux * dy
                                uzao(al, i, w) = dx2 * ux * dz
                                uzao(al + 1, i, w) = dy2 * ux * dz
                                uzao(al + 2, i, w) = (2d0 * u + ux * dz2) * dz
                                u2ao(al, i, w) = 2d0 * u + dx2 * uxx
                                u2ao(al + 1, i, w) = 2d0 * u + dy2 * uxx
                                u2ao(al + 2, i, w) = 2d0 * u + dz2 * uxx

                                u = sqr3 * u                   ! correction of norm for last 3
                                ux = sqr3 * ux
                                uxx = sqr3 * uxx

                                uao(al + 3, i, w) = dx * dy * u
                                uao(al + 4, i, w) = dx * dz * u
                                uao(al + 5, i, w) = dy * dz * u

                                tmp = ux * dx * dy * dz
                                uxao(al + 3, i, w) = (u + ux * dx2) * dy
                                uxao(al + 4, i, w) = (u + ux * dx2) * dz
                                uxao(al + 5, i, w) = tmp
                                uyao(al + 3, i, w) = (u + ux * dy2) * dx
                                uyao(al + 4, i, w) = tmp
                                uyao(al + 5, i, w) = (u + ux * dy2) * dz
                                uzao(al + 3, i, w) = tmp
                                uzao(al + 4, i, w) = (u + ux * dz2) * dx
                                uzao(al + 5, i, w) = (u + ux * dz2) * dy

                                u2ao(al + 3, i, w) = uxx * dx * dy
                                u2ao(al + 4, i, w) = uxx * dx * dz
                                u2ao(al + 5, i, w) = uxx * dy * dz
                                al = al + 6

                            else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! 3f GTO
                                !              // do all 10 F simultaneously (same exponent is required)
                                !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                                !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

                                ispl = 3 * so(bf) - 2
                                u = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                ux = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                uxx = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                                dx = xx - atoms(a)%cx
                                dx2 = dx * dx
                                dy = yy - atoms(a)%cy
                                dy2 = dy * dy
                                dz = zz - atoms(a)%cz
                                dz2 = dz * dz
                                dxyz = dx * dy * dz

                                !              // f_xxx, f_yyy, f_zzz
                                uao(al, i, w) = dx2 * dx * u
                                uao(al + 1, i, w) = dy2 * dy * u
                                uao(al + 2, i, w) = dz2 * dz * u

                                uxao(al, i, w) = (3d0 * u + ux * dx2) * dx2
                                uxao(al + 1, i, w) = dy2 * dy * ux * dx
                                uxao(al + 2, i, w) = dz2 * dz * ux * dx
                                uyao(al, i, w) = dx2 * dx * ux * dy
                                uyao(al + 1, i, w) = (3d0 * u + ux * dy2) * dy2
                                uyao(al + 2, i, w) = dz2 * dz * ux * dy
                                uzao(al, i, w) = dx2 * dx * ux * dz
                                uzao(al + 1, i, w) = dy2 * dy * ux * dz
                                uzao(al + 2, i, w) = (3d0 * u + ux * dz2) * dz2
                                u2ao(al, i, w) = (6d0 * u + dx2 * uxx) * dx
                                u2ao(al + 1, i, w) = (6d0 * u + dy2 * uxx) * dy
                                u2ao(al + 2, i, w) = (6d0 * u + dz2 * uxx) * dz

                                !              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                                u = sqr5 * u                   ! correction of norm
                                ux = sqr5 * ux
                                uxx = sqr5 * uxx

                                uao(al + 3, i, w) = dx2 * dy * u
                                uao(al + 4, i, w) = dx2 * dz * u
                                uao(al + 5, i, w) = dy2 * dx * u
                                uao(al + 6, i, w) = dy2 * dz * u
                                uao(al + 7, i, w) = dz2 * dx * u
                                uao(al + 8, i, w) = dz2 * dy * u
                                tmp = ux * dxyz
                                uxao(al + 3, i, w) = (2d0 * u + ux * dx2) * dx * dy
                                uxao(al + 4, i, w) = (2d0 * u + ux * dx2) * dx * dz
                                uxao(al + 5, i, w) = (u + ux * dx2) * dy2
                                uxao(al + 6, i, w) = tmp * dy
                                uxao(al + 7, i, w) = (u + ux * dx2) * dz2
                                uxao(al + 8, i, w) = tmp * dz
                                uyao(al + 3, i, w) = (u + ux * dy2) * dx2
                                uyao(al + 4, i, w) = tmp * dx
                                uyao(al + 5, i, w) = (2d0 * u + ux * dy2) * dx * dy
                                uyao(al + 6, i, w) = (2d0 * u + ux * dy2) * dy * dz
                                uyao(al + 7, i, w) = tmp * dz
                                uyao(al + 8, i, w) = (u + ux * dy2) * dz2
                                uzao(al + 3, i, w) = tmp * dx
                                uzao(al + 4, i, w) = (u + ux * dz2) * dx2
                                uzao(al + 5, i, w) = tmp * dy
                                uzao(al + 6, i, w) = (u + ux * dz2) * dy2
                                uzao(al + 7, i, w) = (2d0 * u + ux * dz2) * dx * dz
                                uzao(al + 8, i, w) = (2d0 * u + ux * dz2) * dy * dz

                                u2ao(al + 3, i, w) = (2d0 * u + dx2 * uxx) * dy
                                u2ao(al + 4, i, w) = (2d0 * u + dx2 * uxx) * dz
                                u2ao(al + 5, i, w) = (2d0 * u + dy2 * uxx) * dx
                                u2ao(al + 6, i, w) = (2d0 * u + dy2 * uxx) * dz
                                u2ao(al + 7, i, w) = (2d0 * u + dz2 * uxx) * dx
                                u2ao(al + 8, i, w) = (2d0 * u + dz2 * uxx) * dy

                                !              // f_xyz
                                u = sqr3 * u                  ! correction of norm
                                ux = sqr3 * ux
                                uxx = sqr3 * uxx

                                uao(al + 9, i, w) = dxyz * u

                                uxao(al + 9, i, w) = (u + ux * dx2) * dy * dz
                                uyao(al + 9, i, w) = (u + ux * dy2) * dx * dz
                                uzao(al + 9, i, w) = (u + ux * dz2) * dx * dy
                                u2ao(al + 9, i, w) = dxyz * uxx
                                al = al + 10

                            else if (bl(bf)=='F'.and.gaussFOrder) then     ! f GTO
                                !              // do all 10 F simultaneously (same exponent is required)
                                !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                                !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

                                ispl = 3 * so(bf) - 2
                                u = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                ux = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))
                                ispl = ispl + 1
                                uxx = cspla(ispl, j) + df * (csplb(ispl, j)&
                                        + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                                dx = xx - atoms(a)%cx
                                dx2 = dx * dx
                                dy = yy - atoms(a)%cy
                                dy2 = dy * dy
                                dz = zz - atoms(a)%cz
                                dz2 = dz * dz
                                dxyz = dx * dy * dz

                                !              // f_xxx, f_yyy, f_zzz
                                uao(al, i, w) = dx2 * dx * u
                                uao(al + 1, i, w) = dy2 * dy * u
                                uao(al + 2, i, w) = dz2 * dz * u

                                uxao(al, i, w) = (3d0 * u + ux * dx2) * dx2
                                uxao(al + 1, i, w) = dy2 * dy * ux * dx
                                uxao(al + 2, i, w) = dz2 * dz * ux * dx
                                uyao(al, i, w) = dx2 * dx * ux * dy
                                uyao(al + 1, i, w) = (3d0 * u + ux * dy2) * dy2
                                uyao(al + 2, i, w) = dz2 * dz * ux * dy
                                uzao(al, i, w) = dx2 * dx * ux * dz
                                uzao(al + 1, i, w) = dy2 * dy * ux * dz
                                uzao(al + 2, i, w) = (3d0 * u + ux * dz2) * dz2
                                u2ao(al, i, w) = (6d0 * u + dx2 * uxx) * dx
                                u2ao(al + 1, i, w) = (6d0 * u + dy2 * uxx) * dy
                                u2ao(al + 2, i, w) = (6d0 * u + dz2 * uxx) * dz

                                !              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                                u = sqr5 * u                   ! correction of norm
                                ux = sqr5 * ux
                                uxx = sqr5 * uxx

                                uao(al + 4, i, w) = dx2 * dy * u
                                uao(al + 5, i, w) = dx2 * dz * u
                                uao(al + 3, i, w) = dy2 * dx * u
                                uao(al + 8, i, w) = dy2 * dz * u
                                uao(al + 6, i, w) = dz2 * dx * u
                                uao(al + 7, i, w) = dz2 * dy * u
                                tmp = ux * dxyz
                                uxao(al + 4, i, w) = (2d0 * u + ux * dx2) * dx * dy
                                uxao(al + 5, i, w) = (2d0 * u + ux * dx2) * dx * dz
                                uxao(al + 3, i, w) = (u + ux * dx2) * dy2
                                uxao(al + 8, i, w) = tmp * dy
                                uxao(al + 6, i, w) = (u + ux * dx2) * dz2
                                uxao(al + 7, i, w) = tmp * dz
                                uyao(al + 4, i, w) = (u + ux * dy2) * dx2
                                uyao(al + 5, i, w) = tmp * dx
                                uyao(al + 3, i, w) = (2d0 * u + ux * dy2) * dx * dy
                                uyao(al + 8, i, w) = (2d0 * u + ux * dy2) * dy * dz
                                uyao(al + 6, i, w) = tmp * dz
                                uyao(al + 7, i, w) = (u + ux * dy2) * dz2
                                uzao(al + 4, i, w) = tmp * dx
                                uzao(al + 5, i, w) = (u + ux * dz2) * dx2
                                uzao(al + 3, i, w) = tmp * dy
                                uzao(al + 8, i, w) = (u + ux * dz2) * dy2
                                uzao(al + 6, i, w) = (2d0 * u + ux * dz2) * dx * dz
                                uzao(al + 7, i, w) = (2d0 * u + ux * dz2) * dy * dz

                                u2ao(al + 4, i, w) = (2d0 * u + dx2 * uxx) * dy
                                u2ao(al + 5, i, w) = (2d0 * u + dx2 * uxx) * dz
                                u2ao(al + 3, i, w) = (2d0 * u + dy2 * uxx) * dx
                                u2ao(al + 8, i, w) = (2d0 * u + dy2 * uxx) * dz
                                u2ao(al + 6, i, w) = (2d0 * u + dz2 * uxx) * dx
                                u2ao(al + 7, i, w) = (2d0 * u + dz2 * uxx) * dy

                                !              // f_xyz
                                u = sqr3 * u                  ! correction of norm
                                ux = sqr3 * ux
                                uxx = sqr3 * uxx

                                uao(al + 9, i, w) = dxyz * u

                                uxao(al + 9, i, w) = (u + ux * dx2) * dy * dz
                                uyao(al + 9, i, w) = (u + ux * dy2) * dx * dz
                                uzao(al + 9, i, w) = (u + ux * dz2) * dx * dy
                                u2ao(al + 9, i, w) = dxyz * uxx
                                al = al + 10

                            else
                                call abortp('(aosplcalc): wrong GTO')
                            end if       ! bl
                        end if    ! simple GTO, no splines
                    end if        ! STO/GTO
                end do          ! bf-loop over basis functions
            end do             ! i-loop over electrons

        end do ! w-loop over elec configs


    contains

        subroutine internal_GamessOrderGFunctions()

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

            !           // g_xxxx, g_yyyy, g_zzzz
            uao(al, i, w) = uao(al, i, w) + dx2 * dx2 * u
            uao(al + 1, i, w) = uao(al + 1, i, w) + dy2 * dy2 * u
            uao(al + 2, i, w) = uao(al + 2, i, w) + dz2 * dz2 * u

            uxao(al, i, w) = uxao(al, i, w) + (4d0 * u + ux * dx2) * dx2 * dx
            uxao(al + 1, i, w) = uxao(al + 1, i, w) + dy2 * dy2 * ux * dx
            uxao(al + 2, i, w) = uxao(al + 2, i, w) + dz2 * dz2 * ux * dx
            uyao(al, i, w) = uyao(al, i, w) + dx2 * dx2 * ux * dy
            uyao(al + 1, i, w) = uyao(al + 1, i, w) + (4d0 * u + ux * dy2) * dy2 * dy
            uyao(al + 2, i, w) = uyao(al + 2, i, w) + dz2 * dz2 * ux * dy
            uzao(al, i, w) = uzao(al, i, w) + dx2 * dx2 * ux * dz
            uzao(al + 1, i, w) = uzao(al + 1, i, w) + dy2 * dy2 * ux * dz
            uzao(al + 2, i, w) = uzao(al + 2, i, w) + (4d0 * u + ux * dz2) * dz2 * dz
            tmp = (11d0 - 2d0 * alp * r2) * ux
            u2ao(al, i, w) = u2ao(al, i, w) + (12d0 * u + dx2 * tmp) * dx2
            u2ao(al + 1, i, w) = u2ao(al + 1, i, w) + (12d0 * u + dy2 * tmp) * dy2
            u2ao(al + 2, i, w) = u2ao(al + 2, i, w) + (12d0 * u + dz2 * tmp) * dz2

            !           // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
            u = sqr7 * u                   ! correction of norm
            ux = sqr7 * ux

            uao(al + 3, i, w) = uao(al + 3, i, w) + dx2 * dx * dy * u
            uao(al + 4, i, w) = uao(al + 4, i, w) + dx2 * dx * dz * u
            uao(al + 5, i, w) = uao(al + 5, i, w) + dy2 * dy * dx * u
            uao(al + 6, i, w) = uao(al + 6, i, w) + dy2 * dy * dz * u
            uao(al + 7, i, w) = uao(al + 7, i, w) + dz2 * dz * dx * u
            uao(al + 8, i, w) = uao(al + 8, i, w) + dz2 * dz * dy * u

            tmp = ux * dxyz
            uxao(al + 3, i, w) = uxao(al + 3, i, w) + (3d0 * u + ux * dx2) * dx2 * dy
            uxao(al + 4, i, w) = uxao(al + 4, i, w) + (3d0 * u + ux * dx2) * dx2 * dz
            uxao(al + 5, i, w) = uxao(al + 5, i, w) + (u + ux * dx2) * dy2 * dy
            uxao(al + 6, i, w) = uxao(al + 6, i, w) + tmp * dy2
            uxao(al + 7, i, w) = uxao(al + 7, i, w) + (u + ux * dx2) * dz2 * dz
            uxao(al + 8, i, w) = uxao(al + 8, i, w) + tmp * dz2
            uyao(al + 3, i, w) = uyao(al + 3, i, w) + (u + ux * dy2) * dx2 * dx
            uyao(al + 4, i, w) = uyao(al + 4, i, w) + tmp * dx2
            uyao(al + 5, i, w) = uyao(al + 5, i, w) + (3d0 * u + ux * dy2) * dy2 * dx
            uyao(al + 6, i, w) = uyao(al + 6, i, w) + (3d0 * u + ux * dy2) * dy2 * dz
            uyao(al + 7, i, w) = uyao(al + 7, i, w) + tmp * dz2
            uyao(al + 8, i, w) = uyao(al + 8, i, w) + (u + ux * dy2) * dz2 * dz
            uzao(al + 3, i, w) = uzao(al + 3, i, w) + tmp * dx2
            uzao(al + 4, i, w) = uzao(al + 4, i, w) + (u + ux * dz2) * dx2 * dx
            uzao(al + 5, i, w) = uzao(al + 5, i, w) + tmp * dy2
            uzao(al + 6, i, w) = uzao(al + 6, i, w) + (u + ux * dz2) * dy2 * dy
            uzao(al + 7, i, w) = uzao(al + 7, i, w) + (3d0 * u + ux * dz2) * dz2 * dx
            uzao(al + 8, i, w) = uzao(al + 8, i, w) + (3d0 * u + ux * dz2) * dz2 * dy

            tmp = (11d0 - 2d0 * alp * r2) * ux
            u2ao(al + 3, i, w) = u2ao(al + 3, i, w) + (6d0 * u + dx2 * tmp) * dx * dy
            u2ao(al + 4, i, w) = u2ao(al + 4, i, w) + (6d0 * u + dx2 * tmp) * dx * dz
            u2ao(al + 5, i, w) = u2ao(al + 5, i, w) + (6d0 * u + dy2 * tmp) * dy * dx
            u2ao(al + 6, i, w) = u2ao(al + 6, i, w) + (6d0 * u + dy2 * tmp) * dy * dz
            u2ao(al + 7, i, w) = u2ao(al + 7, i, w) + (6d0 * u + dz2 * tmp) * dz * dx
            u2ao(al + 8, i, w) = u2ao(al + 8, i, w) + (6d0 * u + dz2 * tmp) * dz * dy

            !           // g_xxyy, g_xxzz, g_yyzz
            u = sqr5 / sqr3 * u          ! correction of norm
            ux = sqr5 / sqr3 * ux

            uao(al + 9, i, w) = uao(al + 9, i, w) + dx2 * dy2 * u
            uao(al + 10, i, w) = uao(al + 10, i, w) + dx2 * dz2 * u
            uao(al + 11, i, w) = uao(al + 11, i, w) + dy2 * dz2 * u

            tmp = ux * dxyz
            uxao(al + 9, i, w) = uxao(al + 9, i, w) + (2d0 * u + ux * dx2) * dx * dy2
            uxao(al + 10, i, w) = uxao(al + 10, i, w) + (2d0 * u + ux * dx2) * dx * dz2
            uxao(al + 11, i, w) = uxao(al + 11, i, w) + tmp * dy * dz
            uyao(al + 9, i, w) = uyao(al + 9, i, w) + (2d0 * u + ux * dy2) * dy * dx2
            uyao(al + 10, i, w) = uyao(al + 10, i, w) + tmp * dx * dz
            uyao(al + 11, i, w) = uyao(al + 11, i, w) + (2d0 * u + ux * dy2) * dy * dz2
            uzao(al + 9, i, w) = uzao(al + 9, i, w) + tmp * dx * dy
            uzao(al + 10, i, w) = uzao(al + 10, i, w) + (2d0 * u + ux * dz2) * dz * dx2
            uzao(al + 11, i, w) = uzao(al + 11, i, w) + (2d0 * u + ux * dz2) * dz * dy2

            tmp = (11d0 - 2d0 * alp * r2) * ux
            u2ao(al + 9, i, w) = u2ao(al + 9, i, w) + 2d0 * u * (dx2 + dy2) + dx2 * dy2 * tmp
            u2ao(al + 10, i, w) = u2ao(al + 10, i, w) + 2d0 * u * (dx2 + dz2) + dx2 * dz2 * tmp
            u2ao(al + 11, i, w) = u2ao(al + 11, i, w) + 2d0 * u * (dy2 + dz2) + dy2 * dz2 * tmp


            !           // g_xxyz, g_yyxz, g_zzxy
            u = sqr3 * u                  ! correction of norm
            ux = sqr3 * ux

            uao(al + 12, i, w) = uao(al + 12, i, w) + dx * dxyz * u
            uao(al + 13, i, w) = uao(al + 13, i, w) + dy * dxyz * u
            uao(al + 14, i, w) = uao(al + 14, i, w) + dz * dxyz * u

            uxao(al + 12, i, w) = uxao(al + 12, i, w) + (2d0 * u + ux * dx2) * dxyz
            uxao(al + 13, i, w) = uxao(al + 13, i, w) + (u + ux * dx2) * dy2 * dz
            uxao(al + 14, i, w) = uxao(al + 14, i, w) + (u + ux * dx2) * dz2 * dy
            uyao(al + 12, i, w) = uyao(al + 12, i, w) + (u + ux * dy2) * dx2 * dz
            uyao(al + 13, i, w) = uyao(al + 13, i, w) + (2d0 * u + ux * dy2) * dxyz
            uyao(al + 14, i, w) = uyao(al + 14, i, w) + (u + ux * dy2) * dz2 * dx
            uzao(al + 12, i, w) = uzao(al + 12, i, w) + (u + ux * dz2) * dx2 * dy
            uzao(al + 13, i, w) = uzao(al + 13, i, w) + (u + ux * dz2) * dy2 * dx
            uzao(al + 14, i, w) = uzao(al + 14, i, w) + (2d0 * u + ux * dz2) * dxyz

            tmp = (11d0 - 2d0 * alp * r2) * ux
            u2ao(al + 12, i, w) = u2ao(al + 12, i, w) + (2d0 * u + tmp * dx2) * dy * dz
            u2ao(al + 13, i, w) = u2ao(al + 13, i, w) + (2d0 * u + tmp * dy2) * dx * dz
            u2ao(al + 14, i, w) = u2ao(al + 14, i, w) + (2d0 * u + tmp * dz2) * dx * dy

        end subroutine internal_GamessOrderGFunctions


        subroutine internal_GaussianOrderGFunctions()

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

            !           // g_xxxx, g_yyyy, g_zzzz
            uao(al + 14, i, w) = uao(al + 14, i, w) + dx2 * dx2 * u
            uao(al + 4, i, w) = uao(al + 4, i, w) + dy2 * dy2 * u
            uao(al, i, w) = uao(al, i, w) + dz2 * dz2 * u

            uxao(al + 14, i, w) = uxao(al + 14, i, w) + (4d0 * u + ux * dx2) * dx2 * dx
            uxao(al + 4, i, w) = uxao(al + 4, i, w) + dy2 * dy2 * ux * dx
            uxao(al, i, w) = uxao(al, i, w) + dz2 * dz2 * ux * dx
            uyao(al + 14, i, w) = uyao(al + 14, i, w) + dx2 * dx2 * ux * dy
            uyao(al + 4, i, w) = uyao(al + 4, i, w) + (4d0 * u + ux * dy2) * dy2 * dy
            uyao(al, i, w) = uyao(al, i, w) + dz2 * dz2 * ux * dy
            uzao(al + 14, i, w) = uzao(al + 14, i, w) + dx2 * dx2 * ux * dz
            uzao(al + 4, i, w) = uzao(al + 4, i, w) + dy2 * dy2 * ux * dz
            uzao(al, i, w) = uzao(al, i, w) + (4d0 * u + ux * dz2) * dz2 * dz
            tmp = (11d0 - 2d0 * alp * r2) * ux
            u2ao(al + 14, i, w) = u2ao(al + 14, i, w) + (12d0 * u + dx2 * tmp) * dx2
            u2ao(al + 4, i, w) = u2ao(al + 4, i, w) + (12d0 * u + dy2 * tmp) * dy2
            u2ao(al, i, w) = u2ao(al, i, w) + (12d0 * u + dz2 * tmp) * dz2

            !           // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
            u = sqr7 * u                   ! correction of norm
            ux = sqr7 * ux

            uao(al + 13, i, w) = uao(al + 13, i, w) + dx2 * dx * dy * u
            uao(al + 12, i, w) = uao(al + 12, i, w) + dx2 * dx * dz * u
            uao(al + 8, i, w) = uao(al + 8, i, w) + dy2 * dy * dx * u
            uao(al + 3, i, w) = uao(al + 3, i, w) + dy2 * dy * dz * u
            uao(al + 5, i, w) = uao(al + 5, i, w) + dz2 * dz * dx * u
            uao(al + 1, i, w) = uao(al + 1, i, w) + dz2 * dz * dy * u

            tmp = ux * dxyz
            uxao(al + 13, i, w) = uxao(al + 13, i, w) + (3d0 * u + ux * dx2) * dx2 * dy
            uxao(al + 12, i, w) = uxao(al + 12, i, w) + (3d0 * u + ux * dx2) * dx2 * dz
            uxao(al + 8, i, w) = uxao(al + 8, i, w) + (u + ux * dx2) * dy2 * dy
            uxao(al + 3, i, w) = uxao(al + 3, i, w) + tmp * dy2
            uxao(al + 5, i, w) = uxao(al + 5, i, w) + (u + ux * dx2) * dz2 * dz
            uxao(al + 1, i, w) = uxao(al + 1, i, w) + tmp * dz2
            uyao(al + 13, i, w) = uyao(al + 13, i, w) + (u + ux * dy2) * dx2 * dx
            uyao(al + 12, i, w) = uyao(al + 12, i, w) + tmp * dx2
            uyao(al + 8, i, w) = uyao(al + 8, i, w) + (3d0 * u + ux * dy2) * dy2 * dx
            uyao(al + 3, i, w) = uyao(al + 3, i, w) + (3d0 * u + ux * dy2) * dy2 * dz
            uyao(al + 5, i, w) = uyao(al + 5, i, w) + tmp * dz2
            uyao(al + 1, i, w) = uyao(al + 1, i, w) + (u + ux * dy2) * dz2 * dz
            uzao(al + 13, i, w) = uzao(al + 13, i, w) + tmp * dx2
            uzao(al + 12, i, w) = uzao(al + 12, i, w) + (u + ux * dz2) * dx2 * dx
            uzao(al + 8, i, w) = uzao(al + 8, i, w) + tmp * dy2
            uzao(al + 3, i, w) = uzao(al + 3, i, w) + (u + ux * dz2) * dy2 * dy
            uzao(al + 5, i, w) = uzao(al + 5, i, w) + (3d0 * u + ux * dz2) * dz2 * dx
            uzao(al + 1, i, w) = uzao(al + 1, i, w) + (3d0 * u + ux * dz2) * dz2 * dy

            tmp = (11d0 - 2d0 * alp * r2) * ux
            u2ao(al + 13, i, w) = u2ao(al + 13, i, w) + (6d0 * u + dx2 * tmp) * dx * dy
            u2ao(al + 12, i, w) = u2ao(al + 12, i, w) + (6d0 * u + dx2 * tmp) * dx * dz
            u2ao(al + 8, i, w) = u2ao(al + 8, i, w) + (6d0 * u + dy2 * tmp) * dy * dx
            u2ao(al + 3, i, w) = u2ao(al + 3, i, w) + (6d0 * u + dy2 * tmp) * dy * dz
            u2ao(al + 5, i, w) = u2ao(al + 5, i, w) + (6d0 * u + dz2 * tmp) * dz * dx
            u2ao(al + 1, i, w) = u2ao(al + 1, i, w) + (6d0 * u + dz2 * tmp) * dz * dy

            !           // g_xxyy, g_xxzz, g_yyzz
            u = sqr5 / sqr3 * u          ! correction of norm
            ux = sqr5 / sqr3 * ux

            uao(al + 11, i, w) = uao(al + 11, i, w) + dx2 * dy2 * u
            uao(al + 9, i, w) = uao(al + 9, i, w) + dx2 * dz2 * u
            uao(al + 2, i, w) = uao(al + 2, i, w) + dy2 * dz2 * u

            tmp = ux * dxyz
            uxao(al + 11, i, w) = uxao(al + 11, i, w) + (2d0 * u + ux * dx2) * dx * dy2
            uxao(al + 9, i, w) = uxao(al + 9, i, w) + (2d0 * u + ux * dx2) * dx * dz2
            uxao(al + 2, i, w) = uxao(al + 2, i, w) + tmp * dy * dz
            uyao(al + 11, i, w) = uyao(al + 11, i, w) + (2d0 * u + ux * dy2) * dy * dx2
            uyao(al + 9, i, w) = uyao(al + 9, i, w) + tmp * dx * dz
            uyao(al + 2, i, w) = uyao(al + 2, i, w) + (2d0 * u + ux * dy2) * dy * dz2
            uzao(al + 11, i, w) = uzao(al + 11, i, w) + tmp * dx * dy
            uzao(al + 9, i, w) = uzao(al + 9, i, w) + (2d0 * u + ux * dz2) * dz * dx2
            uzao(al + 2, i, w) = uzao(al + 2, i, w) + (2d0 * u + ux * dz2) * dz * dy2

            tmp = (11d0 - 2d0 * alp * r2) * ux
            u2ao(al + 11, i, w) = u2ao(al + 11, i, w) + 2d0 * u * (dx2 + dy2) + dx2 * dy2 * tmp
            u2ao(al + 9, i, w) = u2ao(al + 9, i, w) + 2d0 * u * (dx2 + dz2) + dx2 * dz2 * tmp
            u2ao(al + 2, i, w) = u2ao(al + 2, i, w) + 2d0 * u * (dy2 + dz2) + dy2 * dz2 * tmp


            !           // g_xxyz, g_yyxz, g_zzxy
            u = sqr3 * u                  ! correction of norm
            ux = sqr3 * ux

            uao(al + 10, i, w) = uao(al + 10, i, w) + dx * dxyz * u
            uao(al + 7, i, w) = uao(al + 7, i, w) + dy * dxyz * u
            uao(al + 6, i, w) = uao(al + 6, i, w) + dz * dxyz * u

            uxao(al + 10, i, w) = uxao(al + 10, i, w) + (2d0 * u + ux * dx2) * dxyz
            uxao(al + 7, i, w) = uxao(al + 7, i, w) + (u + ux * dx2) * dy2 * dz
            uxao(al + 6, i, w) = uxao(al + 6, i, w) + (u + ux * dx2) * dz2 * dy
            uyao(al + 10, i, w) = uyao(al + 10, i, w) + (u + ux * dy2) * dx2 * dz
            uyao(al + 7, i, w) = uyao(al + 7, i, w) + (2d0 * u + ux * dy2) * dxyz
            uyao(al + 6, i, w) = uyao(al + 6, i, w) + (u + ux * dy2) * dz2 * dx
            uzao(al + 10, i, w) = uzao(al + 10, i, w) + (u + ux * dz2) * dx2 * dy
            uzao(al + 7, i, w) = uzao(al + 7, i, w) + (u + ux * dz2) * dy2 * dx
            uzao(al + 6, i, w) = uzao(al + 6, i, w) + (2d0 * u + ux * dz2) * dxyz

            tmp = (11d0 - 2d0 * alp * r2) * ux
            u2ao(al + 10, i, w) = u2ao(al + 10, i, w) + (2d0 * u + tmp * dx2) * dy * dz
            u2ao(al + 7, i, w) = u2ao(al + 7, i, w) + (2d0 * u + tmp * dy2) * dx * dz
            u2ao(al + 6, i, w) = u2ao(al + 6, i, w) + (2d0 * u + tmp * dz2) * dx * dy

        end subroutine internal_GaussianOrderGFunctions


    end subroutine aosplcalc


    !================================================
    !BH

    !     -------------------------------
    subroutine ao1calc(ie, x, y, z, rai)
        !     -------------------------------

        ! aocalc calculates all atomic orbitals for position vector x,y,z
        ! without derivatives. For ie>0 update only

        ! on entry: requires rai distance matrix properly dimensioned and set
        !           for input position vector x,y,z.
        !           for ie>0 position vector must be modified only compared to
        !           last call only at electron ie
        ! on exit:  updates AO data structure in aos.h (i.e. AOs)

        ! Version 1.3 (28.4.98)    references to walker eliminated
        ! Version 1.1 (6/26/1996) simultaneous p,d calculation
        ! Version 1.0 (5/29/1996)

        ! input parameters:
        integer ie                 ! if >0 only AO's for electron ie recalculated
        real(r8)   x(:),& ! x,y,z coordinates of position vector&
                       y(:),&
                       z(:), &
                       rai(:, :)    ! r_ai electron-nucleus distances
        ! output (to local data structure (aos.h) in commons block):
        !     uao  : AO-array uao(n,i,w) for nth AO at electron i
        !EH
        ! constants:
        real(r8) sqr3, sqr5, sqr7
        parameter (sqr3 = 1.73205080756887729d0, sqr5 = 2.236067977499789696d0, &
                sqr7 = 2.645751311064591d0)
        ! variables
        !     integer al,bf,a,i,ii,i1,i2,nn,ic
        integer al, bf, a, i, ii, i1, i2, nn, ic, n
        real(r8) xx, yy, zz, rr, alp, nrm, u, dx, dy, dz, r2, dx2, dy2, dz2, dxyz
        !TS
        logical gaussFOrder       ! .t.: Gaussian order for f function
        !                               ! .f.: Gamess==Turbomole order used

        mAOElecConfigs = 1

        !-----Calculation of the AO's

        if (evfmt=='gau' .or. evfmt=='mol') then
            gaussFOrder = .true.
        else
            gaussFOrder = .false.
        end if

        if (ie == 0) then                     ! AO's for all electrons
            i1 = 1
            i2 = ne
        else
            i1 = ie                              ! only AO for electron ie
            i2 = ie
        end if

        do i = i1, i2                              ! loop over electrons

            xx = x(i)
            yy = y(i)
            zz = z(i)

            al = 1

            do n = 1, nbasf                   ! loop over basis functions
                bf = n
                a = bc(bf)                         ! center of AO
                rr = rai(a, i)                      ! r_ai
                nn = bn(bf)                        ! n quantum no. of AO

                if (typ(bf) == 'STO') then      ! STO-type basis function

                    alp = bzet(bf)                    ! orbital exponent
                    nrm = norm(bf)                    ! normalization constant

                    if (nn == 1) then               ! 1s orbital
                        u = nrm * exp(-alp * rr)
                        uao(al, i, 1) = u
                        al = al + 1
                    else if (nn == 2) then
                        if (bl(bf) == 'S') then           ! 2s orbital
                            !                 // 2s orbital
                            u = nrm * exp(-alp * rr)
                            uao(al, i, 1) = rr * u
                            al = al + 1
                        else if (bl(bf) == 'P') then      ! 2p orbital
                            u = nrm * exp(-alp * rr)
                            dx = xx - atoms(a)%cx
                            dy = yy - atoms(a)%cy
                            dz = zz - atoms(a)%cz
                            uao(al, i, 1) = u * dx
                            uao(al + 1, i, 1) = u * dy
                            uao(al + 2, i, 1) = u * dz
                            al = al + 3

                        else
                            call abortp(' ao1calc: n=2 and l > 2')
                        end if

                    else if (nn == 3) then
                        if (bl(bf) == 'S') then
                            !                 // 3s orbital
                            u = nrm * exp(-alp * rr)
                            uao(al, i, 1) = rr * rr * u
                            al = al + 1
                        else if (bl(bf) == 'P') then
                            call abortp('ao1calc: 3p orbitals not implemented')

                        else if (bl(bf) == 'D') then         ! 3d orbital (6D)
                            !                 // Norm is different for d_xx and d_xy !
                            u = nrm * exp(-alp * rr)
                            dx = xx - atoms(a)%cx
                            dx2 = dx * dx
                            dy = yy - atoms(a)%cy
                            dy2 = dy * dy
                            dz = zz - atoms(a)%cz
                            dz2 = dz * dz

                            uao(al, i, 1) = dx2 * u
                            uao(al + 1, i, 1) = dy2 * u
                            uao(al + 2, i, 1) = dz2 * u

                            u = sqr3 * u                   ! correction of norm for last 3
                            uao(al + 3, i, 1) = u * dx * dy
                            uao(al + 4, i, 1) = u * dx * dz
                            uao(al + 5, i, 1) = u * dy * dz
                            al = al + 6
                        end if     ! bl
                    end if        ! nn


                    !          // Contracted GTO's as basis function (AO)
                else

                    r2 = rr * rr

                    !           // only primitive cartesian gaussians: 1s,2p,3d,4f
                    !           // i.e. no r factor. Thus nn is not used here.

                    if (bl(bf) == 'S') then                 ! 1s GTO
                        uao(al, i, 1) = 0d0
                        do ic = 1, ngto(bf)                     ! loop over contraction
                            alp = cntrctn(1, ic, bf)
                            u = cntrctn(2, ic, bf) * exp(-alp * r2)
                            uao(al, i, 1) = uao(al, i, 1) + u
                        end do
                        al = al + 1

                    else if (bl(bf) == 'P') then             ! 2p GTO's
                        !              // do all 3 P simultaneously (same exponent is required)
                        !              // order p_x,p_y,p_z
                        do ii = 0, 2
                            uao(al + ii, i, 1) = 0d0
                        end do
                        do ic = 1, ngto(bf)                      ! loop over contraction
                            alp = cntrctn(1, ic, bf)
                            u = cntrctn(2, ic, bf) * exp(-alp * r2)
                            dx = xx - atoms(a)%cx
                            dy = yy - atoms(a)%cy
                            dz = zz - atoms(a)%cz
                            uao(al, i, 1) = uao(al, i, 1) + dx * u
                            uao(al + 1, i, 1) = uao(al + 1, i, 1) + dy * u
                            uao(al + 2, i, 1) = uao(al + 2, i, 1) + dz * u
                        end do
                        al = al + 3

                    else if (bl(bf) == 'D') then         ! 3d GTO
                        !              // do all 6 D simultaneously (same exponent is required)
                        !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
                        do ii = 0, 5
                            uao(al + ii, i, 1) = 0d0
                        end do
                        do ic = 1, ngto(bf)                      ! loop over contraction
                            alp = cntrctn(1, ic, bf)
                            u = cntrctn(2, ic, bf) * exp(-alp * r2)
                            dx = xx - atoms(a)%cx
                            dx2 = dx * dx
                            dy = yy - atoms(a)%cy
                            dy2 = dy * dy
                            dz = zz - atoms(a)%cz
                            dz2 = dz * dz

                            uao(al, i, 1) = uao(al, i, 1) + dx2 * u
                            uao(al + 1, i, 1) = uao(al + 1, i, 1) + dy2 * u
                            uao(al + 2, i, 1) = uao(al + 2, i, 1) + dz2 * u

                            u = sqr3 * u                   ! correction of norm for last 3

                            uao(al + 3, i, 1) = uao(al + 3, i, 1) + dx * dy * u
                            uao(al + 4, i, 1) = uao(al + 4, i, 1) + dx * dz * u
                            uao(al + 5, i, 1) = uao(al + 5, i, 1) + dy * dz * u
                        end do
                        al = al + 6

                    else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! f GTO
                        !              // do all 10 F simultaneously (same exponent is required)
                        !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                        !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)
                        do ii = 0, 9
                            uao(al + ii, i, 1) = 0d0
                        end do
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

                            !                 // f_xxx, f_yyy, f_zzz
                            uao(al, i, 1) = uao(al, i, 1) + dx2 * dx * u
                            uao(al + 1, i, 1) = uao(al + 1, i, 1) + dy2 * dy * u
                            uao(al + 2, i, 1) = uao(al + 2, i, 1) + dz2 * dz * u

                            !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                            u = sqr5 * u                   ! correction of norm

                            uao(al + 3, i, 1) = uao(al + 3, i, 1) + dx2 * dy * u
                            uao(al + 4, i, 1) = uao(al + 4, i, 1) + dx2 * dz * u
                            uao(al + 5, i, 1) = uao(al + 5, i, 1) + dy2 * dx * u
                            uao(al + 6, i, 1) = uao(al + 6, i, 1) + dy2 * dz * u
                            uao(al + 7, i, 1) = uao(al + 7, i, 1) + dz2 * dx * u
                            uao(al + 8, i, 1) = uao(al + 8, i, 1) + dz2 * dy * u
                            !                 // f_xyz
                            u = sqr3 * u                  ! correction of norm
                            uao(al + 9, i, 1) = uao(al + 9, i, 1) + dxyz * u
                        end do
                        al = al + 10

                    else if (bl(bf)=='F'.and.gaussFOrder) then     ! 3f GTO
                        !              // do all 10 F simultaneously (same exponent is required)
                        !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                        !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)
                        do ii = 0, 9
                            uao(al + ii, i, 1) = 0d0
                        end do
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

                            !                 // f_xxx, f_yyy, f_zzz
                            uao(al, i, 1) = uao(al, i, 1) + dx2 * dx * u
                            uao(al + 1, i, 1) = uao(al + 1, i, 1) + dy2 * dy * u
                            uao(al + 2, i, 1) = uao(al + 2, i, 1) + dz2 * dz * u

                            !                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                            u = sqr5 * u                   ! correction of norm

                            uao(al + 3, i, 1) = uao(al + 3, i, 1) + dy2 * dx * u
                            uao(al + 4, i, 1) = uao(al + 4, i, 1) + dx2 * dy * u
                            uao(al + 5, i, 1) = uao(al + 5, i, 1) + dx2 * dz * u
                            uao(al + 6, i, 1) = uao(al + 6, i, 1) + dz2 * dx * u
                            uao(al + 7, i, 1) = uao(al + 7, i, 1) + dz2 * dy * u
                            uao(al + 8, i, 1) = uao(al + 8, i, 1) + dy2 * dz * u
                            !                 // f_xyz
                            u = sqr3 * u                  ! correction of norm
                            uao(al + 9, i, 1) = uao(al + 9, i, 1) + dxyz * u
                        end do
                        al = al + 10

                    else if (bl(bf)=='G') then     ! 5g GTO
                        !              // do all 15 cartesian G simultaneously (same exponent is required)

                        uao(al:al + 14, i, 1) = 0d0

                        if (gaussFOrder) then
                            call internal_GaussianOrderGFunctions()
                        else
                            call internal_GamessOrderGFunctions()
                        end if

                        al = al + 15

                    else
                        call abortp('(ao1calc): wrong GTO')
                    end if  ! bl
                end if  ! STO/GTO
            end do  ! bf-loop over basis functions
        end do  ! i-loop over electrons


    contains

        subroutine internal_GamessOrderGFunctions()

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

                !              // g_xxxx, g_yyyy, g_zzzz
                uao(al, i, 1) = uao(al, i, 1) + dx2 * dx2 * u
                uao(al + 1, i, 1) = uao(al + 1, i, 1) + dy2 * dy2 * u
                uao(al + 2, i, 1) = uao(al + 2, i, 1) + dz2 * dz2 * u

                !              // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
                u = sqr7 * u                   ! correction of norm

                uao(al + 3, i, 1) = uao(al + 3, i, 1) + dx2 * dx * dy * u
                uao(al + 4, i, 1) = uao(al + 4, i, 1) + dx2 * dx * dz * u
                uao(al + 5, i, 1) = uao(al + 5, i, 1) + dy2 * dy * dx * u
                uao(al + 6, i, 1) = uao(al + 6, i, 1) + dy2 * dy * dz * u
                uao(al + 7, i, 1) = uao(al + 7, i, 1) + dz2 * dz * dx * u
                uao(al + 8, i, 1) = uao(al + 8, i, 1) + dz2 * dz * dy * u

                !              // g_xxyy, g_xxzz, g_yyzz
                u = sqr5 / sqr3 * u          ! correction of norm

                uao(al + 9, i, 1) = uao(al + 9, i, 1) + dx2 * dy2 * u
                uao(al + 10, i, 1) = uao(al + 10, i, 1) + dx2 * dz2 * u
                uao(al + 11, i, 1) = uao(al + 11, i, 1) + dy2 * dz2 * u

                !              // g_xxyz, g_yyxz, g_zzxy
                u = sqr3 * u                  ! correction of norm

                uao(al + 12, i, 1) = uao(al + 12, i, 1) + dx * dxyz * u
                uao(al + 13, i, 1) = uao(al + 13, i, 1) + dy * dxyz * u
                uao(al + 14, i, 1) = uao(al + 14, i, 1) + dz * dxyz * u
            end do

        end subroutine internal_GamessOrderGFunctions


        subroutine internal_GaussianOrderGFunctions()
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

                !              // g_xxxx, g_yyyy, g_zzzz
                uao(al + 14, i, 1) = uao(al + 14, i, 1) + dx2 * dx2 * u
                uao(al + 4, i, 1) = uao(al + 4, i, 1) + dy2 * dy2 * u
                uao(al, i, 1) = uao(al, i, 1) + dz2 * dz2 * u

                !              // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
                u = sqr7 * u                   ! correction of norm

                uao(al + 13, i, 1) = uao(al + 13, i, 1) + dx2 * dx * dy * u
                uao(al + 12, i, 1) = uao(al + 12, i, 1) + dx2 * dx * dz * u
                uao(al + 8, i, 1) = uao(al + 8, i, 1) + dy2 * dy * dx * u
                uao(al + 3, i, 1) = uao(al + 3, i, 1) + dy2 * dy * dz * u
                uao(al + 6, i, 1) = uao(al + 6, i, 1) + dz2 * dz * dx * u
                uao(al + 1, i, 1) = uao(al + 1, i, 1) + dz2 * dz * dy * u

                !              // g_xxyy, g_xxzz, g_yyzz
                u = sqr5 / sqr3 * u          ! correction of norm

                uao(al + 11, i, 1) = uao(al + 11, i, 1) + dx2 * dy2 * u
                uao(al + 9, i, 1) = uao(al + 9, i, 1) + dx2 * dz2 * u
                uao(al + 2, i, 1) = uao(al + 2, i, 1) + dy2 * dz2 * u

                !              // g_xxyz, g_yyxz, g_zzxy
                u = sqr3 * u                  ! correction of norm

                uao(al + 10, i, 1) = uao(al + 10, i, 1) + dx * dxyz * u
                uao(al + 7, i, 1) = uao(al + 7, i, 1) + dy * dxyz * u
                uao(al + 6, i, 1) = uao(al + 6, i, 1) + dz * dxyz * u
            end do
        end subroutine internal_GaussianOrderGFunctions

    end subroutine ao1calc

    ! ==================================================
    !BH

    !     ----------------------------------
    subroutine ao1splcalc(ie, x, y, z, rai)
        !     ----------------------------------

        ! ao1splcalc calculates all atomic orbitals for position vector x,y,z


        ! on entry: requires rai distance matrix properly dimensioned and set
        !           for input position vector x,y,z.
        !           for ie>0 position vector must be modified only compared to
        !           last call to aocalc only at electron ie

        ! input parameters:
        integer ie                 ! if >0 only AO's for electron ie recalculated
        real(r8)   x(:),& ! x,y,z coordinates of position vector&
                       y(:),&
                       z(:), &
                       rai(:, :)     ! r_ai electron-nucleus distances

        !EH
        ! constants:
        real(r8) sqr3, sqr5, sqr7
        parameter (sqr3 = 1.73205080756887729d0, sqr5 = 2.236067977499789696d0, &
                sqr7 = 2.645751311064591d0)
        ! variables
        !     integer bf,a,i,j,i1,i2,nn,al,ispl
        integer al, bf, a, i, j, i1, i2, n, ispl

        real(r8) xx, yy, zz, rr, r2, alp, u, dx, dy, dz, dx2, dy2, dz2, dxyz, df
        !TS
        logical gaussFOrder       ! .t.: Gaussian order for f function
        !                               ! .f.: Gamess==Turbomole order used


        ! bf refers to the degenerate set of cartesian
        ! basis function (S:1,P:3,D:6,F:10) as input, which may be of type STO
        ! or contracted GTO.
        ! al refers to the individual basis function, as used in LCAO-MO's.
        ! (composed in subroutine mdetwf)
        ! i refers to the current electron.
        !-----Calculation of the AO's and their derivatives

        mAOElecConfigs = 1

        if (evfmt=='gau' .or. evfmt=='mol') then
            gaussFOrder = .true.
        else
            gaussFOrder = .false.
        end if

        if (ie == 0) then                     ! AO's for all electrons
            i1 = 1
            i2 = ne
        else
            i1 = ie                              ! only AO for electron ie
            i2 = ie
        end if

        do i = i1, i2                              ! loop over electrons

            xx = x(i)
            yy = y(i)
            zz = z(i)

            al = 1

            !        ! loop over basis functions:
            do n = 1, nbasf
                bf = n
                a = bc(bf)                         ! center of AO
                rr = rai(a, i)                      ! r_ai
                !          ! nn = bn(bf)                        ! n quantum no. of AO

                if (typ(bf) == 'STO') then      ! STO-type basis function
                    call abortp('STOs make no sense in ao1splcalc')

                    !          // Contracted GTO's as basis function (AO)
                    !          // Evaluate with splines
                else
                    if (so(bf)==0) then !only 1 GTO in contraction, no splines used !

                        r2 = rr * rr

                        !           // only primitive cartesian gaussians: 1s,2p,3d,4f
                        !           // i.e. no r factor. Thus nn is not used here.

                        if (bl(bf) == 'S') then                 ! 1s GTO
                            alp = cntrctn(1, 1, bf)
                            u = cntrctn(2, 1, bf) * exp(-alp * r2)
                            uao(al, i, 1) = u
                            al = al + 1
                        else if (bl(bf) == 'P') then             ! 2p GTO's
                            !              // do all 3 P simultaneously (same exponent is required)
                            !              // order p_x,p_y,p_z
                            alp = cntrctn(1, 1, bf)
                            u = cntrctn(2, 1, bf) * exp(-alp * r2)
                            dx = xx - atoms(a)%cx
                            dy = yy - atoms(a)%cy
                            dz = zz - atoms(a)%cz
                            uao(al, i, 1) = dx * u
                            uao(al + 1, i, 1) = dy * u
                            uao(al + 2, i, 1) = dz * u

                            al = al + 3

                        else if (bl(bf) == 'D') then         ! 3d GTO
                            !              // do all 6 D simultaneously (same exponent is required)
                            !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
                            alp = cntrctn(1, 1, bf)
                            u = cntrctn(2, 1, bf) * exp(-alp * r2)
                            dx = xx - atoms(a)%cx
                            dx2 = dx * dx
                            dy = yy - atoms(a)%cy
                            dy2 = dy * dy
                            dz = zz - atoms(a)%cz
                            dz2 = dz * dz
                            uao(al, i, 1) = dx2 * u
                            uao(al + 1, i, 1) = dy2 * u
                            uao(al + 2, i, 1) = dz2 * u

                            u = sqr3 * u                   ! correction of norm for last 3

                            uao(al + 3, i, 1) = dx * dy * u
                            uao(al + 4, i, 1) = dx * dz * u
                            uao(al + 5, i, 1) = dy * dz * u

                            al = al + 6

                        else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! 3f GTO
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

                            !              // f_xxx, f_yyy, f_zzz
                            uao(al, i, 1) = dx2 * dx * u
                            uao(al + 1, i, 1) = dy2 * dy * u
                            uao(al + 2, i, 1) = dz2 * dz * u

                            !              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                            u = sqr5 * u                   ! correction of norm
                            uao(al + 3, i, 1) = dx2 * dy * u
                            uao(al + 4, i, 1) = dx2 * dz * u
                            uao(al + 5, i, 1) = dy2 * dx * u
                            uao(al + 6, i, 1) = dy2 * dz * u
                            uao(al + 7, i, 1) = dz2 * dx * u
                            uao(al + 8, i, 1) = dz2 * dy * u

                            !              // f_xyz
                            u = sqr3 * u                  ! correction of norm
                            uao(al + 9, i, 1) = dxyz * u

                            al = al + 10

                        else if (bl(bf)=='F'.and.gaussFOrder) then         ! 3f GTO
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

                            !              // f_xxx, f_yyy, f_zzz
                            uao(al, i, 1) = dx2 * dx * u
                            uao(al + 1, i, 1) = dy2 * dy * u
                            uao(al + 2, i, 1) = dz2 * dz * u

                            !              // f_xyy, f_xxy, f_xxz, f_xzz, f_yzz, f_yyz
                            u = sqr5 * u                   ! correction of norm
                            uao(al + 3, i, 1) = dy2 * dx * u
                            uao(al + 4, i, 1) = dx2 * dy * u
                            uao(al + 5, i, 1) = dx2 * dz * u
                            uao(al + 6, i, 1) = dz2 * dx * u
                            uao(al + 7, i, 1) = dz2 * dy * u
                            uao(al + 8, i, 1) = dy2 * dz * u

                            !              // f_xyz
                            u = sqr3 * u                  ! correction of norm
                            uao(al + 9, i, 1) = dxyz * u

                            al = al + 10

                        else if (bl(bf)=='G') then     ! 5g GTO
                            !              // do all 15 cartesian G simultaneously (same exponent is required)
                            uao(al:al + 14, i, 1) = 0d0

                            if (gaussFOrder) then
                                call internal_GaussianOrderGFunctions()
                            else
                                call internal_GamessOrderGFunctions()
                            end if

                            al = al + 15

                        else
                            call abortp('(getaos): wrong GTO')
                        end if  ! bl

                    else     !more then 1 GTO in contraction, splines used !
                        r2 = rr * rr
                        j = (csplnpnt - 1) * rr / (csalpha + rr) + 1
                        df = rr - csplx(j)

                        !           // only primitive cartesian gaussians: 1s,2p,3d,4f
                        !           // i.e. no r factor. Thus nn is not used here.

                        if (bl(bf) == 'S') then                 ! 1s GTO
                            ispl = 3 * so(bf) - 2
                            uao(al, i, 1) = cspla(ispl, j) + df * (csplb(ispl, j)&
                                    + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                            al = al + 1

                        else if (bl(bf) == 'P') then             ! 2p GTO's
                            !              // do all 3 P simultaneously (same exponent is required)
                            !              // order p_x,p_y,p_z
                            ispl = 3 * so(bf) - 2
                            u = cspla(ispl, j) + df * (csplb(ispl, j)&
                                    + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                            dx = xx - atoms(a)%cx
                            dy = yy - atoms(a)%cy
                            dz = zz - atoms(a)%cz
                            uao(al, i, 1) = dx * u
                            uao(al + 1, i, 1) = dy * u
                            uao(al + 2, i, 1) = dz * u

                            al = al + 3

                        else if (bl(bf) == 'D') then         ! 3d GTO
                            !              // do all 6 D simultaneously (same exponent is required)
                            !              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
                            ispl = 3 * so(bf) - 2
                            u = cspla(ispl, j) + df * (csplb(ispl, j)&
                                    + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                            dx = xx - atoms(a)%cx
                            dx2 = dx * dx
                            dy = yy - atoms(a)%cy
                            dy2 = dy * dy
                            dz = zz - atoms(a)%cz
                            dz2 = dz * dz
                            uao(al, i, 1) = dx2 * u
                            uao(al + 1, i, 1) = dy2 * u
                            uao(al + 2, i, 1) = dz2 * u

                            u = sqr3 * u                   ! correction of norm for last 3

                            uao(al + 3, i, 1) = dx * dy * u
                            uao(al + 4, i, 1) = dx * dz * u
                            uao(al + 5, i, 1) = dy * dz * u

                            al = al + 6

                        else if (bl(bf)=='F'.and..not.gaussFOrder) then  ! 3f GTO
                            !              // do all 10 F simultaneously (same exponent is required)
                            !              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
                            !              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

                            ispl = 3 * so(bf) - 2
                            u = cspla(ispl, j) + df * (csplb(ispl, j)&
                                    + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                            dx = xx - atoms(a)%cx
                            dx2 = dx * dx
                            dy = yy - atoms(a)%cy
                            dy2 = dy * dy
                            dz = zz - atoms(a)%cz
                            dz2 = dz * dz
                            dxyz = dx * dy * dz

                            !              // f_xxx, f_yyy, f_zzz
                            uao(al, i, 1) = dx2 * dx * u
                            uao(al + 1, i, 1) = dy2 * dy * u
                            uao(al + 2, i, 1) = dz2 * dz * u

                            !              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                            u = sqr5 * u                   ! correction of norm

                            uao(al + 3, i, 1) = dx2 * dy * u
                            uao(al + 4, i, 1) = dx2 * dz * u
                            uao(al + 5, i, 1) = dy2 * dx * u
                            uao(al + 6, i, 1) = dy2 * dz * u
                            uao(al + 7, i, 1) = dz2 * dx * u
                            uao(al + 8, i, 1) = dz2 * dy * u

                            !              // f_xyz
                            u = sqr3 * u                  ! correction of norm

                            uao(al + 9, i, 1) = dxyz * u

                            al = al + 10

                        else if (bl(bf)=='F'.and.gaussFOrder) then  ! 3f GTO
                            !              // do all 10 F simultaneously (same exponent is required)
                            !              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
                            !              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

                            ispl = 3 * so(bf) - 2
                            u = cspla(ispl, j) + df * (csplb(ispl, j)&
                                    + df * (csplc(ispl, j) + df * cspld(ispl, j)))

                            dx = xx - atoms(a)%cx
                            dx2 = dx * dx
                            dy = yy - atoms(a)%cy
                            dy2 = dy * dy
                            dz = zz - atoms(a)%cz
                            dz2 = dz * dz
                            dxyz = dx * dy * dz

                            !              // f_xxx, f_yyy, f_zzz
                            uao(al, i, 1) = dx2 * dx * u
                            uao(al + 1, i, 1) = dy2 * dy * u
                            uao(al + 2, i, 1) = dz2 * dz * u

                            !              // f_xyy, f_xxy, f_xxz, f_xzz, f_yzz, f_yyz
                            u = sqr5 * u                   ! correction of norm

                            uao(al + 3, i, 1) = dy2 * dx * u
                            uao(al + 4, i, 1) = dx2 * dy * u
                            uao(al + 5, i, 1) = dx2 * dz * u
                            uao(al + 6, i, 1) = dz2 * dx * u
                            uao(al + 7, i, 1) = dz2 * dy * u
                            uao(al + 8, i, 1) = dy2 * dz * u

                            !              // f_xyz
                            u = sqr3 * u                  ! correction of norm

                            uao(al + 9, i, 1) = dxyz * u

                            al = al + 10

                        else
                            call abortp('(getaos): wrong GTO')
                        end if ! bl
                    end if    ! simple GTO, no splines
                end if  ! STO/GTO
            end do  ! bf-loop over basis functions
        end do  ! i-loop over electrons

    contains

        subroutine internal_GamessOrderGFunctions()

            alp = cntrctn(1, 1, bf)
            u = cntrctn(2, 1, bf) * exp(-alp * r2)
            dx = xx - atoms(a)%cx
            dx2 = dx * dx
            dy = yy - atoms(a)%cy
            dy2 = dy * dy
            dz = zz - atoms(a)%cz
            dz2 = dz * dz
            dxyz = dx * dy * dz

            !           // g_xxxx, g_yyyy, g_zzzz
            uao(al, i, 1) = uao(al, i, 1) + dx2 * dx2 * u
            uao(al + 1, i, 1) = uao(al + 1, i, 1) + dy2 * dy2 * u
            uao(al + 2, i, 1) = uao(al + 2, i, 1) + dz2 * dz2 * u

            !           // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
            u = sqr7 * u                   ! correction of norm

            uao(al + 3, i, 1) = uao(al + 3, i, 1) + dx2 * dx * dy * u
            uao(al + 4, i, 1) = uao(al + 4, i, 1) + dx2 * dx * dz * u
            uao(al + 5, i, 1) = uao(al + 5, i, 1) + dy2 * dy * dx * u
            uao(al + 6, i, 1) = uao(al + 6, i, 1) + dy2 * dy * dz * u
            uao(al + 7, i, 1) = uao(al + 7, i, 1) + dz2 * dz * dx * u
            uao(al + 8, i, 1) = uao(al + 8, i, 1) + dz2 * dz * dy * u

            !           // g_xxyy, g_xxzz, g_yyzz
            u = sqr5 / sqr3 * u          ! correction of norm

            uao(al + 9, i, 1) = uao(al + 9, i, 1) + dx2 * dy2 * u
            uao(al + 10, i, 1) = uao(al + 10, i, 1) + dx2 * dz2 * u
            uao(al + 11, i, 1) = uao(al + 11, i, 1) + dy2 * dz2 * u

            !           // g_xxyz, g_yyxz, g_zzxy
            u = sqr3 * u                  ! correction of norm

            uao(al + 12, i, 1) = uao(al + 12, i, 1) + dx * dxyz * u
            uao(al + 13, i, 1) = uao(al + 13, i, 1) + dy * dxyz * u
            uao(al + 14, i, 1) = uao(al + 14, i, 1) + dz * dxyz * u

        end subroutine internal_GamessOrderGFunctions


        subroutine internal_GaussianOrderGFunctions()

            alp = cntrctn(1, 1, bf)
            u = cntrctn(2, 1, bf) * exp(-alp * r2)
            dx = xx - atoms(a)%cx
            dx2 = dx * dx
            dy = yy - atoms(a)%cy
            dy2 = dy * dy
            dz = zz - atoms(a)%cz
            dz2 = dz * dz
            dxyz = dx * dy * dz

            !           // g_xxxx, g_yyyy, g_zzzz
            uao(al + 14, i, 1) = uao(al + 14, i, 1) + dx2 * dx2 * u
            uao(al + 4, i, 1) = uao(al + 4, i, 1) + dy2 * dy2 * u
            uao(al, i, 1) = uao(al, i, 1) + dz2 * dz2 * u

            !           // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
            u = sqr7 * u                   ! correction of norm

            uao(al + 13, i, 1) = uao(al + 13, i, 1) + dx2 * dx * dy * u
            uao(al + 12, i, 1) = uao(al + 12, i, 1) + dx2 * dx * dz * u
            uao(al + 8, i, 1) = uao(al + 8, i, 1) + dy2 * dy * dx * u
            uao(al + 3, i, 1) = uao(al + 3, i, 1) + dy2 * dy * dz * u
            uao(al + 5, i, 1) = uao(al + 5, i, 1) + dz2 * dz * dx * u
            uao(al + 1, i, 1) = uao(al + 1, i, 1) + dz2 * dz * dy * u

            !           // g_xxyy, g_xxzz, g_yyzz
            u = sqr5 / sqr3 * u          ! correction of norm

            uao(al + 11, i, 1) = uao(al + 11, i, 1) + dx2 * dy2 * u
            uao(al + 9, i, 1) = uao(al + 9, i, 1) + dx2 * dz2 * u
            uao(al + 2, i, 1) = uao(al + 2, i, 1) + dy2 * dz2 * u

            !           // g_xxyz, g_yyxz, g_zzxy
            u = sqr3 * u                  ! correction of norm

            uao(al + 10, i, 1) = uao(al + 10, i, 1) + dx * dxyz * u
            uao(al + 7, i, 1) = uao(al + 7, i, 1) + dy * dxyz * u
            uao(al + 6, i, 1) = uao(al + 6, i, 1) + dz * dxyz * u

        end subroutine internal_GaussianOrderGFunctions


    end subroutine ao1splcalc


end module aos_m
