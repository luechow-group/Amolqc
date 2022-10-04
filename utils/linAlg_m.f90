! Copyright (C) 1996, 2006-2007 Arne Luechow
! Copyright (C) 2005 Tony C. Scott
! Copyright (C) 2011, 2013 Alexander Sturm
!
! SPDX-License-Identifier: GPL-3.0-or-later

! linAlg.f90: contains miscellaneous linear algebra subroutines
!
! depends on BLAS and on UMFPACK

module linAlg_m

use kinds_m, only : r8
use globalUtils_m, only : useLAlib
implicit none

contains
!     ---------------------------------------
subroutine invcupd(col, i, n, nmax, a1, det)
    !     ---------------------------------------

    ! invupd updates the inverse matrix a1 and the determinant det
    ! for a matrix with a new column i col.
    !
    ! Both BLAS and NON-BLAS versions in one routine
    !
    ! input:
    !     col(nmax) : new column in matrix
    !             i : column index
    !             n : actual dimension of matrix a1 and vector col
    !          nmax : definition of matrix and vector
    !  a1(nmax,nmax): old inverse matrix
    !           det : old determinant

    ! output:
    !  a1(nmax,nmax): updated inverse matrix
    !           det : updated determinant
    !TS

    integer n, nmax, i, j, k, l
    real(r8) alpha, f(nmax)
    real(r8) col(nmax), a1(nmax, nmax), det, r, tmp, one, zero
    !
    if (useLAlib) then
        one = 1.0d0
        zero = 0.0d0
        call DGEMV('N', n, n, one, a1, nmax, col, 1, zero, f, 1)
        r = f(i)
        det = r * det
        alpha = one / r

        !       // update row i of inverse
        !       a1(1,i) is address of first element of i-th column
        !       In FORTRAN, matrices are stored column-wise as arrays
        call DSCAL(n, alpha, a1(i, 1), nmax)

        !       // update all other columns
        do j = 1, n
            if (j==i) goto 202
            alpha = -f(j)
            call DAXPY(n, alpha, a1(i, 1), nmax, a1(j, 1), nmax)
            202       continue
        enddo
    else
        r = 0d0
        do k = 1, n
            r = r + col(k) * a1(i, k)
        enddo
        det = r * det
        !       // update row i
        do k = 1, n
            a1(i, k) = a1(i, k) / r
        enddo

        !       // update all other rows
        do j = 1, n
            if (j==i) goto 102
            tmp = 0d0
            do l = 1, n
                tmp = tmp + col(l) * a1(j, l)
            enddo
            do k = 1, n
                a1(j, k) = a1(j, k) - a1(i, k) * tmp
            enddo
            102       continue
        enddo
    endif

end

!     ---------------------------------------
subroutine invdetcalc(col, n, nmax, a1, det)
    !     ---------------------------------------
    ! calculates the determinant for an updated inverse matrix (following the
    ! matrix determinant lemma). note that no matrices are modified by this,
    ! only the new determinant is calculated

    integer n, nmax, i
    real(r8) col(nmax), a1(nmax), det
    real(r8) ddot
    !
    if (useLAlib) then
        det = DDOT(n, col(1:n), 1, a1(1:n), 1) * det
    else
        det = dot_product(col(1:n), a1(1:n)) * det
    endif

end

!================================================

!     ---------------------------------------
subroutine invrupd(row, i, n, nmax, a1, det)
    !     ---------------------------------------

    ! invupd updates the inverse matrix a1 and the determinant det
    ! for a matrix with a new row i row(.).
    !
    ! Both BLAS and NON-BLAS versions in one routine
    !
    ! input:
    !     row(nmax) : new row in matrix
    !             i : row index
    !             n : actual dimension of matrix a1 and vector col
    !          nmax : definition of matrix and vector
    !  a1(nmax,nmax): old inverse matrix
    !           det : old determinant

    ! output:
    !  a1(nmax,nmax): updated inverse matrix
    !           det : updated determinant
    !
    !TS
    integer n, nmax, i, j, k, l
    real(r8) alpha, f(nmax), row(nmax), a1(nmax, nmax), det, r, tmp, zero, one
    !
    if (useLAlib) then
        zero = 0.d0
        one = 1.d0
        call DGEMV('T', n, n, one, a1, nmax, row, 1, zero, f, 1)
        r = f(i)
        det = r * det
        alpha = one / r

        !       // update column i of inverse
        !       a1(1,i) is address of first element of i-th column
        !       In FORTRAN, matrices are stored column-wise as arrays
        call DSCAL(n, alpha, a1(1, i), 1)

        !       // update all other columns
        do j = 1, n
            if (j==i) goto 202
            alpha = -f(j)
            call DAXPY(n, alpha, a1(1, i), 1, a1(1, j), 1)
            202      continue
        enddo

    else
        r = 0d0
        do k = 1, n
            r = r + row(k) * a1(k, i)
        enddo
        det = r * det

        !       // update column i of inverse
        do k = 1, n
            a1(k, i) = a1(k, i) / r
        enddo

        !       // update all other columns
        do j = 1, n
            if (j==i) goto 302
            tmp = 0d0
            do l = 1, n
                tmp = tmp + row(l) * a1(l, j)
            enddo
            do k = 1, n
                a1(k, j) = a1(k, j) - a1(k, i) * tmp
            enddo
            302       continue
        enddo
    endif

end


!================================================

!     --------------------------------
subroutine inv1(a, n, nmax, a1, det)
    !     --------------------------------

    ! calculates the inverse matrix of a by updating the unit matrix

    ! input:
    !     col(nmax) : new column in matrix
    !             i : column index
    !             n : actual dimension of matrix a1 and vector col
    !          nmax : definition of matrix and vector
    !  a1(nmax,nmax): old inverse matrix
    !           det : old determinant

    ! output:
    !  a1(nmax,nmax): updated inverse matrix
    !           det : updated determinant

    integer n, nmax, i, j, k, l
    real(r8) a(nmax, nmax), a1(nmax, nmax), det, r, tmp

    !     // unit matrix is self-inverse with det=1
    do j = 1, n
        do i = 1, n
            a1(i, j) = 0d0
        enddo
        a1(j, j) = 1d0
    enddo
    det = 1d0

    !     // loop over columns i for updating unit matrix
    do i = 1, n

        r = 0d0
        do k = 1, n
            r = r + a(k, i) * a1(i, k)
        enddo
        det = r * det

        !        // update row i
        do k = 1, n
            a1(i, k) = a1(i, k) / r
        enddo

        !        // update all other rows
        do j = 1, n
            if (j==i) goto 102
            tmp = 0d0
            do l = 1, n
                tmp = tmp + a(l, i) * a1(j, l)
            enddo
            do k = 1, n
                a1(j, k) = a1(j, k) - a1(i, k) * tmp
            enddo
            102        continue
        enddo

    enddo

end

end module linAlg_m
