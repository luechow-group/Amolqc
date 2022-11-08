! Copyright (C) 2017-2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module fGradPsi2_m
    use, intrinsic :: IEEE_ARITHMETIC, only: IEEE_IS_NAN
    use kinds_m, only: r8
    use error_m
    use elocData_m, only: elPhi, elU, elxDrift, elyDrift, elzDrift
    use eloc_m, only: eloc, ELOC_X_INF_ERR
    use eConfigs_m
    use fctn_module
    implicit none

    type, extends(Function_t) :: fctn_gradpsi2
    contains
        procedure :: eval => fgradpsi2_f
    end type fctn_gradpsi2

contains

    function fgradpsi2_f(this, x) result(f)
        class(fctn_gradpsi2), intent(in)  :: this
        real(r8), intent(in)  :: x(:)      ! electron vector (without elecs at cores)
        real(r8)              :: f         ! function value
        real(r8), allocatable :: g(:)      ! gradient
        type(eConfigArray)    :: eca
        real(r8), allocatable :: xx(:), yy(:), zz(:)
        integer n, i, error_code

        if (asserts) then
            call assert(mod(size(x),3) == 0, "fgradpsi2_f: illegal size")
            call assert(all(abs(x)<huge(1.d0)), "fgradpsi2_f: illegal x coords")
        end if

        allocate(g(SIZE(x)))

        n = size(x)/3
        allocate(xx(n), yy(n), zz(n))
        do i = 1, n
            xx(i) = x(3*i-2)
            yy(i) = x(3*i-1)
            zz(i) = x(3*i)
        end do

        call eConfigArray_new(eca, n, 1)
        call eConfigArray_set(eca, 1, xx, yy, zz)
        call eloc(0, eca, 'none', error_code)

        if (debug) then
            if (error_code == ELOC_X_INF_ERR) then
                call abortp("fpsi2_fg: eloc returns ELOC_X_INF_ERR")
            end if
        end if

        do i = 1, n
            g(3*i-2) = - 2.d0 * elxDrift(i, 1)
            g(3*i-1) = - 2.d0 * elyDrift(i, 1)
            g(3*i)   = - 2.d0 * elzDrift(i, 1)
        end do

        ! sing correction
        where (IEEE_IS_NAN(g)) g = 0._r8

        f = NORM2(g)

        call eConfigArray_destroy(eca)
    end function fgradpsi2_f

end module fGradPsi2_m
