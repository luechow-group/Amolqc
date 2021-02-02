! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module boltzmann_m
    use kinds_m, only: i4, r8
    use nl2sol_i, only: nl2sno

    implicit none

    private
    public :: Boltzmann_fit

    contains
        function Boltzmann_fit(x, y, size) result(p)
            integer(i4)              :: size
            real(r8),    intent(in) :: x(size), y(size)
            integer (i4), parameter  :: nvar = 4

            real(r8)                :: p(nvar)

            integer(i4)              :: iv(60+nvar), uiparm(1), meqn
            real(r8)                :: urparm(2 * size)
            real(r8)                :: v(93+size*nvar+3*size+nvar*(3*nvar+33)/2)

            meqn = size
            urparm(1 : size) = x
            urparm(size + 1 : 2 * size) = y
            p = 1d0
            iv(1) = 0

            call nl2sno ( meqn, nvar, p, Boltzmann_residuals, iv, v, uiparm, urparm, ufparm )

        end function Boltzmann_fit


        subroutine Boltzmann_residuals ( meqn, nvar, p, nf, r, uiparm, urparm, ufparm )
            integer(i4) :: meqn, nvar

            integer(i4) :: nf, uiparm(*), i
            external       :: ufparm
            real(r8)   :: urparm(2 * meqn), p(nvar), r(meqn), x(meqn), y(meqn)

            x = urparm(1 : meqn)
            y = urparm(meqn + 1 : 2 * meqn)

            do i = 1, meqn
                r(i) = Boltzmann_value(p, x(i)) - y(i)
            end do

        end subroutine Boltzmann_residuals

        pure function Boltzmann_value(p, x) result(f)
            real(r8), intent(in) :: p(4), x
            real(r8)             :: f

            f = p(2) + (p(1) - p(2))/(1 + exp((x - p(3)) / p(4) ) )

        end function Boltzmann_value

        subroutine ufparm ( meqn, nvar, x )
            integer(i4) :: meqn
            integer(i4) :: nvar

            real(r8)   :: x(nvar)

        end subroutine ufparm

end module boltzmann_m
