! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later


module nl2sol_i
    use kinds_m, only: r8, i4

    implicit none

    ! nl2sol.f90
    interface
        subroutine nl2sno ( n, p, x, calcr, iv, v, uiparm, urparm, ufparm )
            import :: r8, i4
            integer(i4) :: n
            integer(i4) :: p
            real(r8) :: x(p)
            interface                     ! external function
                subroutine calcr ( meqn, nvar, p, nf, r, uiparm, urparm, ufparm )
                    import :: r8, i4
                    integer(i4) :: meqn
                    integer(i4) :: nvar
                    real(r8) :: p(nvar)
                    integer(i4) :: nf
                    real(r8) :: r(meqn)
                    integer(i4) ::  uiparm(*)
                    real(r8) :: urparm(2 * meqn)
                    external       :: ufparm
                end subroutine calcr
            end interface
            integer(i4) :: iv(60+p)
            real(r8) :: v(93 + n*p + 3*n + (p*(3*p+33))/2)
            integer(i4) :: uiparm(*)
            real(r8) :: urparm(*)
            interface                     ! external function
                subroutine ufparm ( meqn, nvar, x )
                    import :: r8, i4
                    integer(i4) :: meqn
                    integer(i4) :: nvar
                    real(r8) :: x(nvar)
                end subroutine ufparm
            end interface
        end subroutine nl2sno
    end interface

    interface
        subroutine dfault ( iv, v )
            import :: r8, i4
            integer(i4) :: iv(25)
            real(r8) :: v(45)
        end subroutine dfault
    end interface

    interface
        subroutine itsmry ( d, iv, p, v, x )
            import :: r8, i4
            integer(i4) :: p ! placed here because of d(p)
            real(r8) :: d(p)
            integer(i4) :: iv(*)
            real(r8) :: v(*)
            real(r8) :: x(p)
        end subroutine itsmry
    end interface

    interface
        subroutine nl2itr ( d, iv, j, n, nn, p, r, v, x )
            import :: r8, i4
            integer(i4) :: p ! placed here because of d(p)
            real(r8) :: d(p)
            integer(i4) :: iv(60+p)
            integer(i4) :: nn ! placed here because of j(nn,p)
            real(r8) :: j(nn,p)
            integer(i4) :: n
            real(r8) :: r(n)
            real(r8) :: v(93 + 2*n + (p*(3*p+31))/2)
            real(r8) :: x(p)
        end subroutine nl2itr
    end interface

end module nl2sol_i
