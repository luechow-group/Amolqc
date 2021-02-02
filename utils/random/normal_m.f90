! Copyright (C) 2020 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

! this module provides standard (normal) distributed numbers
module normal_m
    use kinds_m, only: r8
    use random_m, only: myran, myran_is_initialized
    use error_m, only: assert

    implicit none
    private
    public :: normal_ran, normal, box_muller

    real(r8), parameter :: PI = 2 * ATAN(1._r8)
    logical :: first_m = .True.
    real(r8) :: u1_m, u2_m

contains
    ! using the Box-Muller transform
    function normal_ran() result(random_number)
        real(r8) :: random_number

        call assert(myran_is_initialized(), 'normal_ran: myran has to be initialized')

        if (first_m) then
            u1_m = myran()
            u2_m = myran()
        end if

        random_number = box_muller(u1_m, u2_m, first_m)

        first_m = .not. first_m
    end function normal_ran

    function box_muller(u1, u2, first) result(random_number)
        real(r8), intent(in) :: u1, u2
        logical, intent(in) :: first
        real(r8) :: random_number

        real(r8) :: theta

        random_number = SQRT(-2*LOG(u1))

        theta = 2 * PI * u2

        if (first) then
            random_number = random_number * SIN(theta)
        else
            random_number = random_number * COS(theta)
        end if
    end function box_muller

    ! this is just the standard normal distribution
    function normal(x) result(f)
        real(r8), intent(in) :: x
        real(r8) :: f

        f = 1._r8 / SQRT(2*PI) * EXP(-x**2/2._r8)
    end function normal

end module normal_m