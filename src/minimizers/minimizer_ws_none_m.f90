! Copyright (C) 2020 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_ws_none_module

    ! dummy minimizer doing nothing
    use kinds_m, only: r8
    use fctn_module, only: Function_t
    use minimizer_w_sing_module, only: minimizer_w_sing
    use, intrinsic:: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    implicit none

    private
    public :: minimizer_ws_none


    type, extends(minimizer_w_sing) :: minimizer_ws_none
    contains
        procedure :: minimize => minimizer_ws_none_minimize
        procedure :: write_params_minimizer_ws_none
        procedure :: write_params => write_params_minimizer_ws_none
    end type minimizer_ws_none

    interface minimizer_ws_none
        module procedure constructor1
    end interface minimizer_ws_none

contains

    function constructor1()
        type(minimizer_ws_none), pointer :: constructor1
        allocate(constructor1)
    end function constructor1

    subroutine minimizer_ws_none_minimize(this, fn, x)
        class(minimizer_ws_none), intent(inout) :: this
        class(Function_t), intent(in)                 :: fn
        real(r8), intent(inout)                       :: x(:)

        real(r8)                                      :: f
        integer                                       :: verbose, iul

        verbose = this%verbose()
        iul = this%verbose_unit()

        if (verbose > 0) write(iul,"(a)") " *** NONE minimizer ***"
        call this%reset()

        f = fn%eval(x)

        if (verbose > 0) then
            write(iul,"(a,g20.10)") " initial position with function value f=", f
        end if

        x = ieee_value(0._r8,ieee_quiet_nan)
        call this%set_converged(.true.)
        call this%set_value(ieee_value(0._r8,ieee_quiet_nan))

        if (verbose > 0) write(iul,"(a)") " *** end NONE minimizer ***"
    end subroutine minimizer_ws_none_minimize

    subroutine write_params_minimizer_ws_none(this, iu)
        class(minimizer_ws_none), intent(in) :: this
        integer, intent(in), optional :: iu
        integer iull

        iull = this%verbose_unit()
        if (present(iu)) iull = iu
        write(iull,*) " * * * No minimization! This is a dummy =) * * *"
    end subroutine write_params_minimizer_ws_none

end module minimizer_ws_none_module
