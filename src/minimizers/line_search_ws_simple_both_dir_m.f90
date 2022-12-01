! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module line_search_ws_simple_both_dir_m
    use kinds_m, only: r8
    use error_m, only: assert, asserts, debug, abortp
    use fctn_module, only: Function_t
    use line_search_ws_m, only: line_search_ws
    use singularityCorrection_m, only: singularity_correction
    use singularityParticles_m, only: singularity_particles

    implicit none

    private

    type, extends(line_search_ws) :: line_search_ws_simple_both_dir
    contains
        procedure   :: find_step => find_step_ws_simple_both_dir
    end type line_search_ws_simple_both_dir

    public :: line_search_ws_simple_both_dir

    ! default constructor

contains

    subroutine find_step_ws_simple_both_dir(this, ff, sc, x0, falpha, grad0, d, xalpha, gradalpha, nfeval, is_corrected&
            &, sp)
        class(line_search_ws_simple_both_dir), intent(in) :: this
        class(Function_t), intent(in)             :: ff  ! function object
        type(singularity_correction), intent(inout) :: sc     ! out only for debugging: increase verbosity
        real(r8), intent(in)    :: x0(:)   ! position
        real(r8), intent(inout) :: falpha    ! function value (in: at x, out: at xalpha)
        real(r8), intent(in)    :: grad0(:)   ! gradient
        real(r8), intent(in)    :: d(:)   ! search direction
        real(r8), intent(inout) :: xalpha(:)   ! new position
        real(r8), intent(inout) :: gradalpha(:)   ! gradient at new position
        integer, intent(inout):: nfeval    ! in: maxiter, out: function evaluations
        logical, intent(inout):: is_corrected ! indicates whether step is scaled direction or a corrected scaled direction
        !logical, intent(out) :: is_backwards
        type(singularity_particles), intent(inout) :: sp
        real(r8) :: p(size(x0)), xOld(size(x0)), gOld(size(x0))
        real(r8) alpha, f, dir_deriv, fOld, thresh, rho, new_alpha
        integer i, maxiter
        type(singularity_particles) :: sp_old

        if (asserts) then
            call assert(size(x0) == size(grad0), "find_step_ws_simple: illegal size")
            call assert(all(abs(x0)<huge(1.d0)), "find_step_ws_simple: illegal x coords")
        end if

        f = falpha     ! initial function value
        alpha = 1._r8
        rho = 2._r8/9._r8

        ! restriction of step length
        if (this%delta_max_ > 0._r8) then
            do i = 1, SIZE(d), 3
                new_alpha = MIN(alpha, this%delta_max_ / NORM2(d(i:i+2)) )
            end do
            rho = rho * new_alpha/alpha
            alpha = new_alpha
        end if

        maxiter = 9
        thresh = 1e-16_r8

        xalpha = x0
        p = alpha * d

        if (this%step_max_ > 0.d0) call this%restrict_particle_step(p, this%step_max_)

        sp_old = sp

        call sc%correct_for_singularities(xalpha, p, sp, is_corrected, correction_only=.false.)

        fOld = ff%eval(xalpha)
        xOld = xalpha

        do nfeval = 1, maxiter

            xalpha = x0
            alpha = alpha - rho
            sp = sp_old

            p = alpha * d

            if (this%step_max_ > 0.d0) call this%restrict_particle_step(p, this%step_max_)

            if (debug) then
                if (.not.all(abs(xalpha)<huge(1.d0))) then
                    call abortp("find_step_ws_simple: illegal x coords before correct_for_singularities")
                end if
                if (.not.all(abs(p)<huge(1.d0))) then
                    call abortp("find_step_ws_simple: illegal p coords before correct_for_singularities")
                end if
            end if

            call sc%correct_for_singularities(xalpha, p, sp, is_corrected, correction_only=.false.)

            if (debug) call internal_test_xalpha()

            falpha =  ff%eval(xalpha)
            print*, falpha, alpha, rho, dot_product(grad0, p)

            if (f - fOld > thresh .and. falpha - fOld > thresh) then ! f - fOld > thresh .and.
                call ff%eval_fg(xalpha, falpha, gradalpha)
                falpha = fOld
                xalpha = xOld
                gradalpha = gOld
                exit
            else if (nfeval == maxiter) then
                call ff%eval_fg(xalpha, falpha, gradalpha)
                falpha = fOld
                xalpha = xOld
                gradalpha = gOld
                exit
            else
                fOld = falpha
                xOld = xalpha
            end if

        end do

        !!!print*, "DBG:LSS2: ", nfeval, alpha

    contains

        subroutine internal_test_xalpha()
            if (.not.all(abs(xalpha)<huge(1.d0))) then  ! repeat last call
                xalpha = x0
                sp = sp_old
                call sc%set_verbose(5)
                call sc%correct_for_singularities(xalpha, p, sp, is_corrected, correction_only=.false.)
                call abortp("find_step_ws_simple: illegal x coords after correct_for_singularities")
            end if
        end subroutine internal_test_xalpha

    end subroutine find_step_ws_simple_both_dir


end module line_search_ws_simple_both_dir_m
