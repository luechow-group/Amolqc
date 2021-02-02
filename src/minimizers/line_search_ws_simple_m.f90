! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module line_search_ws_simple_module
   use kinds_m, only: r8
   use error_m, only: assert, asserts, debug, abortp
   use fctn_module, only: Function_t
   use singularityCorrection_m, only: singularity_correction
   use singularityParticles_m, only: singularity_particles
   use verbosity_m, only: Verbosity_t

   implicit none

   private

   type, extends(Verbosity_t) :: line_search_ws_simple
      !!!private              ! with private the default constructor is impossible
      real(r8)      :: alpha_ = 1.d0   ! line search parameter
      real(r8)      :: c_ = 0.1d0
      real(r8)      :: rho_ = 0.33d0
      real(r8)      :: delta_max_ = 0.d0   ! maximum step length (scaling search direction)
      real(r8)      :: step_max_  = 0.d0   ! maximum step length per particle (modifying search direction)
   contains
      procedure   :: find_step => find_step_ws_simple
      procedure   :: write_params => line_search_ws_simple_write_params
   end type line_search_ws_simple

   public :: line_search_ws_simple

   ! default constructor

contains

   subroutine find_step_ws_simple(this, ff, sc, x0, falpha, grad0, d, xalpha, gradalpha, nfeval, is_corrected, sp)
      class(line_search_ws_simple), intent(in) :: this
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
      type(singularity_particles), intent(inout) :: sp
      real(r8) :: p(size(x0))
      real(r8) alpha, f, dir_deriv
      integer i, maxiter
      type(singularity_particles) :: sp_old

      if (asserts) then 
         call assert(size(x0) == size(grad0), "find_step_ws_simple: illegal size")
         call assert(all(abs(x0)<huge(1.d0)), "find_step_ws_simple: illegal x coords")
      end if

      f = falpha     ! old function value
      alpha = this%alpha_

      ! restriction of step length
      if (this%delta_max_ > 0._r8) then
         do i = 1, SIZE(d), 3
            alpha = MIN(alpha, this%delta_max_ / NORM2(d(i:i+2)) )
         end do
      end if

      maxiter = nfeval

      dir_deriv = dot_product(grad0, d)    ! directional derivative

      if (dir_deriv >= 0.d0) then
         nfeval = -1
         return
      end if

      do nfeval = 1, maxiter

         xalpha = x0

         p = alpha * d

         if (this%step_max_ > 0.d0) call restrict_particle_step(p, this%step_max_)

         sp_old = sp

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

         call ff%eval_fg(xalpha, falpha, gradalpha)

      if (falpha < f + this%c_ * alpha * dir_deriv .or. nfeval == maxiter) exit

         alpha = this%rho_ * alpha
         sp = sp_old

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

   end subroutine find_step_ws_simple

   subroutine line_search_ws_simple_write_params(this, iu)
      class(line_search_ws_simple), intent(in)    :: this
      integer, intent(in) , optional  :: iu
      integer iull

      iull = this%verbose_unit()
      if (present(iu)) iull = iu
      write(iull, "(/a/a)") " -- simple line search with step correction", " with:"
      write(iull,"(a,g12.3)") " alpha               =", this%alpha_
      write(iull,"(a,g12.3)") " c                   =", this%c_
      write(iull,"(a,g12.3)") " rho                 =", this%rho_
      if (this%delta_max_ > 0.d0)  &
         write(iull, "(a,g12.3)")  " max_distance        =", this%delta_max_
      if (this%step_max_ > 0.d0)   &
         write(iull, "(a,g12.3)")  " max_step            =", this%step_max_
   end subroutine line_search_ws_simple_write_params

   subroutine restrict_particle_step(p, step_max)
      real(r8), intent(inout)  :: p(:)
      real(r8), intent(in)     :: step_max
      real(r8) :: step_length
      integer i

      call assert(mod(size(p),3) == 0, "restrict_particle_step: illegal size")
      do i = 1, size(p), 3
         step_length = NORM2(p(i:i+2))
         if (step_length > step_max) then
            p(i:i+2) = p(i:i+2) * step_max / step_length
         end if
      end do
   end subroutine restrict_particle_step


end module line_search_ws_simple_module
