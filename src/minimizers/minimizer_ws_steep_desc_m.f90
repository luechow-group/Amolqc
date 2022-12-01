! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_ws_steep_desc_module

   ! steepest descent using particle singularity information
   use kinds_m, only: r8
   use error_m, only: assert, asserts
   use fctn_module, only: Function_t
   use singularityCorrection_m, only: singularity_correction
   use singularityParticles_m, only: singularity_particles
   use line_search_ws_m, only: line_search_ws
   use minimizer_w_sing_module, only: minimizer_w_sing
   implicit none

   private
   public :: minimizer_ws_steep_desc


   type, extends(minimizer_w_sing) :: minimizer_ws_steep_desc
      class(line_search_ws), allocatable  :: lss_
      real(r8)                 :: dt_ = 1.d0         ! gradient prefactor
      integer                      :: uphill_steps_ = 0
   contains
      procedure :: minimize => minimizer_ws_steep_desc_minimize
      procedure :: reset_uphill_steps => minimizer_ws_steep_desc_reset_uphill_steps
      procedure :: uphill_steps => minimizer_ws_steep_desc_uphill_steps
      procedure :: write_params_minimizer_ws_steep_desc
      procedure :: write_params => write_params_minimizer_ws_steep_desc
   end type minimizer_ws_steep_desc 

   interface minimizer_ws_steep_desc
      module procedure constructor1
   end interface minimizer_ws_steep_desc

contains

   function constructor1(lss, sc, dt)
      class(line_search_ws), intent(in)        :: lss
      type(singularity_correction), intent(in)        :: sc
      real(r8), intent(in)                              :: dt
      type(minimizer_ws_steep_desc), pointer :: constructor1
      allocate(constructor1)
      constructor1%lss_ = lss
      constructor1%sc_ = sc
      constructor1%dt_ = dt
   end function constructor1

   subroutine minimizer_ws_steep_desc_reset_uphill_steps(this)
      class(minimizer_ws_steep_desc), intent(inout) :: this
      this%uphill_steps_ = 0
   end subroutine minimizer_ws_steep_desc_reset_uphill_steps

   function minimizer_ws_steep_desc_uphill_steps(this) result(res)
      class(minimizer_ws_steep_desc), intent(inout) :: this
      integer                                       :: res
      res = this%uphill_steps_
   end function minimizer_ws_steep_desc_uphill_steps

   subroutine minimizer_ws_steep_desc_minimize(this, fn, x)
      class(minimizer_ws_steep_desc), intent(inout) :: this
      class(Function_t), intent(in)                  :: fn
      real(r8), intent(inout)                         :: x(:)
      real(r8)                     :: f, f_old
      real(r8)                     :: g(size(x)), g_new(size(x)), x_new(size(x))
      real(r8)                     :: delta_x(size(x))
      real(r8)                     :: gmax, dt
      integer                    :: iter, verbose, iul, ls_iter, n_eval
      type(singularity_particles):: sp
      logical                    :: is_corrected

      if (asserts) call assert(mod(size(x),3) == 0, "minimizer_ws_steep_desc_minimize: illegal size")

      call sp%create(size(x)/3, this%sc_%n_singularities())

      verbose = this%verbose()
      iul = this%verbose_unit()
      if (verbose > 0) write(iul,"(a)") " *** steepest descent minimizer with singularity information ***"
      call this%reset()

      call fn%eval_fg(x, f, g)
      n_eval = 1

      ! before the first step, check for singularities, but with zero step and set particles to singularities
      delta_x = 0._r8
      call this%restrict_gradient(g)
      call this%sc_%correct_for_singularities(x, delta_x, sp, is_corrected, correction_only=.false.)
      where ( sp%At_singularity() ) g = 0.d0

      if (verbose > 0) then 
         write(iul,"(a,g20.10)") " initial position with function value f=", f
         call this%write_params(iul)
         if (verbose > 1) then
            write(iul,"(9g13.5)") x
            write(iul,"(9g13.5)") g
         end if
      end if

      if (this%do_write_opt_path()) call this%write_opt_path_entry(1, x, f)

      do iter = 1, this%max_iterations()

         if (verbose > 3) write(iul,"(a,i6)") " iter=", iter

         dt = this%dt_
         call this%restrict_gradient(g)
         delta_x = - dt * g
         f_old = f
         ls_iter = 5
         call this%lss_%find_step(fn, this%sc_, x, f, g, delta_x, x_new, g_new, ls_iter, is_corrected, sp)

         if (verbose > 4) write(iul,"(2i6,3g18.10)") iter, ls_iter, f, f_old, dt
         n_eval = n_eval + ls_iter

         if (this%do_write_opt_path()) call this%write_opt_path_entry(iter + 1, x_new, f)

         x = x_new
         g = g_new
         call this%restrict_gradient(g)
         where ( sp%At_singularity() ) g = 0.d0

         gmax = maxval(abs(g))

         if (verbose > 2) then
            write(iul,"(2i6,2(a,g18.10),2(a,g12.3))") iter, ls_iter, " f=", f, " gmax=", gmax, " dt=", dt
            if (verbose > 3) then
               write(iul,"(9g13.5)") x
               write(iul,"(9g13.5)") g
            end if
         end if

         if (this%is_gradient_converged(gmax) .or. (sp%n_sing() == size(x)/3) .or. this%is_value_converged(f)) then
            call this%set_converged(.true.)
            exit
         end if

      end do

      if (verbose > 0) then
         write(iul,"(a,g20.10)") " final position with function value f=", f
         write(iul,"(a,l5,2(a,i8))") " converged=",this%is_converged(), " with iter=", iter, " fctn_eval=", n_eval
         if (verbose > 1) then
            write(iul,"(9g13.5)") x
            write(iul,"(9g13.5)") g
         end if
      end if

      call this%set_value(f)
      call this%set_gradient(g)
      call this%set_iterations(iter)
      call this%set_function_evaluations(n_eval)

   end subroutine minimizer_ws_steep_desc_minimize

   subroutine write_params_minimizer_ws_steep_desc(this, iu)
      class(minimizer_ws_steep_desc), intent(in) :: this
      integer, intent(in), optional :: iu
      integer iull

      iull = this%verbose_unit()
      if (present(iu)) iull = iu
      write(iull,*) " * * * Steepest Descent minimizer with singularity correction * * *"

      if (present(iu)) then
         call this%write_params_minimizer_w_sing(iu) ! call parent class
         call this%lss_%write_params(iu)
      else
         call this%write_params_minimizer_w_sing()
         call this%lss_%write_params()
      end if
      write(iull,"(a,g13.5)")  " step=", this%dt_
   end subroutine write_params_minimizer_ws_steep_desc



end module minimizer_ws_steep_desc_module
