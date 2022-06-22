! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_steep_desc_module

   use kinds_m, only: r8
   use fctn_module, only: Function_t
   use minimizer_module, only: minimizer
   use line_search_simple_module, only: line_search_simple
   implicit none

   private
   public :: minimizer_steep_desc

   type, extends(minimizer) :: minimizer_steep_desc
      type(line_search_simple)  :: lss_
      real(r8)                    :: dt_ = 1.d0         ! step size in front of gradient
   contains
      procedure    :: minimize => minimizer_steep_desc_minimize
      procedure    :: write_params => write_params_minimizer_steep_desc
   end type minimizer_steep_desc

   interface minimizer_steep_desc
      module procedure constructor
   end interface minimizer_steep_desc

contains

   function constructor(lss, dt)
      type(line_search_simple)  :: lss
      real(r8), intent(in) :: dt
      type(minimizer_steep_desc), pointer :: constructor
      allocate(constructor)
      constructor%lss_ = lss
      constructor%dt_ = dt
   end function constructor

   subroutine minimizer_steep_desc_minimize(this, fn, x)
      class(minimizer_steep_desc), intent(inout) :: this
      class(Function_t), intent(in)               :: fn
      real(r8), intent(inout)                      :: x(:)
      real(r8)                     :: f
      real(r8)                     :: g(size(x))
      real(r8)                     :: delta_x(size(x)), x_new(size(x)), g_new(size(x))
      real(r8)                     :: gmax
      integer                    :: iter, verbose, iul, n_eval, ls_iter

      verbose = this%verbose()
      iul = this%verbose_unit()
      if (verbose > 0) write(iul,"(a)") " *** steepest descent minimizer ***"
      call this%reset()

      call fn%eval_fg(x, f, g)

      if (verbose > 0) then 
         write(iul,"(a,g20.10)") " initial position with function value f=", f
         call this%write_params(iul)
         if (verbose > 1) then
            write(iul,"(9g13.5)") x
         end if
      end if

      n_eval = 1

      do iter = 1, this%max_iterations()

         ! steepest descent step
         delta_x = - this%dt_ * g

         ls_iter = 5
         call this%lss_%find_step(fn, x, f, g, delta_x, x_new, g_new, ls_iter)
         n_eval = n_eval + ls_iter

         x = x_new
         g = g_new

         call fn%eval_fg(x, f, g)
         gmax = maxval(abs(g))

         if (this%verbose() > 2) then
            write(iul,"(i6,3(a,g18.10))") iter, " f=", f, " gmax=", gmax
            if (verbose > 3) then
               write(iul,"(9g13.5)") x
               write(iul,"(9g13.5)") g
            end if
         end if

         if (this%is_gradient_converged(gmax) .or. this%is_value_converged(f)) then
            call this%set_converged(.true.)
            exit
         end if

      end do

      call this%set_value(f)
      call this%set_gradient(g)
      call this%set_iterations(iter)
      call this%set_function_evaluations(n_eval)
   end subroutine minimizer_steep_desc_minimize

   subroutine write_params_minimizer_steep_desc(this, iu)
      class(minimizer_steep_desc), intent(in) :: this
      integer, intent(in), optional :: iu
      integer iull

      if (present(iu)) then
         call this%write_params_minimizer(iu) ! call parent class
      else
         call this%write_params_minimizer()
      end if
      if (present(iu)) then
         iull = iu
      else
         iull = this%verbose_unit()
      end if
      write(iull,"(a,g13.5)")  " step_size           =", this%dt_
   end subroutine write_params_minimizer_steep_desc


end module minimizer_steep_desc_module
