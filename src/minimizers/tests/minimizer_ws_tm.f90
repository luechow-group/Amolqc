! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_ws_tm
   use kinds_m, only: r8
   use error_m, only: assertEqualRelative, assert
   use myfctn_module, only: mysingFunction_t
   use fctn_module, only: Function_t
   use minimizer_w_sing_module, only: minimizer_w_sing
   use minimizer_ws_factory_module, only: create_ws_minimizer
   implicit none

contains
   subroutine minimizer_ws_test()


      real(r8), allocatable       :: x(:), g(:)
      real(r8)                    :: f
      integer                   :: n, m
      type(mysingFunction_t), target :: fg
      class(Function_t), pointer :: fg_p => null()
      class(minimizer_w_sing), pointer :: minws_p => null()
      !!class(minimizer), pointer :: minimizer_p => null()
      character(len=120)        :: lines(7)
      real(r8), allocatable       :: singularities(:,:)


      n = 6
      allocate(x(n), g(n))
      x = [3.d0, -2.d0, 0.1d0, 0.d0, 4.d0, 1.1d0]

      allocate(singularities(3, 1))

      fg = mysingFunction_t(a=2.0d0, b=2.0d0)
      fg_p => fg

      call fg_p%eval_fg(x, f, g)
      call assertEqualRelative(f, 73.423725770523959d0, 1.d-14, "minimizer_ws_mt: singfctn eval_fg function value test (1) failed")
      call assertEqualRelative(sum(g), 4.5290299503941833d0, 1.d-14, "minimizer_ws_mt: singfctn eval_fg grad sum test (1) failed")

      !!! testing steepest descent with cut step correction

      lines(1) = "method=steepest_descent, verbose=0"
      lines(2) = "step_size=0.1, max_iter=500, no_scaling"
      lines(3) = "singularity_threshold=1.d-3, correction_threshold=0.5d0, correction_mode=cut_one"
      lines(4) = "convergence_gradient=1.d-4"

      minws_p => create_ws_minimizer(lines)
      call assert(associated(minws_p), "minimizer_ws_mt: association (1) failed")

      singularities(:, 1) = [0.d0, 0.d0, -1.d0]
      call minws_p%set_singularities(singularities)

      m = minws_p%max_iterations()
      call assert(m == 500, "minimizer_ws_mt: max_iterations test (1) failed")

      call minws_p%minimize(fg, x)

      call assertEqualRelative(minws_p%value(), 1.228848398d-5, 1.d-6, &
         "minimizer_ws_mt: minimize value test (1) failed")
      call assert(minws_p%is_converged(), "minimizer_ws_mt: convergence test (1) failed")
      call assert(minws_p%iterations() == 13, "minimizer_ws_mt: iteration test (1) failed")
      call assert(minws_p%function_evaluations() == 14, "minimizer_ws_mt: function evaluation test (1) failed")

      deallocate(minws_p)

      !!! testing steepest descent with umrigar correction

      lines(1) = "method=steepest_descent, verbose=0"
      lines(3) = "singularity_threshold=1.d-3, correction_threshold=0.5d0, correction_mode=umr_one"
      minws_p => create_ws_minimizer(lines)
      call assert(associated(minws_p), "minimizer_ws_mt: association (2) failed")

      call minws_p%set_singularities(singularities)

      x = [3.d0, -2.d0, 0.1d0, 0.d0, 4.d0, 1.1d0]

      call minws_p%minimize(fg, x)

      call assertEqualRelative(minws_p%value(), 1.228848398d-5, 1.d-6, &
         "minimizer_ws_mt: minimize value test (2) failed")
      call assert(minws_p%is_converged(), "minimizer_ws_mt: convergence test (2) failed")
      call assert(minws_p%iterations() == 13, "minimizer_ws_mt: iteration test (2) failed")
      call assert(minws_p%function_evaluations() == 14, "minimizer_ws_mt: function evaluation test (2) failed")

      deallocate(minws_p)

      !!! testing fire with cut step correction

      lines(1) = "method=fire, verbose=0"
      lines(2) = "tau_init=0.1, max_iter=200, no_scaling"
      lines(3) = "convergence_gradient=1.d-4"
      lines(4) = "singularity_threshold=1.d-3, correction_threshold=0.5d0, correction_mode=cut_one"

      minws_p => create_ws_minimizer(lines)
      call assert(associated(minws_p), "minimizer_ws_mt: association (3) failed")

      call minws_p%set_singularities(singularities)

      x = [3.d0, -2.d0, 0.1d0, 0.d0, 4.d0, 1.1d0]

      call minws_p%minimize(fg, x)

      call assertEqualRelative(minws_p%value(), 1.228904072d-5, 1.d-6, &
         "minimizer_ws_mt: minimize value test (3) failed")
      call assert(minws_p%is_converged(), "minimizer_ws_mt: convergence test (3) failed")
      call assert(minws_p%iterations() == 43, "minimizer_ws_mt: iteration test (3) failed")
      call assert(minws_p%function_evaluations() == 44, "minimizer_ws_mt: function evaluation test (3) failed")

      deallocate(minws_p)

      !!! testing bfgs with cut step correction

      lines(1) = "method=bfgs, verbose=0"
      lines(2) = "tau_init=0.1, max_iter=200, no_scaling"
      lines(3) = "convergence_gradient=1.d-4"
      lines(4) = "singularity_threshold=1.d-3, correction_threshold=0.5d0, correction_mode=cut_one"

      minws_p => create_ws_minimizer(lines)
      call assert(associated(minws_p), "minimizer_ws_mt: association (4) failed")

      call minws_p%set_singularities(singularities)

      x = [3.d0, -2.d0, 0.1d0, 0.d0, 4.d0, 1.1d0]

      call minws_p%minimize(fg, x)

      call assertEqualRelative(minws_p%value(), 1.228842503d-5, 1.d-5, &
         "minimizer_ws_mt: minimize value test (4) failed")
      call assert(minws_p%is_converged(), "minimizer_ws_mt: convergence test (4) failed")
      call assert(minws_p%iterations() == 13, "minimizer_ws_mt: iteration test (4) failed")
      call assert(minws_p%function_evaluations() == 14, "minimizer_ws_mt: function evaluation test (4) failed")

      deallocate(minws_p)

   end subroutine minimizer_ws_test
end module minimizer_ws_tm