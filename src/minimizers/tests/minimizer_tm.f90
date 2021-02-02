! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_tm
   use kinds_m, only: r8
   use error_m, only: assert, assertEqualAbsolute, assertEqualRelative
   use myfctn_module, only: myFunction_t
   use fctn_module, only: Function_t
   use minimizer_module, only: minimizer
   use minimizer_factory_module, only: create_minimizer
   implicit none

contains
   subroutine minimizer_test()

      real(r8), allocatable       :: x(:), g(:)
      real(r8)                    :: f
      integer                   :: n, m
      type(myFunction_t), target :: fg
      class(Function_t), pointer :: fg_p => null()
      class(minimizer), pointer :: minimizer_p => null()
      character(len=120)        :: lines(4) = ["", "", "", ""]

      n = 3
      allocate(x(n), g(n))
      x = [1, 2, 3]

      fg = myFunction_t(a=2.0d0)
      fg_p => fg

      call fg_p%eval_fg(x, f, g)
      call assertEqualRelative(f, 196.d0, 1.d-14, "minimizer_mt: function test (1) failed")

      g = g - [8.d0, 64.d0, 216.d0]
      f = dot_product(g, g)
      call assertEqualAbsolute(f, 0.d0, 1.d-13, "minimizer_mt: grad test (1) failed")

      !!! testing steepest descent

      lines(1) = "method=steepest_descent, verbose=0"
      lines(2) = "step_size=1.0, max_distance=0.1, max_iter=500"
      lines(3) = "convergence_gradient=1.d-4"
      minimizer_p => create_minimizer(lines)
      call assert(associated(minimizer_p), "minimizer_mt: association (1) failed")

      m = minimizer_p%max_iterations()
      call assert(m == 500, "minimizer_mt: max_iterations test (1) failed")

      call minimizer_p%minimize(fg, x)

      call assertEqualRelative(minimizer_p%value(), 1.724078157d-6, 1.d-4, &
         "minimizer_mt: minimize value test (1) failed")

      call assert(minimizer_p%is_converged(), "minimizer_mt: convergence test (1) failed")
      call assert(minimizer_p%iterations() == 147, "minimizer_mt: iterations (1) failed")
      call assert(minimizer_p%function_evaluations() == 148, "minimizer_mt: function_evaluations (1) failed")

      deallocate(minimizer_p)   !!! does this really free the memory of minimizer_p??? yes, also calls final if it exists

      !!! testing FIRE

      x = [1, 2, 3]

      lines(1) = "method=fire, verbose=0"
      lines(2) = "tau_init=0.1, max_iter=500"
      lines(3) = "convergence_gradient=1.d-4"
      minimizer_p => create_minimizer(lines)
      call assert(associated(minimizer_p), "minimizer_mt: association (2) failed")

      call minimizer_p%minimize(fg, x)

      call assertEqualRelative(minimizer_p%value(), 1.361773583693D-6, 1.d-4, &
         "minimizer_mt: minimize value test (2) failed")

      call assert(minimizer_p%is_converged(), "minimizer_mt: convergence test (2) failed")
      call assert(minimizer_p%iterations() == 69, "minimizer_mt: iterations (2) failed")
      call assert(minimizer_p%function_evaluations() == 70, "minimizer_mt: function_evaluations (2) failed")

      deallocate(minimizer_p)

      !!! testing BFGS

      x = [1, 2, 3]

      lines(1) = "method=bfgs, verbose=0"
      lines(2) = "max_iter=100"
      lines(3) = "convergence_gradient=1.d-4"
      minimizer_p => create_minimizer(lines)
      call assert(associated(minimizer_p), "minimizer_mt: association (3) failed")

      call minimizer_p%minimize(fg, x)

      call assertEqualRelative(minimizer_p%value(), 6.463772698576d-7, 1.d-4, &
         "minimizer_mt: minimize value test (3) failed")

      call assert(minimizer_p%is_converged(), "minimizer_mt: convergence test (3) failed")
      call assert(minimizer_p%iterations() == 24, "minimizer_mt: iterations (3) failed")
      call assert(minimizer_p%function_evaluations() == 43, "minimizer_mt: function_evaluations (3) failed")

      deallocate(minimizer_p)

   end subroutine minimizer_test

end module minimizer_tm