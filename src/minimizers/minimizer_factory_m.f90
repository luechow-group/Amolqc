! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_factory_module

   use kinds_m, only: r8
   use error_m, only: error
   use parsing_m, only: getdbla, getstra, getinta
   use minimizer_module, only: minimizer
   use minimizer_steep_desc_module, only: minimizer_steep_desc
   use minimizer_fire_module, only: minimizer_fire, fire_parameters, BACKTRACK
   use minimizer_bfgs_module, only: minimizer_bfgs
   use line_search_simple_module, only: line_search_simple
   use line_search_weak_wolfe_module, only: line_search_weak_wolfe
   implicit none

   private
   public :: create_minimizer

contains 

   function create_minimizer(lines) result(minimizer_p)
      character(len=*) :: lines(:)
      class(minimizer), pointer :: minimizer_p
      real(r8) step_size, delta_max, gradient, alpha, cc, rho, c1, c2, val
      integer iflag, max_iter, v, nlines, iflag1
      logical yn
      character(len=20) value, str
      type(line_search_simple)     :: lss
      type(line_search_weak_wolfe) :: lsww
      type(fire_parameters) :: params

      nlines = size(lines)

      alpha = 1.d0
      call getdbla(lines, nlines, "alpha=", alpha, iflag)
      cc = 0.1d0
      call getdbla(lines, nlines, "cc=", cc, iflag)
      rho = 0.33d0
      call getdbla(lines, nlines, "rho=", rho, iflag)
      c1 = 0.35d0
      call getdbla(lines, nlines, "c1=", c1, iflag)
      c2 = 0.55d0
      call getdbla(lines, nlines, "c2=", c2, iflag)

      call getstra(lines, nlines, "method=", value, iflag)
      if (iflag /= 0) call error("minimizer: missing method key")

      if (value == "steepest_descent") then
         call getdbla(lines, nlines, "step_size=", step_size, iflag)
         if (iflag /= 0) step_size = 0.1d0
         call getdbla(lines, nlines, "max_distance=", delta_max, iflag)
         if (iflag == 0) then
            lss = line_search_simple(alpha_=alpha, c_=cc, rho_=rho, delta_max_=delta_max)     ! note: alpha = 1.d0 refers to (- step_size * gradient)
         else
            lss = line_search_simple(alpha_=alpha, c_=cc, rho_=rho)
         end if

         minimizer_p => minimizer_steep_desc(lss, step_size)

      else if (value == "fire") then
         call getdbla(lines, nlines, "tau_init=", params%tau_init, iflag)
         if (iflag /= 0) call error("method=fire: tau_init required")
         call getdbla(lines, nlines, "tau_max=", params%tau_max, iflag)
         if (iflag /= 0) params%tau_max = 10*params%tau_init
         call getdbla(lines, nlines, "alpha0=", params%alpha0, iflag)
         call getdbla(lines, nlines, "f_alpha=", params%f_alpha, iflag)
         call getdbla(lines, nlines, "f_inc=", params%f_inc, iflag)
         call getdbla(lines, nlines, "f_dec=", params%f_dec, iflag)
         call getinta(lines, nlines, "latency=", params%n_min, iflag)
         call getstra(lines, nlines, "overshoot", str, iflag)
         if (str == "backtrack") params%overshoot=BACKTRACK

         minimizer_p => minimizer_fire(params)

      else if (value == "bfgs") then
         call getstra(lines, nlines, "scale_initial_H=", str, iflag)
         if (iflag == 0 .and. str == "y") then
            yn = .true.
         else
            yn = .false.
         endif
         call getdbla(lines, nlines, "max_distance=", delta_max, iflag1)
         if (iflag1 == 0) then
            lsww = line_search_weak_wolfe(c1_=c1, c2_=c2, delta_max_=delta_max)
         else
            lsww = line_search_weak_wolfe(c1_=c1, c2_=c2)
         end if         

         minimizer_p => minimizer_bfgs(lsww, yn)
      else
         call error(" create_minimizer: illegal argument")
         minimizer_p => null()
      end if

      if (associated(minimizer_p)) then
         call getinta(lines, nlines, "max_iter=", max_iter, iflag)
         if (iflag == 0) call minimizer_p%set_max_iterations(max_iter)
         call getdbla(lines, nlines, "convergence_gradient=", gradient, iflag)
         if (iflag == 0) call minimizer_p%set_convergence_gradient(gradient)
         val = -HUGE(0._r8)
         call getdbla(lines, nlines, "convergence_value=", val, iflag)
         call minimizer_p%set_convergence_value(val)
         call getinta(lines, nlines, "verbose=", v, iflag)
         if (iflag == 0) call minimizer_p%set_verbose(v)
      end if
   end function create_minimizer

end module minimizer_factory_module
