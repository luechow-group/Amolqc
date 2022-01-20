! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_ws_factory_module

   use kinds_m, only: r8
   use error_m, only: error
   use parsing_m, only: getdbla, getstra, getinta, finda, getintarra
   use singularityCorrection_m, only: singularity_correction, NONE, CUTSTEP, UMRIGAR
   use line_search_ws_simple_module, only: line_search_ws_simple
   use line_search_weak_wolfe_module, only: line_search_weak_wolfe
   use minimizer_w_sing_module, only: minimizer_w_sing
   use minimizer_ws_steep_desc_module, only: minimizer_ws_steep_desc
   use minimizer_ws_fire_module, only: minimizer_ws_fire, fire_parameters, BACKTRACK
   use minimizer_ws_bfgs_module, only: minimizer_ws_bfgs
   use minimizer_ws_bfgst_module, only: minimizer_ws_bfgst
   use minimizer_ws_newton_module, only: minimizer_ws_newton
   use minimizer_ws_none_module, only: minimizer_ws_none
   implicit none

   private
   public :: create_ws_minimizer

contains 

   function create_ws_minimizer(lines) result(minimizer_p)
      character(len=*) :: lines(:)
      class(minimizer_w_sing), pointer :: minimizer_p
      real(r8) step_size, delta_max, step_max, gradient, sing_thresh, corr_thresh, alpha, cc, rho, c1, c2
      integer iflag, iflag1, iflag2, nlines, max_iter, mode, latency, switch_step
      logical yn, scaling
      character(len=20) value, string, str
      integer, allocatable :: not_to_minimize(:)
      type(fire_parameters) :: params
      type(line_search_ws_simple)     :: lss
      type(line_search_weak_wolfe)    :: lsww
      type(singularity_correction), pointer    :: sc_p

      minimizer_p => null()
      nlines = size(lines)

      call getdbla(lines, nlines, "singularity_threshold=", sing_thresh, iflag)
      if (iflag /= 0) sing_thresh = 0.005d0

      call getdbla(lines, nlines, "correction_threshold=", corr_thresh, iflag)

      if (iflag /= 0) corr_thresh = 0.1d0
      call getstra(lines, nlines, "correction_mode=", string, iflag)
      mode = CUTSTEP + 10
      if (iflag == 0) then
         if (string(1:3) == "cut") then
            mode = CUTSTEP
         else if (string(1:3) == "umr") then
            mode = UMRIGAR
         else if (string == "none") then
            mode = NONE
         else
            call error("illegal correction_mode")
         end if

         if (mode == CUTSTEP .or. mode == UMRIGAR) then
            if (string(4:7) == "_one") then
               continue
            else if (string(4:7) /= "   ") then
               call error("illegal correction_mode")
            else  ! cut the whole 3ND vector
               mode = mode + 10
            end if
         end if
      end if

      scaling = .not. finda(lines, nlines, "no_scaling")
      sc_p => singularity_correction(sing_thresh, corr_thresh, mode, scaling=scaling)

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

      call getdbla(lines, nlines, "max_distance=", delta_max, iflag1)
      call getdbla(lines, nlines, "max_distance_one=", step_max, iflag2)

      if (iflag1 == 0 .and. iflag2 == 0) call &
              error(" minimizer input: both max_distance and max_distance_one not allowed")

      value = "bfgs"
      call getstra(lines, nlines, "method=", value, iflag)

      if (value == "steepest_descent") then

         call getdbla(lines, nlines, "step_size=", step_size, iflag)
         if (iflag /= 0) step_size = 0.1d0
         if (iflag1 == 0) then
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho, delta_max_=delta_max)     ! note: alpha = 1.d0 refers to (- step_size * gradient)
         else if (iflag2 == 0) then
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho, step_max_=step_max)       ! note: alpha = 1.d0 refers to (- step_size * gradient)
         else
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho)                             ! note: alpha = 1.d0 refers to (- step_size * gradient)
         end if
         minimizer_p => minimizer_ws_steep_desc(lss, sc_p, step_size)
      else if (value == "newton") then

         call getdbla(lines, nlines, "step_size=", step_size, iflag)
         if (iflag /= 0) step_size = .02_r8

         minimizer_p => minimizer_ws_newton(sc_p, step_size, delta_max)
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
         call getstra(lines, nlines, "overshoot", string, iflag)
         if (string == "backtrack") params%overshoot=BACKTRACK

         if (iflag1 == 0) then
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho, delta_max_=delta_max)     ! note: alpha = 1.d0 refers to (- step_size * gradient)
         else if (iflag2 == 0) then
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho, step_max_=step_max)       ! note: alpha = 1.d0 refers to (- step_size * gradient)
         else
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho)                             ! note: alpha = 1.d0 refers to (- step_size * gradient)
         end if
         minimizer_p => minimizer_ws_fire(lss, sc_p, params)

      else if (value == "bfgst") then

         call getstra(lines, nlines, "scale_initial_H=", str, iflag)
         if (iflag == 0 .and. str == "y") then
            yn = .true.
         else
            yn = .false.
         endif

         if (iflag2 == 0) call error(" minimizer input: max_distance_one not implemented for bfgst")
         if (iflag1 == 0) then
            lsww = line_search_weak_wolfe(c1_=c1, c2_=c2, delta_max_=delta_max)
         else
            lsww = line_search_weak_wolfe(c1_=c1, c2_=c2)
         end if
         minimizer_p => minimizer_ws_bfgst(lsww, sc_p, yn)

      else if (value == "bfgs") then

         call getstra(lines, nlines, "scale_initial_H=", str, iflag)
         if (iflag == 0 .and. str == "y") then
            yn = .true.
         else
            yn = .false.
         endif

         if (iflag1 == 0) then
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho, delta_max_=delta_max)     ! note: alpha = 1.d0 refers to (- step_size * gradient)
         else if (iflag2 == 0) then
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho, step_max_=step_max)       ! note: alpha = 1.d0 refers to (- step_size * gradient)
         else
            lss = line_search_ws_simple(alpha_=alpha, c_=cc, rho_=rho)                             ! note: alpha = 1.d0 refers to (- step_size * gradient)
         end if
         step_size = 0.1d0
         call getdbla(lines, nlines, "step_size=", step_size, iflag)
         call getinta(lines, nlines, "latency=", latency, iflag2)
         switch_step = 50
         call getinta(lines, nlines, "switch_step=", switch_step, iflag)
         if (iflag2 == 0) then
            minimizer_p => minimizer_ws_bfgs(lss, sc_p, yn, step_size, latency=latency)
         else if (switch_step == 0) then
            minimizer_p => minimizer_ws_bfgs(lss, sc_p, yn)
         else
            minimizer_p => minimizer_ws_bfgs(lss, sc_p, yn, step_size, switch_step=switch_step)
         end if
      else if (value == "none") then
         minimizer_p => minimizer_ws_none()
      else
         call error("minimizer: method not implemented")
         minimizer_p => null()
      end if


      if (associated(minimizer_p)) then
         max_iter = 1000
         call getinta(lines, nlines, "max_iter=", max_iter, iflag)
         call minimizer_p%set_max_iterations(max_iter)
         gradient = 1.d-4
         call getdbla(lines, nlines, "convergence_gradient=", gradient, iflag)
         call minimizer_p%set_convergence_gradient(gradient)
         call getintarra(lines, nlines, "not_to_minimize=", not_to_minimize, iflag)
         if (iflag == 0) call minimizer_p%set_not_to_minimize(not_to_minimize)
      end if

   end function create_ws_minimizer

end module minimizer_ws_factory_module
