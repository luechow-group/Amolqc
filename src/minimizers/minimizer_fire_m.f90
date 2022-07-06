! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_fire_module

   use kinds_m, only: r8
   use fctn_module, only: Function_t
   use minimizer_module, only: minimizer
   implicit none

   private
   public :: minimizer_fire, fire_parameters

   integer, parameter, public  :: PASS=1, BACKTRACK=2

   type fire_parameters
      real(r8) :: tau_init  = 0.d0    ! FIRE parameters to be set
      real(r8) :: tau_max   = 0.d0
      real(r8) :: alpha0    = 0.1d0   ! default FIRE parameters
      real(r8) :: f_alpha   = 0.99d0
      real(r8) :: f_inc     = 1.1d0
      real(r8) :: f_dec     = 0.5d0
      integer :: n_min    = 5
      integer :: overshoot = PASS
   end type fire_parameters

   type, extends(minimizer) :: minimizer_fire
      type(fire_parameters) :: params_
   contains
      procedure      :: minimize => minimizer_fire_minimize
      procedure      :: write_params_minimizer_fire
      procedure      :: write_params => write_params_minimizer_fire
   end type minimizer_fire

   interface minimizer_fire
      module procedure constructor1
   end interface minimizer_fire

contains

   function constructor1(params)
      type(fire_parameters), intent(in) :: params
      type(minimizer_fire), pointer :: constructor1
      allocate(constructor1)
      constructor1%params_ = params
   end function constructor1

   subroutine minimizer_fire_minimize(this, fn, x)
      class(minimizer_fire), intent(inout) :: this
      class(Function_t), intent(in)         :: fn
      real(r8), intent(inout)                :: x(:)
      real(r8)                     :: f, f_old
      real(r8)                     :: delta_x(size(x))
      real(r8)                     :: g(size(x)), velocity(size(x)), force(size(x))
      real(r8)                     :: x_old(size(x)), force_old(size(x)), velocity_old(size(x))
      real(r8)                     :: force0(size(x))
      real(r8)                     :: alpha, tau, mass, gmax, vel, proj
      integer                    :: iter, verbose, iul, latency

      verbose = this%verbose()
      iul = this%verbose_unit()
      if (verbose > 0) write(iul,"(a)") " *** FIRE minimizer ***"
      call this%reset()

      call fn%eval_fg(x, f, g)

      force = -g
      velocity = 0
      alpha = this%params_%alpha0
      tau = this%params_%tau_init
      latency = 0
      mass = 1.d0    ! fictitious mass

      if (verbose > 0) then 
         write(iul,"(a,g20.10)") " initial position with function value f=", f
         call this%write_params(iul)
         if (verbose > 1) then
            write(iul,"(9g13.5)") x
         end if
      end if

      do iter = 1, this%max_iterations()

         ! save for backtracking
         x_old = x; f_old = f
         force_old = force
         velocity_old = velocity

         ! velocity Verlet step, a = F/m
         delta_x = velocity * tau + 0.5d0 * force/mass * tau**2

         x = x + delta_x

         call fn%eval_fg(x, f, g)

         gmax = maxval(abs(g))

         if (verbose > 2) then
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

         force = -g
         velocity = velocity + 0.5d0 * (force + force_old) / mass * tau

         ! FIRE F1
         proj = dot_product(force, velocity)

         ! FIRE F2
         vel = sqrt( dot_product(velocity, velocity) )
         force0 = force / sqrt( dot_product(force, force) )
         velocity = (1.d0 - alpha) * velocity + alpha * vel * force0

         if (verbose > 3) then
            write(iul,'(a,2g18.8)') "P, |v| = ", proj, vel
         end if

         ! FIRE F3/4

         if (proj > 0.d0 .and. latency > this%params_%n_min) then
            tau = min( this%params_%f_inc * tau, this%params_%tau_max )
            alpha = this%params_%f_alpha * alpha
            if (verbose > 3) then
               write(iul,'(a,i5,2g18.8)') " incr tau; latency, tau, alpha = ", latency, tau, alpha
            end if
         end if

         !!!if (proj <= 0.d0 .or. isSet) then  ! uphill or new elec at nuc: freeze velocity and reset alpha
         if (proj <= 0.d0) then  ! uphill: freeze velocity and reset alpha
            tau = this%params_%f_dec * tau
            velocity = 0.d0
            alpha = this%params_%alpha0
            latency = 0
            if (verbose > 3) then
               !!!if (isSet) then
               !!!   write(iull,'(a,2g18.8)') " freeze at nuc; tau, alpha = ", tau, alpha
               !!!else
                  write(iul,'(a,2g18.8)') " freeze uphill; tau, alpha = ", tau, alpha
               !!!end if
            end if
            !!!if (.not. isSet .and. params%overshoot == BACKTRACK) then
            if (this%params_%overshoot == BACKTRACK) then
               x = x_old; force = force_old; velocity = velocity_old; f = f_old
            end if
         else
            latency = latency + 1
         end if

      end do

      call this%set_value(f)
      call this%set_gradient(g)
      call this%set_iterations(iter)
      call this%set_function_evaluations(iter + 1)
   end subroutine minimizer_fire_minimize

   subroutine write_params_minimizer_fire(this, iu)
      class(minimizer_fire), intent(in) :: this
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
      write(iull,"(a,g13.5,a,i3)")  " tau_init=", this%params_%tau_init, " overshoot=",this%params_%overshoot
      write(iull,"(a,g13.5,a,g13.5)")  " tau_max=", this%params_%tau_max, " alpha0=", this%params_%alpha0
      write(iull,"(a,g13.5,a,g13.5)")  " f_alpha=", this%params_%f_alpha, " f_inc=", this%params_%f_inc
      write(iull,"(a,g13.5,a,i3)")  " f_dec=", this%params_%f_dec, " latency=", this%params_%n_min
   end subroutine write_params_minimizer_fire


end module minimizer_fire_module
