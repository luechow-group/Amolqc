! Copyright (C) 2020 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_ws_newton_module

   use kinds_m, only: r8
   use error_m, only: assert, asserts
   use fctn_module, only: Function_t
   use singularityCorrection_m, only: singularity_correction
   use singularityParticles_m, only: singularity_particles
   use minimizer_w_sing_module, only: minimizer_w_sing
   implicit none

   private
   public :: minimizer_ws_newton


   type, extends(minimizer_w_sing) :: minimizer_ws_newton
      real(r8) :: dt = 1._r8
      real(r8) :: delta_max = -1._r8
   contains
      procedure :: minimize => minimizer_ws_newton_minimize
      procedure :: write_params_minimizer_ws_newton
      procedure :: write_params => write_params_minimizer_ws_newton
   end type minimizer_ws_newton

   interface minimizer_ws_newton
      module procedure constructor1
   end interface minimizer_ws_newton

contains

   function constructor1(sc, dt, delta_max) result(minimizer)
      type(singularity_correction), intent(in)          :: sc
      real(r8), optional, intent(in) :: dt
      real(r8), optional, intent(in) :: delta_max
      type(minimizer_ws_newton), pointer :: minimizer
      allocate(minimizer)
      minimizer%sc_ = sc

      if (PRESENT(dt)) minimizer%dt = dt
      if (PRESENT(delta_max)) minimizer%delta_max = delta_max
   end function constructor1

   subroutine minimizer_ws_newton_minimize(this, fn, x)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_normal
      class(minimizer_ws_newton), intent(inout) :: this
      class(Function_t), intent(in)                  :: fn
      real(r8), intent(inout)                         :: x(:)
      real(r8)                   :: f, f_old
      real(r8)                   :: g(SIZE(x)), g_new(SIZE(x)), x_new(SIZE(x))
      real(r8)                   :: delta_x(SIZE(x)), H(SIZE(x),SIZE(x)), work(SIZE(x)**2)
      real(r8)                   :: ipiv(SIZE(x)), lambda(SIZE(X)), work2(3*SIZE(x)-1)
      real(r8)                   :: gmax, delta_x_scale
      integer                    :: iter, verbose, iul, i, j, k
      integer                    :: n, info, lwork
      type(singularity_particles):: sp
      logical                    :: is_corrected

      n = SIZE(x)

      if (verbose > 0) write(iul,"(a)") " *** newton minimizer with singularity information ***"
      call this%reset()

      if (asserts) call assert(mod(SIZE(x),3) == 0, "minimizer_ws_steep_desc_minimize: illegal size")

      verbose = this%verbose()
      iul = this%verbose_unit()

      call sp%create(SIZE(x)/3, this%sc_%n_singularities())

      ! before the first step, check for singularities, but with zero step.
      delta_x = 0._r8
      call this%sc_%correct_for_singularities(x, delta_x, sp, is_corrected, correction_only=.false.)

      ! these calculations are just for printing the initial position
      call fn%eval_fgh(x, f, g, H)
      call sp%Fix_gradients(g, H)

      if (verbose > 0) then
         write(iul,"(a,g20.10)") " initial position with function value f=", f
         write(iul,"(a,g20.10)") "                           and exp(-f) =", EXP(-f)
         call this%write_params(iul)
         if (verbose > 1) then
            write(iul,"(9g13.5)") x
            write(iul,"(9g13.5)") g
         end if
      end if

      if (this%do_write_opt_path()) call this%write_opt_path_entry(1, x, f)

      do iter = 1, this%max_iterations()
         if (verbose > 3) write(iul,"(a,i6)") " iter=", iter
         f_old = f

         ! calculating step
         call fn%eval_fgh(x, f, g, H)
         call sp%Fix_gradients(g, H)

         ! calculate delta_x = -H^-1 * g by solving H*delta_x = -g
         delta_x = - this%dt*g
         ! ?SYSV info: https://software.intel.com/en-us/node/468912

         if (.not. (ALL(ieee_is_normal(delta_x)) .and. ALL(ieee_is_normal(H)))) then
            call this%set_converged(.false.)
            exit
         end if

         call DSYSV('U',n,1,H,n,ipiv,delta_x,n,work,n**2,info)
         !call assert(info==0, 'newton_minimize: inversion failed')

         if (info /= 0 .or. .not. all(abs(delta_x)<huge(1.d0))) then
            call this%set_converged(.false.)
            exit
         end if

         ! saving position and value before step
         if (this%do_write_opt_path()) call this%write_opt_path_entry(iter, x, f)

         ! restriction of step length
         if (this%delta_max > 0._r8) then
            delta_x_scale = 1._r8
            do i = 1, SIZE(delta_x), 3
               delta_x_scale = MIN(delta_x_scale, this%delta_max / NORM2(delta_x(i:i+2)))
            end do
            delta_x = delta_x * delta_x_scale
         end if

         ! doing step (with singularity correction)
         call this%sc_%correct_for_singularities(x, delta_x, sp, is_corrected, correction_only=.false.)

         if (verbose > 4) write(iul,"(2i6,3g18.10)") iter, f, f_old

         gmax = maxval(abs(g))
         if (verbose > 2) then
            write(iul,"(i6,2(a,g18.10),(a,g12.3))") iter, " f=", f, " gmax=", gmax, " dt=", this%dt
            if (verbose > 3) then
               write(iul,"(9g13.5)") x
               write(iul,"(9g13.5)") g
            end if
         end if

         if (this%is_gradient_converged(gmax) .or. (sp%n_sing() == size(x)/3)) then
            call this%set_converged(.true.)
            exit
         end if

      end do

      if (verbose > 0) then
         write(iul,"(a,g20.10)") " final position with function value f=", f
         write(iul,"(a,g20.10)") "                         and exp(-f) =", EXP(-f)
         write(iul,"(a,l5,2(a,i8))") " converged=",this%is_converged(), " with iter=", iter
         if (verbose > 1) then
            write(iul,*) 'position:'
            write(iul,"(9g13.5)") x
            call fn%eval_fgh(x, f, g, H)
            call sp%Fix_gradients(g, H)
            write(iul,*) 'gradient:'
            write(iul,"(9g13.5)") g
            ! Getting eigenvalues and -vectors of Hessian
            lwork = 3*SIZE(X)-1
            call DSYEV('V', 'U', n, H, n, lambda, work2, lwork, info)
            !call assert(info==0, 'newton_minimize: diagonalization failed')
            write(iul,*) 'hessian eigenvalues and -vectors:'
            do i = 1, n
               write(iul,*) i, lambda(i)
               do k = 1, n/3
                  write(iul,'(i5,3f15.6)') k, ( H( j, i ), j = 3*k-2, 3*k )
               end do
            end do
         end if
      end if

      call this%set_value(f)
      call this%set_gradient(g)
      call this%set_iterations(iter)
      !call this%set_function_evaluations(n_eval)

   end subroutine minimizer_ws_newton_minimize

   subroutine write_params_minimizer_ws_newton(this, iu)
      class(minimizer_ws_newton), intent(in) :: this
      integer, intent(in), optional :: iu
      integer iull

      iull = this%verbose_unit()
      if (present(iu)) iull = iu
      write(iull,*) " * * * Newton minimizer with singularity correction * * *"

      if (present(iu)) then
         call this%write_params_minimizer_w_sing(iu) ! call parent class
      else
         call this%write_params_minimizer_w_sing()
      end if
   end subroutine write_params_minimizer_ws_newton

end module minimizer_ws_newton_module
