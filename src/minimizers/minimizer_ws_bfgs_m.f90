! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_ws_bfgs_module

   use kinds_m, only: r8
   use error_m, only: assert, asserts, debug, abortp
   use fctn_module, only: Function_t
   use singularityCorrection_m, only: singularity_correction
   use singularityParticles_m, only: singularity_particles
   use line_search_ws_m, only: line_search_ws
   use minimizer_w_sing_module, only: minimizer_w_sing
   implicit none

   private
   public :: minimizer_ws_bfgs

   type, extends(minimizer_w_sing) :: minimizer_ws_bfgs
      private
      class(line_search_ws), allocatable     :: lss_
      logical                         :: scale_initial_H_ = .false.  ! .t. recommended by Nocedal
      real(r8)                    :: step_size_ = 0.d0
      integer                         :: latency_ = 0
      integer                         :: switch_step_ = 0
   contains
      procedure          :: minimize => minimizer_ws_bfgs_minimize
      procedure          :: write_params_minimizer_ws_bfgs
      procedure          :: write_params => write_params_minimizer_ws_bfgs
   end type minimizer_ws_bfgs

   interface minimizer_ws_bfgs
      module procedure constructor
   end interface minimizer_ws_bfgs

contains

   function constructor(lss, sc, yn, step_size, latency, switch_step)
      class(line_search_ws), intent(in)        :: lss
      type(singularity_correction), intent(in)        :: sc
      logical, intent(in)                             :: yn
      real(r8), optional, intent(in)                  :: step_size
      integer, optional, intent(in)                   :: latency
      integer, optional, intent(in)                   :: switch_step
      type(minimizer_ws_bfgs), pointer :: constructor
      allocate(constructor)
      constructor%lss_ = lss
      constructor%sc_ = sc
      constructor%scale_initial_H_ = yn
      if (present(step_size)) constructor%step_size_ = step_size
      if (present(latency)) constructor%latency_ = latency
      if (present(switch_step)) constructor%switch_step_ = switch_step
   end function constructor

   subroutine minimizer_ws_bfgs_minimize(this, fn, x)
      class(minimizer_ws_bfgs), intent(inout) :: this
      class(Function_t), intent(in)         :: fn
      real(r8), intent(inout)                :: x(:)
      real(r8)                     :: f
      real(r8)                     :: g(size(x)), g_new(size(x))
      real(r8)                     :: d, gmax
      real(r8)                     :: H(size(x),size(x))    ! inverse Hessian approximation
      real(r8)                     :: p(size(x))    ! search direction
      real(r8)                     :: y(size(x))    ! notation see Nocedal/Wright book
      real(r8)                     :: x_new(size(x))
      real(r8)                     :: s(size(x))
      real(r8)                     :: H1(size(x),size(x))        ! aux matrix
      real(r8)                     :: v(size(x)), v1(size(x))    ! aux vectors
      real(r8)                     :: curv, rho
      real(r8), allocatable        :: sings(:,:)
      real(r8)                     :: xyz(SIZE(x)/3, 3)
      integer                    :: i, iter, verbose, iul, n, ls_iter, info, nr_eval, max_d_flag
      integer                    :: steps_since_correct
      type(singularity_particles):: sp
      logical                    :: gmask(size(x)), mask(size(x))   ! .true. for components of particles at a singularity
      logical is_corrected, new_sing

      if (asserts) then
         call assert(all(abs(x)<huge(1.d0)), "minimizer_ws_bfgs_minimize: illegal x coords")
         call assert(size(x) == 3*(size(x)/3), "minimizer_ws_bfgs_minimize: illegal x size")
      end if

      call sp%create(size(x)/3, this%sc_%n_singularities())

      verbose = this%verbose()
      iul = this%verbose_unit()
      call this%sc_%set_verbose(verbose)
      call this%sc_%set_verbose_unit(iul)

      if (verbose > 0) write(iul,"(a)") " *** BFGS minimizer ***"
      call this%reset()

      n = size(x)

      nr_eval = 1

      ! before the first step, check for singularities, but with zero step and set particles to singularities
      p = 0._r8
      call this%sc_%correct_for_singularities(x, p, sp, is_corrected, correction_only=.true.)
      mask = sp%At_singularity()
      call fn%eval_fg(x, f, g, mask)

      where ( sp%At_singularity() ) g = 0.d0

      call this%restrict_gradient(g)

      if (verbose > 0) then 
         write(iul,"(a,g20.10)") " initial position with function value f=", f
         if (verbose > 1) then
            write(iul,"(a,g13.6,a,i8)") " max_iterations = ", this%max_iterations()
            call this%write_params(iul)
            write(iul,"(9g13.5)") x
         end if
      end if

      if (this%do_write_opt_path()) call this%write_opt_path_entry(1, x, f)

      ! initial inverse Hessian is unit matrix
      H = 0.d0
      do i = 1, n
         H(i, i) = 1.d0
      end do

      gmask = .false.
      new_sing = .false.
      is_corrected = .false.
      steps_since_correct = 0

      do iter = 1, this%max_iterations()

         ! pseudo Newton direction
         p = - matmul(H, g)

         if (dot_product(p, g) > 0  &
            .or. steps_since_correct <= this%latency_  &
            .or. iter <= this%switch_step_) then 
            ! reset inverse Hessian and search direction
            if (verbose > 3) write(iul,*) "resetting inverse Hessian: p*g =", dot_product(p, g), new_sing
            H = 0.d0
            do i = 1, n
               H(i, i) = 1.d0
            end do
            p = - g
         end if

         ls_iter = 5


         where (gmask) p = 0.d0

         ! before 'latency' continue with steepest descent using step_size
         ! after 'latency' go smoothly from steepest descent to full BFGS
         ! alternatively use steepest descent up to 'switch_step' then switch smoothly to BFGS
         if (this%step_size_ > 0.d0) then
            if (this%switch_step_ > 0) then
               if (iter <= this%switch_step_) then
                  p = - this%step_size_ * g
               else if (iter <= this%switch_step_ + 4) then
                  d = iter - this%switch_step_
                  p = (1.d0 - d/4.d0) * (- this%step_size_ * g) + d/4.d0 * p
               end if
            else if (this%latency_ > 0) then
               if (steps_since_correct <= this%latency_) then
                  p = - this%step_size_ * g
               else if (steps_since_correct < this%latency_ + 4) then
                  d = steps_since_correct - this%latency_
                  p = (1.d0 - d/4.d0) * (- this%step_size_ * g) + d/4.d0 * p
               end if
            end if
         end if

         if (debug) then
            if (.not.all(abs(x)<huge(1.d0))) then
               call abortp("minimizer_ws_bfgs_minimize: illegal x coords before find_step")
            end if
            if (.not.all(abs(p)<huge(1.d0))) then
               call abortp("minimizer_ws_bfgs_minimize: illegal p coords before find_step")
            end if
         end if

         !! line search with distance restriction
         call this%restrict_gradient(g)

         call this%lss_%find_step(fn, this%sc_, x, f, g, p, x_new, g_new, ls_iter, is_corrected, sp)
         call this%restrict_gradient(g_new)


         if (debug) call assert(all(abs(x)<huge(1.d0)), "minimizer_ws_bfgs_minimize: illegal x coords after find_step")

         if (ls_iter < 0) then
            info = - ls_iter
            if (verbose > 1) write(iul,"(a)") " line search unsuccessful! Giving up"
            !!!!exit
         else
            nr_eval = nr_eval + ls_iter
         end if

         if (is_corrected) then
            steps_since_correct = 0
         else
            steps_since_correct = steps_since_correct + 1
         end if

         gmask = sp%At_singularity()
         where (gmask) g_new = 0.d0

         if (verbose > 3) then
            write(iul,"(2(a,i6),a,g18.10,a,i4)") " iter=", iter, " ls_iter=", ls_iter, " f=", f, " n_sing=", sp%n_sing()
            if (verbose > 4) then
               write(iul,*) " x_new, g_new:"
               write(iul,'(9g16.8)') x_new
               write(iul,'(9g16.8)') g_new
            end if
         end if

         if (this%do_write_opt_path()) call this%write_opt_path_entry(iter + 1, x_new, f)

         gmax = maxval(abs(g_new))

         if (verbose > 2) then
            write(iul,"(i6,a,g20.11,a,g15.5,2(a,i3),a,l6,a,i7)") iter, " f=", f, " gmax=", gmax, " ls_iter=", ls_iter,  &
                " n_sing=", sp%n_sing(), " is_corrected=", is_corrected, " steps_since_correct=", steps_since_correct
         end if


         if (this%is_gradient_converged(gmax)) then
            if (this%value_convergence_) then
               call this%set_converged(.false.)
               exit
            else
               x = x_new
               g = g_new
               call this%set_converged(.true.)
               exit
            end if
         end if

         if (this%is_value_converged(f) .and. this%value_convergence_) then
            x = x_new
            g = g_new
            call this%set_converged(.true.)
            exit
         end if


         ! BFGS update of inverse Hessian (see Nocedal/Wright book)

         s = x_new - x
         y = g_new - g

         curv = dot_product(y, s)

         ! check curvature condition
         if (curv <= 0.d0) then
            if (verbose > 2) write(iul,"(a,g20.5)") "curvature condition failed (ignored!):", curv
            info = 2
            !!!!cycle
         end if

         ! as recommended by Nocedal (8.20) scale initial inverse Hessian
         if (iter == 1 .and. this%scale_initial_H_) then
            H = curv/dot_product(y, y) * H
         end if

         ! BFGS update (8.16 in Nocedal/Wright book)
         ! avoiding O(n^3) matmul (expanding 8.16)

         rho = 1.d0 / curv
         if (verbose > 3) then
            write(iul,"(a,g18.10)") " rho=", rho
            if (verbose > 4) then
               write(iul,"(a)") "s,y after correction step:"
               write(iul,'(9g16.8)') s
               write(iul,'(9g16.8)') y
            end if
         end if

         v = matmul(H, y)

         ! outer product (H*y)*s^T
         do i = 1, n
            H1(:, i) = rho * v(:) * s(i)
         end do

         ! save y^T*H
         v = matmul(y, H)

         ! save y^T*H1 = y^T * (H*y)*s^T
         v1 = matmul(y, H1)

         H = H - H1

         ! outer product s*(y^T*H)
         do i = 1, n
            H1(:, i) = rho * s(:) * v(i)
         end do

         H = H - H1

         ! outer product s*(y^T*H1)
         do i = 1, n
            H1(:, i) = rho * s(:) * v1(i)
         end do

         H = H + H1

         ! outer product s*s^T
         do i = 1, n
            H1(:, i) = rho * s(:) * s(i)
         end do

         H = H + H1

         x = x_new
         g = g_new

         max_d_flag = 0

         if (MAXVAL(ABS(x_new)) > MINVAL(this%max_electron_distance())) then
            xyz = RESHAPE(x_new, [SIZE(x_new)/3,3])
            do i = 1, 3
               if (MAXVAL(ABS(xyz(:,i))) > MINVAL(this%max_electron_distance())) then
                  max_d_flag = 1
                  call this%set_converged(.false.)
                  exit
               end if
            end do
            if (max_d_flag == 1) exit
         end if

      end do

      if (verbose > 0) then
         write(iul,"(a,g22.14)") " final position with function value f=", f
         write(iul,"(a,l5,2(a,i8))") " converged=",this%is_converged(), " with iter=", iter, " fctn_eval=", nr_eval
         if (verbose > 1) then
            write(iul, *) " electron positions:"
            write(iul,"(3g22.14)") x
            write(iul, *) " nuclei positions:"
            sings = this%singularities()
            do i = 1, size(sings, dim=2)
               write(iul,"(3g22.14)") sings(1:3, i)
            end do
            write(iul, *) " gradient (0 at singularity):"
            write(iul,"(9g16.8)") g
         end if
      end if

      call this%set_value(f)
      call this%set_gradient(g)
      call this%set_iterations(iter)
      call this%set_function_evaluations(nr_eval)
   end subroutine minimizer_ws_bfgs_minimize

   subroutine write_params_minimizer_ws_bfgs(this, iu)
      class(minimizer_ws_bfgs), intent(in) :: this
      integer, intent(in), optional :: iu
      integer iull

      iull = this%verbose_unit()
      if (present(iu)) iull = iu
      write(iull,"(a)") " -- BFGS minimizer with singularity correction"

      if (present(iu)) then
         call this%write_params_minimizer_w_sing(iu) ! call parent class
         call this%lss_%write_params(iu)
      else
         call this%write_params_minimizer_w_sing()
         call this%lss_%write_params()
      end if
      write(iull, "(/a)") " BFGS parameters:"
      write(iull, "(a,l5)")  " scale_initial_H     =", this%scale_initial_H_
      if (this%step_size_ > 0.d0) write(iull, "(a,g12.5)")  &
                             " step                =", this%step_size_
      if (this%latency_ > 0) write(iull, "(a,i12)")   &
                             " latency             =", this%latency_
      if (this%switch_step_ > 0) write(iull, "(a,i8)")  &
                             " switch_step         =", this%switch_step_
   end subroutine write_params_minimizer_ws_bfgs


end module minimizer_ws_bfgs_module
