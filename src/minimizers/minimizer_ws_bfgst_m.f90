! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_ws_bfgst_module

   use kinds_m, only: r8
   use fctn_module, only: Function_t
   use singularityCorrection_m, only: singularity_correction
   use singularityParticles_m, only: singularity_particles
   use line_search_weak_wolfe_module, only: line_search_weak_wolfe
   use minimizer_w_sing_module, only: minimizer_w_sing
   use error_m, only: abortp
   implicit none

   private
   public :: minimizer_ws_bfgst

   type, extends(minimizer_w_sing) :: minimizer_ws_bfgst
      private
      type(line_search_weak_wolfe)    :: lsww_
      logical                         :: scale_initial_H_ = .false.  ! .t. recommended by Nocedal
   contains
      procedure          :: minimize => minimizer_ws_bfgst_minimize
      procedure          :: write_params_minimizer_ws_bfgst
      procedure          :: write_params => write_params_minimizer_ws_bfgst
   end type minimizer_ws_bfgst

   interface minimizer_ws_bfgst
      module procedure constructor1
   end interface minimizer_ws_bfgst

contains

   function constructor1(lsww, sc, yn)
      class(line_search_weak_wolfe), intent(in)    :: lsww
      type(singularity_correction), intent(in)        :: sc
      logical, intent(in)                             :: yn
      type(minimizer_ws_bfgst), pointer :: constructor1
      allocate(constructor1)
      constructor1%lsww_ = lsww
      constructor1%sc_ = sc
      constructor1%scale_initial_H_ = yn
   end function constructor1

   subroutine minimizer_ws_bfgst_minimize(this, fn, x)
      class(minimizer_ws_bfgst), intent(inout) :: this
      class(Function_t), intent(in)         :: fn
      real(r8), intent(inout)                :: x(:)
      real(r8)                     :: f
      real(r8)                     :: g(size(x)), g_new(size(x))
      real(r8)                     :: gmax
      real(r8)                     :: H(size(x),size(x))    ! inverse Hessian approximation
      real(r8)                     :: p(size(x))    ! search direction
      real(r8)                     :: y(size(x))    ! notation see Nocedal/Wright book
      real(r8)                     :: x_new(size(x))
      real(r8)                     :: s(size(x))
      real(r8)                     :: H1(size(x),size(x))        ! aux matrix
      real(r8)                     :: v(size(x)), v1(size(x))    ! aux vectors
      real(r8)                     :: curv, rho
      real(r8)                     :: dx(3), ds(3)
      real(r8), allocatable        :: sings(:,:)
      integer                    :: i, iter, verbose, iul, n, ls_iter, info, nr_eval, n_sing_old
      type(singularity_particles):: sp
      logical                    :: gmask(size(x))                  ! .true. for components of particles at a singularity
      logical                    :: is_corrected, new_sing

      if (ALLOCATED(this%not_to_minimize_) .or. ALLOCATED(this%to_minimize_)) then
         call abortp('not_to_minimize and minimize_this are not implemented for bfgst')
      end if
      call sp%create(size(x)/3, this%sc_%n_singularities())

      verbose = this%verbose()
      iul = this%verbose_unit()
      if (verbose > 0) write(iul,"(a)") " *** BFGS minimizer ***"
      call this%reset()

      n = size(x)

      call fn%eval_fg(x, f, g)
      nr_eval = 1

      if (verbose > 0) then 
         write(iul,"(a,g20.10)") " initial position with function value f=", f
         if (verbose > 1) then
            call this%write_params(iul)
            call this%lsww_%write_params(iul)
            write(iul,"(9g13.5)") x
         end if
      end if

      ! initial inverse Hessian is unit matrix
      H = 0.d0
      do i = 1, n
         H(i, i) = 1.d0
      end do

      gmask = .false.
      new_sing = .false.

      do iter = 1, this%max_iterations()

         ! pseudo Newton direction
         p = - matmul(H, g)

         if (dot_product(p, g) > 0 .or. new_sing) then ! reset inverse Hessian and search direction
            if (verbose > 3) write(iul,*) "resetting inverse Hessian: p*g =", dot_product(p, g), new_sing
            H = 0.d0
            do i = 1, n
               H(i, i) = 1.d0
            end do
            p = - g
         end if

         ls_iter = 20

         !! line search with distance restriction
         call this%lsww_%find_step(fn, x, f, g, p, x_new, g_new, ls_iter, gmask)

         if (ls_iter < 0) then
            info = - ls_iter
            exit
         else
            nr_eval = nr_eval + ls_iter
         end if

         if (verbose > 3) then
            write(iul,"(2(a,i6),a,g18.10)") " iter=", iter, " ls_iter=", ls_iter, " f=", f
            if (verbose > 4) then
               write(iul,*) " x_new, g_new:"
               write(iul,'(9g16.8)') x_new
               write(iul,'(9g16.8)') g_new
            end if
         end if


         ! correction step with same length (per particle) as Newton step, but along new gradient direction
         p = - g_new
         s = x_new - x   ! BFGS (notation as in Nocedal/Wright): ignore singularity correction step
         y = g_new - g   ! BFGS: ignore singularity correction step

         do i = 1, size(x), 3
            dx = p(i:i+2)
            ds = s(i:i+2)
            p(i:i+2) = dx * sqrt( dot_product(ds, ds) / dot_product(dx, dx) )
         end do

         n_sing_old = sp%n_sing()
         call this%sc_%correct_for_singularities(x_new, p, sp, is_corrected, correction_only=.true.)
         new_sing = (sp%n_sing() /= n_sing_old)

         if (is_corrected) then
            call fn%eval_fg(x_new, f, g_new)
            nr_eval = nr_eval + 1
            gmask = sp%At_singularity()
            where (gmask) g_new = 0.d0
         end if

         if (verbose > 3) then
            write(iul,"(a)") "after correction step:"
            write(iul,"(a,g18.10,a,i4)") " with f=", f, " n_sing=", sp%n_sing()
            if (verbose > 4) then
               write(iul,'(9g16.8)') x_new
               write(iul,'(9g16.8)') g_new
            end if
         end if

         gmax = maxval(abs(g_new))

         if (verbose > 2) then
            write(iul,"(i6,a,g20.11,a,g15.5,2(a,i3),a,l6)") iter, " f=", f, " gmax=", gmax, " ls_iter=", ls_iter,  &
                " n_sing=", sp%n_sing(), " is_corrected=", is_corrected
         end if

         if (this%is_gradient_converged(gmax) .or. this%is_value_converged(f)) then
            x = x_new
            g = g_new
            call this%set_converged(.true.)
            exit
         end if


         ! BFGS update of inverse Hessian (see Nocedal/Wright book)

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
   end subroutine minimizer_ws_bfgst_minimize

   subroutine write_params_minimizer_ws_bfgst(this, iu)
      class(minimizer_ws_bfgst), intent(in) :: this
      integer, intent(in), optional :: iu
      integer iull

      iull = this%verbose_unit()
      if (present(iu)) iull = iu
      write(iull,*) " * * * BFGS minimizer with singularity correction (Dahl version) * * *"

      if (present(iu)) then
         call this%write_params_minimizer_w_sing(iu) ! call parent class
         call this%lsww_%write_params(iu)
      else
         call this%write_params_minimizer_w_sing()
         call this%lsww_%write_params()
      end if
      write(iull, "(a,l5)")  " scale_initial_H     =", this%scale_initial_H_
   end subroutine write_params_minimizer_ws_bfgst


end module minimizer_ws_bfgst_module
