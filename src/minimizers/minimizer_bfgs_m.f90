! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_bfgs_module

   use kinds_m, only: r8
   use fctn_module, only: Function_t
   use minimizer_module, only: minimizer
   use line_search_weak_wolfe_module, only: line_search_weak_wolfe
   implicit none

   private
   public :: minimizer_bfgs

   type, extends(minimizer) :: minimizer_bfgs
      private
      type(line_search_weak_wolfe) :: lnsrch_                   ! linesearch function object
      logical                    :: scale_initial_H_ = .false.  ! .t. recommended by Nocedal
   contains
      procedure                  :: minimize => minimizer_bfgs_minimize
      procedure                  :: write_params_minimizer_bfgs
      procedure                  :: write_params => write_params_minimizer_bfgs
   end type minimizer_bfgs

   interface minimizer_bfgs
      module procedure constructor
   end interface minimizer_bfgs

contains

   function constructor(ls, yn)
      type(line_search_weak_wolfe), intent(in)    :: ls
      logical, intent(in)          :: yn
      type(minimizer_bfgs), pointer :: constructor
      allocate(constructor)
      constructor%lnsrch_ = ls
      constructor%scale_initial_H_ = yn
   end function constructor

   subroutine minimizer_bfgs_minimize(this, fn, x)
      class(minimizer_bfgs), intent(inout) :: this
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
      real(r8) curv, rho
      integer i, iter, verbose, iul, n, ls_iter, info, nr_eval

      verbose = this%verbose()
      iul = this%verbose_unit()
      if (verbose > 0) write(iul,"(a)") " *** BFGS minimizer ***"
      call this%reset()

      n = size(x)

      call fn%eval_fg(x, f, g)
      nr_eval = 1

      if (verbose > 0) then 
         write(iul,"(a,g20.10)") " initial position with function value f=", f
         call this%write_params(iul)
         if (verbose > 1) then
            write(iul,"(9g13.5)") x
         end if
      end if

      ! initial inverse Hessian is unit matrix
      H = 0.d0
      do i = 1, n
         H(i, i) = 1.d0
      end do

      do iter = 1, this%max_iterations()

         ! pseudo Newton direction
         p = - matmul(H, g)

         ls_iter = 20

         ! line search with distance restriction
         call this%lnsrch_%find_step(fn, x, f, g, p, x_new, g_new, ls_iter)

         if (ls_iter < 0) then
            info = - ls_iter
            exit
         else
            nr_eval = nr_eval + ls_iter
         end if

         gmax = maxval(abs(g_new))

         if (verbose > 2) then
            write(iul,"(i6,a,i3,a,i6,2(a,g18.10))") iter, " line search iter=", ls_iter, " total fctn eval=", nr_eval, &
               " f=", f, " gmax=", gmax
            if (verbose > 3) then
               write(iul,'(9g13.5)') x_new
               write(iul,'(9g13.5)') g_new
            end if
         end if

         if (this%is_gradient_converged(gmax)) then
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
            if (verbose > 2) write(iul,*) "curvature condition failed"
            info = 2
            return
         end if

         ! as recommended by Nocedal (8.20) scale initial inverse Hessian
         if (iter == 1 .and. this%scale_initial_H_) then
            H = curv/dot_product(y, y) * H
         end if

         ! BFGS update (8.16 in Nocedal/Wright book)
         ! avoiding O(n^3) matmul (expanding 8.16)

         rho = 1.d0 / curv

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

      call this%set_value(f)
      call this%set_gradient(g)
      call this%set_iterations(iter)
      call this%set_function_evaluations(nr_eval)
   end subroutine minimizer_bfgs_minimize

   subroutine write_params_minimizer_bfgs(this, iu)
      class(minimizer_bfgs), intent(in) :: this
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
      write(iull,"(a,l5)") " scale_initial_H=", this%scale_initial_H_
   end subroutine write_params_minimizer_bfgs



end module minimizer_bfgs_module
