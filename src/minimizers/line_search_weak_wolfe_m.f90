! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module line_search_weak_wolfe_module

      ! fortran version of original matlab linesch_ww from hanso package
      ! implementing inexact find_step until weak Wolfe conditions are satisfied
      ! HANSO 2.0: GPL Michael Overton
   use kinds_m, only: r8
   use fctn_module, only: Function_t
   use line_search_m, only: line_search
   implicit none

   type, extends(line_search) :: line_search_weak_wolfe
      real(r8)    :: c1_ = 0.35d0   ! c1 and c2 values from weak Wolfe line search
      real(r8)    :: c2_ = 0.55d0   ! original values: c1 = 1.d-4, c2 = 0.9d0
   contains
      procedure :: find_step => find_step_hanso_ww
      procedure :: write_params => line_search_weak_wolfe_write_params
   end type line_search_weak_wolfe

   !   default constructor is sufficient
contains

   subroutine find_step_hanso_ww(this, ff, x0, falpha, grad0, d, xalpha, gradalpha, nfeval, gmask)
      ! fortran version of original matlab linesch_ww from hanso package
      ! implementing inexact find_step until weak Wolfe conditions are satisfied
      ! HANSO 2.0: GPL Michael Overton
      class(line_search_weak_wolfe), intent(in) :: this
      class(Function_t), intent(in)            :: ff  ! function object
      real(r8), intent(in)    :: x0(:)   ! position
      real(r8), intent(inout) :: falpha  ! function value (in: at x, out: at x_new)
      real(r8), intent(in)    :: grad0(:)! gradient
      real(r8), intent(in)    :: d(:)    ! search direction
      real(r8), intent(inout) :: xalpha(:)    ! new position
      real(r8), intent(inout) :: gradalpha(:) ! gradient at new position
      integer, intent(inout):: nfeval  ! out: required function evaluations, or error code (<0)
      logical, intent(in), optional :: gmask(:)  ! mask to force gradient component to zero (at singularities)
      real(r8)  :: alpha, beta   ! interval bracketing step length
      real(r8)  :: t             ! current steplength in x0 + t*d
      real(r8)  :: f0            ! old function value (f(x0))
      real(r8)  :: f             ! function value at x0 + t*d
      real(r8)  g0, dnorm, gtd
      real(r8), allocatable  :: x(:), grad(:)     ! trial position and gradient
      real(r8), allocatable  :: gradbeta(:)    ! gradient right bracket: [alpha, beta]
      integer i, nbisect, nexpand, nbisectmax, nexpandmax
      logical found

      allocate(x(size(x0)), grad(size(x0)))
      allocate(gradbeta(size(x0)))

      alpha = 0.d0    ! lower bound on steplength conditions
      xalpha = x0   
      f0 = falpha
      gradalpha = grad0
      beta = huge(f0)   ! upper bound on steplength satisfying weak Wolfe conditions
      gradbeta = huge(f0) 
      g0 = dot_product(grad0, d) 
      if (g0 >= 0) then
         nfeval = -1
         return
      end if

      dnorm = NORM2(d)
      if (dnorm == 0) then
         nfeval = -2 ! find_step_hanso_ww: direction vector is zero
         return
      end if

      t = 1.d0   ! important to try steplength one first
      if (this%delta_max_ > 0._r8) then
         do i = 1, SIZE(d), 3
            t = MIN(t, this%delta_max_ / NORM2(d(i:i+2)) )
         end do
      end if

      nfeval = 0
      nbisect = 0
      nexpand = 0
      ! the following limits are rather arbitrary
      ! nbisectmax = 30; % 50 is TOO BIG, because of rounding errors
      nbisectmax = max(30, nint(log(1e3*dnorm)/log(2.d0))) ! log_2! allows more if ||d|| big
      nexpandmax = max(10, nint(log(1e3/dnorm)/log(2.d0))) ! log_2! allows more if ||d|| small

      found = .false.
      do
         x = x0 + t*d
         nfeval = nfeval + 1
         call ff%eval_fg(x, f, grad)
         if (present(gmask)) then
            where (gmask) grad = 0.d0
         end if
         gtd = dot_product(grad, d)

         ! the first condition must be checked first. NOTE THE >=.
         if (f >= f0 + this%c1_ * t * g0) then   ! .or. IEEE_IS_NAN(f) !! first condition violated, gone too far
            beta = t
            gradbeta = grad            ! discard f
         ! now the second condition.  NOTE THE <=
         else if (gtd <= this%c2_ * g0) then   ! IEEE_IS_NAN(gtd) !! second condition violated, not gone far enough
            alpha = t
            xalpha = x
            falpha = f
            gradalpha = grad
         else                          ! quit, both conditions are satisfied
            alpha = t
            xalpha = x
            falpha = f
            gradalpha = grad
            beta = t
            gradbeta = grad
            found = .true.
            exit
         end if
         ! setup next function evaluation
         if (beta < huge(f0)) then
            if (nbisect < nbisectmax) then
               nbisect = nbisect + 1
               t = (alpha + beta) / 2.d0     ! bisection
            else
               exit
            end if
         else
            if (nexpand < nexpandmax) then
               nexpand = nexpand + 1
               t = 2.d0 * alpha     ! still in expansion mode
            else
               exit
            end if
         end if
      end do

      if (.not.found) then
         ! Wolfe conditions not satisfied: there are two cases
         if (beta == huge(f0)) then
            nfeval = -3 ! Line search failed to bracket point satisfying weak Wolfe conditions, function may be unbounded below
         else
            nfeval = -4 ! Line search failed to satisfy weak Wolfe conditions although point satisfying conditions was bracketed
         end if
      end if
   end subroutine find_step_hanso_ww

   subroutine line_search_weak_wolfe_write_params(this, iu)
      class(line_search_weak_wolfe), intent(in)    :: this
      integer, intent(in) , optional  :: iu
      integer iull

      if (present(iu)) then
         iull = iu
      else
         iull = this%verbose_unit()
      end if

      call this%line_search_write_params(iull)
      write(iull,"(3(a,g12.3))") " c1=", this%c1_, " c2=", this%c2_

   end subroutine line_search_weak_wolfe_write_params

end module line_search_weak_wolfe_module
