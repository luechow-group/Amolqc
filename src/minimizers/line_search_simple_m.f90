! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module line_search_simple_module
   use kinds_m, only: r8
   use fctn_module, only: Function_t
   use line_search_m, only: line_search
   use error_m, only: assert
   implicit none

   type, extends(line_search) :: line_search_simple
      real(r8)      :: alpha_ = 1.d0   ! line search parameter
      real(r8)      :: c_ = 0.1d0
      real(r8)      :: rho_ =0.33d0
   contains
      procedure :: find_step => find_step_simple
   end type line_search_simple

contains

   subroutine find_step_simple(this, ff, x0, falpha, grad0, d, xalpha, gradalpha, nfeval, gmask)
      class(line_search_simple), intent(in) :: this
      class(Function_t), intent(in)          :: ff  ! function object
      real(r8), intent(in)    :: x0(:)   ! position
      real(r8), intent(inout) :: falpha    ! function value (in: at x, out: at xalpha)
      real(r8), intent(in)    :: grad0(:)   ! gradient
      real(r8), intent(in)    :: d(:)   ! search direction
      real(r8), intent(inout) :: xalpha(:)   ! new position
      real(r8), intent(inout) :: gradalpha(:)   ! gradient at new position
      integer, intent(inout):: nfeval    ! in: maxiter, out: function evaluations or error code
      logical, intent(in), optional :: gmask(:)  ! should not be present here!

      real(r8) alpha, f, dir_deriv
      integer i, maxiter



      call assert(.not. PRESENT(gmask), 'find_step_simple: gmask should not be present here.')

      f = falpha     ! old function value
      alpha = this%alpha_

      ! restriction of step length
      if (this%delta_max_ > 0._r8) then
         do i = 1, SIZE(d), 3
            alpha = MIN(alpha, this%delta_max_ / NORM2(d(i:i+2)) )
         end do
      end if

      maxiter = nfeval

      dir_deriv = dot_product(grad0, d)    ! directional derivative

      if (dir_deriv >= 0.d0) then
         nfeval = -1
         return
      end if

      do nfeval = 1, maxiter

         xalpha = x0 + alpha * d

         call ff%eval_fg(xalpha, falpha, gradalpha)

      if (falpha < f + this%c_ * alpha * dir_deriv) exit

         alpha = this%rho_ * alpha

      end do

      !!!print*, "DBG:LSS2: ", nfeval, alpha

   end subroutine find_step_simple

end module line_search_simple_module

