! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module fctn_module

   use kinds_m, only: r8
   use error_m, only: asserts, assert
   implicit none

   private
   public :: Function_t

   type, abstract :: Function_t
      real(r8) :: h = 1.e-5_r8  ! denominator for numerical derivatives
   contains
      procedure(eval_interface), deferred    :: eval
      procedure                              :: eval_fg => eval_fg_num
      procedure                              :: eval_fgh => eval_fgh_num
      procedure, non_overridable             :: eval_fg_num
      procedure, non_overridable             :: eval_fgh_num
      procedure, non_overridable             :: get_numerical_denominator
      procedure, non_overridable             :: set_numerical_denominator
   end type Function_t

   abstract interface
      function eval_interface(this, x) result(f)
         import :: r8
         import Function_t
         class(Function_t), intent(in)  :: this
         real(r8), intent(in) :: x(:)       ! coordinate
         real(r8)             :: f       ! function value
      end function eval_interface
   end interface


contains
   subroutine set_numerical_denominator(this, h)
      class(Function_t), intent(inout) :: this
      real(r8), intent(in)            :: h
      this%h = h
   end subroutine set_numerical_denominator

   function get_numerical_denominator(this) result(h)
      class(Function_t), intent(inout) :: this
      real(r8)                        :: h
      h = this%h
   end function get_numerical_denominator

   subroutine eval_fg_num(this, x, f, g)
      ! calculates numerical gradients using this%eval(x)
      ! currently simple symmetric 2-point formula with fixed h
      class(Function_t), intent(in)   :: this
      real(r8), intent(in)    :: x(:)    ! coordinate
      real(r8), intent(inout) :: f       ! function value
      real(r8), intent(inout) :: g(:)    ! gradient
      real(r8) :: xx(SIZE(x))
      real(r8) f1, f2
      integer i

      if (asserts) call assert(size(x) == size(g), "eval_fg_num: inconsistent sizes")

      xx = x
      do i = 1, size(x)
         xx(i) = x(i) + this%h
         f1 = this%eval(xx)
         xx(i) = x(i) - this%h
         f2 = this%eval(xx)
         g(i) = (f1 - f2) / (2.d0 * this%h)
         xx(i) = x(i)
      end do

      f = this%eval(x)
   end subroutine eval_fg_num

   subroutine eval_fgh_num(this, x, f, g, H)
      ! calculates numerical hessian using this%eval_fg(x)
      ! currently simple symmetric 2-point formula with fixed h
      class(Function_t), intent(in)   :: this
      real(r8), intent(in)    :: x(:)    ! coordinate
      real(r8), intent(inout) :: f       ! function value
      real(r8), intent(inout) :: g(:)    ! gradient
      real(r8), intent(inout) :: H(:,:)  ! hessian
      real(r8) :: xx(SIZE(x))
      real(r8), dimension(SIZE(g)) :: g1, g2
      real(r8) :: f_temp
      integer :: i

      if (asserts) then
         call assert(SIZE(x) == SIZE(g), "eval_fgh_num: inconsistent sizes")
         call assert(SIZE(x) == SIZE(H,1), "eval_fgh_num: inconsistent sizes")
         call assert(SIZE(x) == SIZE(H,2), "eval_fgh_num: inconsistent sizes")
      end if

      xx = x
      do i = 1, SIZE(x)
         xx(i) = x(i) + this%h
         call this%eval_fg(xx,f_temp,g1)
         xx(i) = x(i) - this%h
         call this%eval_fg(xx,f_temp,g2)
         H(i,:) = (g1 - g2) / (2.d0 * this%h)
         xx(i) = x(i)
      end do

      ! enforcing symmetry
      H = (H + TRANSPOSE(H)) / 2

      call this%eval_fg(x,f,g)
   end subroutine eval_fgh_num

end module fctn_module