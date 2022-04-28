! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module myfctn_module

   use kinds_m, only: r8
   use error_m, only: assert, asserts
   use fctn_module, only: Function_t
   implicit none

   type, extends(Function_t) :: myFunction_t
      real(r8) :: a
   contains
      procedure :: eval => myfctn_eval
      procedure :: eval_fg => myFunction_t_eval
      procedure :: set => myFunction_t_set
      procedure :: get => myFunction_t_get
   end type myFunction_t

   type, extends(Function_t) :: mysingFunction_t
      real(r8) :: a, b
      real(r8) :: za = -1.d0
      real(r8) :: zb = +1.d0
   contains
      procedure :: eval => mysingfctn_eval
      procedure :: eval_fg => mysingFunction_t_eval
      procedure :: set => mysingFunction_t_set
      procedure :: get => mysingFunction_t_get
   end type mysingFunction_t

contains

   !!! f(x) = sum_{i=1}^n a * x_i^4

   subroutine myFunction_t_set(this, a)
      class(myFunction_t), intent(inout) :: this
      real(r8), intent(in) :: a
      this%a = a
      call assert(a > 0, "myfctn: a > 0 required")
   end subroutine myFunction_t_set

   subroutine myFunction_t_get(this, a)
      class(myFunction_t), intent(in) :: this
      real(r8), intent(inout) :: a
      a = this%a
   end subroutine myFunction_t_get

   subroutine myFunction_t_eval(this, x, f, g)
      class(myFunction_t), intent(in) :: this
      real(r8), intent(in) :: x(:)       ! coordinate
      real(r8), intent(out) :: f       ! function value
      real(r8), intent(out) :: g(:)    ! gradient
      integer i
      f = 0.d0
      do i = 1, size(x)
         f = f + this%a * x(i)**4
         g(i) = 4.d0 * this%a * x(i)**3
      end do
   end subroutine myFunction_t_eval

   function myfctn_eval(this, x) result(f)
      class(myFunction_t), intent(in) :: this
      real(r8), intent(in) :: x(:)       ! coordinate
      real(r8)             :: f          ! function value
      integer i
      f = 0.d0
      do i = 1, size(x)
         f = f + this%a * x(i)**4
      end do
   end function myfctn_eval


   !!! six-dimensional function with 1s-like singularity at r_A = (0,0,-1)
   !!! and gaussian behaviour at r_B = (0,0,1). Symmetric wrt particle interchange
   !!! Squared antisymmetric function
   !!! f(r_1, r_2) = - ln h^2 with h =  g(r_1)*h(r_2) - h(r_1)*g(r_2)
   !!! g(r) = exp(-a*(r - rA)) and h(r) =  exp(-b*(r - rB)^2)

   subroutine mysingFunction_t_set(this, a, b)
      class(mysingFunction_t), intent(inout) :: this
      real(r8), intent(in) :: a, b
      call assert(a > 0 .and. b > 0, "mysingfctn: a > 0 and b > 0 required")
      this%a = a
      this%b = b
   end subroutine mysingFunction_t_set

   subroutine mysingFunction_t_get(this, a, b)
      class(mysingFunction_t), intent(in) :: this
      real(r8), intent(inout) :: a, b
      a = this%a
      b = this%b
   end subroutine mysingFunction_t_get

   subroutine mysingFunction_t_eval(this, x, f, g)
      class(mysingFunction_t), intent(in) :: this
      real(r8), intent(in) :: x(:)       ! coordinate
      real(r8), intent(out) :: f       ! function value
      real(r8), intent(out) :: g(:)    ! gradient
      real(r8) g1d(3), g2d(3), h1d(3), h2d(3), f1d(3), f2d(3)
      real(r8) ra1, ra2, rb12, rb22, g1, g2, h1, h2
      if (asserts) call assert(size(x) == size(g) .and. size(x) == 6, "mysingfctn: size must be 6")
      f = 0.d0

      ra1 = sqrt(x(1)**2 + x(2)**2 + (x(3) - this%za)**2)
      rb12 = x(1)**2 + x(2)**2 + (x(3) - this%zb)**2
      ra2 = sqrt(x(4)**2 + x(5)**2 + (x(6) - this%za)**2)
      rb22 = x(4)**2 + x(5)**2 + (x(6) - this%zb)**2

      g1 = exp(-this%a * ra1)
      g2 = exp(-this%a * ra2)
      h1 = exp(-this%b * rb12)
      h2 = exp(-this%b * rb22)

      g1d = - this%a * g1 * ( x(1:3) - [0.d0, 0.d0, this%za] ) / ra1
      g2d = - this%a * g2 * ( x(4:6) - [0.d0, 0.d0, this%za] ) / ra2
      h1d = - 2.d0 * this%b * h1 * ( x(1:3) - [0.d0, 0.d0, this%zb] )
      h2d = - 2.d0 * this%b * h2 * ( x(4:6) - [0.d0, 0.d0, this%zb] )

      f = g1 * h2 - h1 * g2
      f1d = g1d * h2 - h1d * g2
      f2d = g1 * h2d - h1 * g2d

      g(1:3) = - 2.d0 / f * f1d
      g(4:6) = - 2.d0 / f * f2d

      f = - log(f**2)

   end subroutine mysingFunction_t_eval

   function mysingfctn_eval(this, x) result(res)
      class(mysingFunction_t), intent(in) :: this
      real(r8), intent(in) :: x(:)       ! coordinate
      real(r8)             :: res        ! function value
      real(r8) ra1, ra2, rb1, rb2, g1, g2, h1, h2
      if (asserts) call assert(size(x) == 6, "mysingfctn: size must be 6")

      ra1 = sqrt(x(1)**2 + x(2)**2 + (x(3) - this%za)**2)
      rb1 = sqrt(x(1)**2 + x(2)**2 + (x(3) - this%zb)**2)
      ra2 = sqrt(x(4)**2 + x(5)**2 + (x(6) - this%za)**2)
      rb2 = sqrt(x(4)**2 + x(5)**2 + (x(6) - this%zb)**2)

      g1 = exp(-this%a * ra1)
      g2 = exp(-this%a * ra2)
      h1 = exp(-this%b * rb1**2)
      h2 = exp(-this%b * rb2**2)

      res = - 2.d0 * log ( abs( g1 * h2 - h1 * g2 ) )

   end function mysingfctn_eval


end module myfctn_module