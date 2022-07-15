! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_module

   use kinds_m, only: r8
   use globalUtils_m, only: getMyTaskId
   use fctn_module, only: Function_t
   use verbosity_m, only: Verbosity_t
   implicit none

   private
   public :: minimizer

   type, extends(Verbosity_t), abstract :: minimizer
      !!!private
      real(r8)                        :: convergence_max_distance_ = 0.d0
      real(r8)                        :: convergence_max_gradient_ = 0.d0
      real(r8)                        :: convergence_max_gradnorm_ = 0.d0
      real(r8)                        :: convergence_value_ = -HUGE(0._r8)
      real(r8)                        :: value_ = 0.d0
      real(r8), allocatable           :: gradient_(:)
      real(r8)                        :: max_abs_gradient_ = 0.1d0
      real(r8)                        :: max_mean_gradient_ = 0.1d0
      real(r8)                        :: max_electron_distance_(3) = HUGE(1._r8)
      integer, allocatable            :: not_to_minimize_(:)
      logical                       :: converged_ = .false.
      logical                       :: value_convergence_ = .false.
      integer                       :: current_iterations_ = 0
      integer                       :: max_iterations_ = 100
      integer                       :: current_fctn_eval_ = 0
      integer                       :: path_count_ = 0
   contains
      procedure(minimizer_interface), deferred :: minimize
      procedure                                :: reset       ! resets current values only
      procedure                                :: is_converged
      procedure                                :: set_converged
      procedure                                :: set_convergence_distance
      procedure                                :: set_convergence_gradient
      procedure                                :: set_convergence_gradnorm
      procedure                                :: set_convergence_value
      procedure                                :: set_max_electron_distance
      procedure                                :: set_not_to_minimize
      procedure                                :: set_value_convergence
      procedure                                :: convergence_distance
      procedure                                :: value_convergence
      procedure                                :: convergence_gradient
      procedure                                :: convergence_gradnorm
      procedure                                :: convergence_value
      procedure                                :: max_electron_distance
      procedure                                :: is_distance_converged
      procedure                                :: is_gradient_converged
      procedure                                :: is_gradnorm_converged
      procedure                                :: is_value_converged
      procedure                                :: value
      procedure                                :: set_value
      procedure                                :: gradient
      procedure                                :: gradnorm
      procedure                                :: set_gradient
      procedure                                :: restrict_gradient
      procedure                                :: iterations
      procedure                                :: set_iterations
      procedure                                :: max_iterations
      procedure                                :: set_max_iterations
      procedure                                :: function_evaluations
      procedure                                :: set_function_evaluations
      procedure, non_overridable               :: write_params_minimizer
      procedure                                :: write_params => write_params_minimizer
      procedure                                :: write_opt_path
      procedure                                :: do_write_opt_path
      procedure                                :: write_opt_path_entry
   end type minimizer 

   abstract interface
      subroutine minimizer_interface(this, fn, x)
         import :: minimizer, Function_t, r8
         class(minimizer), intent(inout)  :: this
         class(Function_t), intent(in)     :: fn
         real(r8), intent(inout)            :: x(:)
      end subroutine minimizer_interface
   end interface

contains

   subroutine reset(this)
      class(minimizer)      :: this
      this%converged_ = .false.
      this%value_ = 0.d0
      if (allocated(this%gradient_)) this%gradient_ = 0.d0
      this%current_iterations_ = 0
      this%current_fctn_eval_ = 0
   end subroutine reset

   function is_converged(this) result (res)
      class(minimizer), intent(in) :: this
      logical                      :: res
      res = this%converged_
   end function is_converged

   subroutine set_converged(this, tf)
      class(minimizer), intent(inout) :: this
      logical, intent(in)             :: tf
      this%converged_ = tf
   end subroutine set_converged

   subroutine set_convergence_distance(this, distance)
      class(minimizer), intent(inout) :: this
      real(r8), intent(in)              :: distance
      this%convergence_max_distance_ = distance
   end subroutine set_convergence_distance

   subroutine set_convergence_gradient(this, gradient)
      class(minimizer), intent(inout) :: this
      real(r8), intent(in)              :: gradient
      this%convergence_max_gradient_ = gradient
   end subroutine set_convergence_gradient

   subroutine set_value_convergence(this, value_convergence)
      class(minimizer), intent(inout) :: this
      logical, intent(in)             :: value_convergence
      this%value_convergence_ = value_convergence
   end subroutine set_value_convergence

   subroutine set_convergence_gradnorm(this, gradnorm)
      class(minimizer), intent(inout) :: this
      real(r8), intent(in)              :: gradnorm
      this%convergence_max_gradnorm_ = gradnorm
   end subroutine set_convergence_gradnorm

   subroutine set_convergence_value(this, value)
      class(minimizer), intent(inout) :: this
      real(r8), intent(in)              :: value
      this%convergence_value_ = value
   end subroutine set_convergence_value

   subroutine set_max_electron_distance(this, value)
      class(minimizer), intent(inout) :: this
      real(r8), intent(in)              :: value(3)
      this%max_electron_distance_ = value
   end subroutine set_max_electron_distance

   subroutine set_not_to_minimize(this, not_to_minimize)
      class(minimizer), intent(inout) :: this
      integer, intent(in)             :: not_to_minimize(:)
      this%not_to_minimize_ = not_to_minimize
   end subroutine set_not_to_minimize

   subroutine restrict_gradient(this, gradient, hessian)
      class(minimizer), intent(inout)   :: this
      real(r8), intent(inout)           :: gradient(:)
      real(r8), intent(inout), optional :: hessian(:,:)
      integer                           :: i, j

      if (ALLOCATED(this%not_to_minimize_)) then
         do i=1, SIZE(this%not_to_minimize_)
            j = this%not_to_minimize_(i)
            gradient(j) = 0._r8
            if (PRESENT(hessian)) then
               hessian(j,:) = 0._r8
               hessian(:,j) = 0._r8
               hessian(j,j) = 1._r8
            end if
         end do
      end if
   end subroutine restrict_gradient

   function convergence_distance(this) result(distance)
      class(minimizer), intent(inout) :: this
      real(r8)                          :: distance
      distance = this%convergence_max_distance_ 
   end function convergence_distance

   function convergence_gradient(this) result(gradient)
      class(minimizer), intent(inout) :: this
      real(r8)                          :: gradient
      gradient = this%convergence_max_gradient_
   end function convergence_gradient

   function convergence_gradnorm(this) result(gradnorm)
      class(minimizer), intent(inout) :: this
      real(r8)                          :: gradnorm
      gradnorm = this%convergence_max_gradnorm_ 
   end function convergence_gradnorm

   function convergence_value(this) result(value)
      class(minimizer), intent(inout) :: this
      real(r8)                          :: value
      value = this%convergence_value_
   end function convergence_value

   function max_electron_distance(this) result(value)
      class(minimizer), intent(inout) :: this
      real(r8)                          :: value(3)
      value = this%max_electron_distance_
   end function max_electron_distance

   function value_convergence(this) result(value)
      class(minimizer), intent(inout) :: this
      logical                         :: value
      value = this%value_convergence_
   end function value_convergence

   function is_distance_converged(this, distance) result (res)
      class(minimizer), intent(in) :: this
      real(r8), intent(in)           :: distance
      logical                      :: res
      res = distance < this%convergence_max_distance_
   end function is_distance_converged

   function is_gradient_converged(this, gradient) result (res)
      class(minimizer), intent(in) :: this
      real(r8), intent(in)           :: gradient
      logical                      :: res
      res = gradient < this%convergence_max_gradient_
   end function is_gradient_converged

   function is_gradnorm_converged(this, gradnorm) result (res)
      class(minimizer), intent(in) :: this
      real(r8), intent(in)           :: gradnorm
      logical                      :: res
      res = gradnorm < this%convergence_max_gradnorm_
   end function is_gradnorm_converged

   function is_value_converged(this, value) result (res)
      class(minimizer), intent(in) :: this
      real(r8), intent(in)           :: value
      logical                      :: res
      res = value < this%convergence_value_
   end function is_value_converged

   function value(this) result (res)
      class(minimizer), intent(in) :: this
      real(r8)                       :: res
      res = this%value_
   end function value

   subroutine set_value(this, val)
      class(minimizer), intent(inout) :: this
      real(r8), intent(in)              :: val
      this%value_ = val
   end subroutine set_value

   function gradient(this) result (res)
      class(minimizer), intent(in) :: this
      real(r8), allocatable          :: res(:)
      if (allocated(this%gradient_)) then
         res = this%gradient_
      else
         allocate(res(1))
         res(1) = 0.d0
      end if
   end function gradient

   function gradnorm(this) result (res)
      class(minimizer), intent(in) :: this
      real(r8)                       :: res
      res = sqrt( dot_product( this%gradient_, this%gradient_ ) )
   end function gradnorm

   subroutine set_gradient(this, grad)
      class(minimizer), intent(inout) :: this
      real(r8), intent(in)              :: grad(:)
      this%gradient_ = grad               ! automatic allocation 
   end subroutine set_gradient

   function iterations(this) result (res)
      class(minimizer), intent(in) :: this
      integer                      :: res
      res = this%current_iterations_
   end function iterations

   subroutine set_iterations(this, val)
      class(minimizer), intent(inout) :: this
      integer, intent(in)             :: val
      this%current_iterations_ = val 
   end subroutine set_iterations

   function max_iterations(this) result (res)
      class(minimizer), intent(in) :: this
      integer                      :: res
      res = this%max_iterations_
   end function max_iterations

   subroutine set_max_iterations(this, val)
      class(minimizer), intent(inout) :: this
      integer, intent(in)             :: val
      this%max_iterations_ = val
   end subroutine set_max_iterations

   function function_evaluations(this) result (res)
      class(minimizer), intent(in) :: this
      integer                      :: res
      res = this%current_fctn_eval_
   end function function_evaluations

   subroutine set_function_evaluations(this, val)
      class(minimizer), intent(inout) :: this
      integer, intent(in)             :: val
      this%current_fctn_eval_ = val 
   end subroutine set_function_evaluations

   subroutine write_params_minimizer(this, iu)
      class(minimizer), intent(in) :: this
      integer, intent(in) , optional :: iu
      integer iull
      if (present(iu)) then
         iull = iu
      else
         iull = this%verbose_unit()
      end if

      write(iull,"(a)") " minimizer parameters:"
      write(iull,"(a,i8)") " max_iterations      =", this%max_iterations_
      write(iull,"(a,i3)") " verbose             =", this%verbose()
      if (this%convergence_max_distance_ > 0.d0)  &
         write(iull, "(a,g13.5)") " convergence_distance=", this%convergence_max_distance_
      if (this%convergence_max_gradient_ > 0.d0)  &
         write(iull, "(a,g13.5)") " convergence_gradient=", this%convergence_max_gradient_
      if (this%convergence_max_gradnorm_ > 0.d0)  &
         write(iull, "(a,g13.5)") " convergence_gradnorm=", this%convergence_max_gradnorm_
      if (ALLOCATED(this%not_to_minimize_))  &
         write(iull, "(a,g13.5)") " not_to_minimize=", this%not_to_minimize_
   end subroutine write_params_minimizer

   subroutine write_opt_path(this, counter)
      class(minimizer), intent(inout) :: this
      integer, intent(in)                    :: counter
      this%path_count_ = counter
   end subroutine write_opt_path

   function do_write_opt_path(this) result(res)
      class(minimizer), intent(inout) :: this
      logical                                :: res
      res = this%path_count_ > 0
   end function do_write_opt_path

   subroutine write_opt_path_entry(this, k, x, f)
      class(minimizer), intent(inout) :: this
      integer, intent(in)                    :: k
      real(r8), intent(in)                     :: x(:)
      real(r8), intent(in)                     :: f
      integer iu, i

      iu = 300 + getMyTaskId()
      write(iu,'(a, 2i6, a, f14.5, a)') "OPTPATH:", this%path_count_, k, "   F(Max):",  f, "  found: 0 0 0 0 0"
      write(iu,'(i5)') size(x) / 3
      do i = 1, size(x) / 3
         write(iu,'(3f14.7)') x(3*i-2), x(3*i-1), x(3*i)
      enddo
   end subroutine write_opt_path_entry

end module minimizer_module
