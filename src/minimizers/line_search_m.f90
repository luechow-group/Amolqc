! Copyright (C) 2018 Arne Luechow
! Copyright (C) 2020 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module line_search_m
   use kinds_m, only: r8
   use fctn_module, only: Function_t
   use singularityCorrection_m, only: singularity_correction
   use singularityParticles_m, only: singularity_particles
   use verbosity_m, only: Verbosity_t

   implicit none

   private
   public :: line_search

   type, extends(Verbosity_t), abstract :: line_search
      real(r8) :: delta_max_ = 0.d0   ! maximum step length (scaling search direction)
      real(r8) :: step_max_  = 0.d0   ! maximum step length per particle (modifying search direction)
   contains
      procedure(Find_step_i), deferred :: find_step
      procedure, non_overridable :: line_search_write_params
      procedure :: write_params => line_search_write_params
   end type line_search

   abstract interface
      subroutine Find_step_i(this, ff, x0, falpha, grad0, d, xalpha, gradalpha, nfeval, gmask)
         import :: line_search, Function_t, singularity_correction, singularity_particles, r8
         class(line_search), intent(in) :: this
         class(Function_t), intent(in) :: ff  ! function object

         real(r8), intent(in) :: x0(:)   ! position
         real(r8), intent(inout) :: falpha  ! function value (in: at x, out: at x_new)
         real(r8), intent(in) :: grad0(:)! gradient
         real(r8), intent(in) :: d(:)    ! search direction
         real(r8), intent(inout) :: xalpha(:)    ! new position
         real(r8), intent(inout) :: gradalpha(:) ! gradient at new position
         integer, intent(inout) :: nfeval  ! out: required function evaluations, or error code (<0)

         logical, intent(in), optional :: gmask(:)  ! mask to force gradient component to zero (at singularities)
      end subroutine Find_step_i
   end interface

contains

   subroutine line_search_write_params(this, iu)
   class(line_search), intent(in)    :: this
   integer, intent(in) , optional  :: iu
   integer iull

   iull = this%verbose_unit()
   if (present(iu)) iull = iu
   write(iull, "(a/a)") " -- line search", " with:"
   if (this%delta_max_ > 0._r8) write(iull, "(a,g13.5)") " max_distance        =", this%delta_max_
   if (this%step_max_ > 0._r8) write(iull, "(a,g13.5)") " max_step        =", this%step_max_
end subroutine line_search_write_params

end module line_search_m
