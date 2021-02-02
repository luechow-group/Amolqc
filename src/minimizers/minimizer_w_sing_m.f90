! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module minimizer_w_sing_module

   ! this data structure is suitable only for collections of 3d particles (dim=3n)
   ! with non-smooth minima (singularities) at certain particle positions 
   ! (array of 3d positions)
   ! By non-smooth minima we mean discontinuous derivatives in 3 of 3n gradient
   ! components.
   ! The type allows calculation of distances and vectors to the nearest 3d (singular)
   ! position as well as setting particles at a singular position.
   ! The position data structure is assumed to be r = (x_1, y_1, z_1, x_2, ... z_n)

   use kinds_m, only: r8
   use minimizer_module, only: minimizer
   use error_m, only: asserts, assert
   use singularityCorrection_m, only: singularity_correction
   implicit none

   private
   public :: minimizer_w_sing

   type, abstract, extends(minimizer) :: minimizer_w_sing
      type(singularity_correction)    :: sc_
   contains
      procedure                               :: set_singularities
      procedure                               :: singularities
      procedure, non_overridable              :: write_params_minimizer_w_sing
      procedure                               :: write_params => write_params_minimizer_w_sing
   end type minimizer_w_sing 

contains
   subroutine set_singularities(this, s, last_a, scalings)
      class(minimizer_w_sing), intent(inout) :: this
      real(r8), intent(in) :: s(:,:)
      integer, intent(in), optional :: last_a
      real(r8), intent(in), optional :: scalings(:)
      call this%sc_%set_singularities(s)
      if (present(last_a)) call this%sc_%set_last_a(last_a)
      if (present(scalings)) then
         if (asserts) call assert(SIZE(s, dim=2) == SIZE(scalings), 'set_singularities: dimensions do not match')
         call this%sc_%set_scalings(scalings)
      end if
   end subroutine set_singularities

   function singularities(this) result (s)
      class(minimizer_w_sing), intent(in) :: this
      real(r8), allocatable :: s(:,:)
      allocate (s, source=this%sc_%singularities())
   end function singularities

   subroutine write_params_minimizer_w_sing(this, iu)
      class(minimizer_w_sing), intent(in) :: this
      integer, intent(in) , optional :: iu

      if (present(iu)) then
         call this%write_params_minimizer(iu) ! call parent class
         call this%sc_%write_params(iu)
      else
         call this%write_params_minimizer()
         call this%sc_%write_params()
      end if
   end subroutine write_params_minimizer_w_sing

end module minimizer_w_sing_module
