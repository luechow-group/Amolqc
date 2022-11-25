! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module singularityParticles_m
   use kinds_m, only: r8
   use error_m, only: assert, asserts
   implicit none

   private
   public :: singularity_particles, assignment(=)

   interface assignment(=)
      module procedure singularity_particles_assign
   end interface

   type :: singularity_particles
      integer, allocatable :: slist(:)  ! contains the singularity each particle is in (0 for not in singularity)
      integer, allocatable :: nuc_list(:)  ! contains the particles that are in each singularity
      real(r8) :: sing_diag_H = 1._r8  ! sets the diagonal element of the Hessian for particles in the singularities
                                       ! for optimization, this should be 1. For EVAnalysis some large number.
   contains
      procedure :: create => singularity_particles_create
      procedure :: n_sing => singularity_particles_n_sing
      procedure :: At_singularity => singularity_particles_At_singularity  ! creates a mask of size 3*SIZE(this%slist)
      procedure :: Fix_gradients => singularity_particles_Fix_gradients  ! sets gradients to zero and add unit matrix
                                                                         ! blocks to Hessian
   end type singularity_particles
contains

   subroutine singularity_particles_create(this, n_particles, n_singularities, sing_diag_H)
      class(singularity_particles), intent(inout) :: this
      integer, intent(in) :: n_particles, n_singularities
      real(r8), intent(in), optional :: sing_diag_H
      allocate(this%slist(n_particles), this%nuc_list(n_singularities))
      this%slist = 0
      this%nuc_list = 0
      if (PRESENT(sing_diag_H)) this%sing_diag_H = sing_diag_H
   end subroutine singularity_particles_create

   function singularity_particles_n_sing(this) result (res)
      class(singularity_particles), intent(in)  :: this
      integer                            :: res
      ! alpha: nuc_list(i) = 3
      ! beta: nuc_list(i) = 4
      ! both: nuc_list(i) = 7
      res = COUNT(this%nuc_list /= 0) + COUNT(this%nuc_list == 7)
   end function singularity_particles_n_sing

   subroutine singularity_particles_assign(copy, src)
      class(singularity_particles), intent(inout) :: copy
      class(singularity_particles), intent(in)    :: src
      if (.not.allocated(copy%slist)) then
         allocate(copy%slist(size(src%slist)), copy%nuc_list(size(src%nuc_list)))
      end if
      copy%slist = src%slist
      copy%nuc_list = src%nuc_list
   end subroutine singularity_particles_assign

   function singularity_particles_At_singularity(this) result(mask)
      class(singularity_particles), intent(in) :: this
      logical :: mask(3 * SIZE(this%slist))
      integer :: i

      mask = [ (this%slist(i/3+1) /= 0, i=0, 3*SIZE(this%slist)-1) ]
   end function singularity_particles_At_singularity

   subroutine singularity_particles_Fix_gradients(this, g, H)
      class(singularity_particles), intent(in) :: this
      real(r8), intent(inout) :: g(3*SIZE(this%slist))
      real(r8), intent(inout), optional :: H(3*SIZE(this%slist),3*SIZE(this%slist))

      logical :: mask(3 * SIZE(this%slist))
      integer :: i

      mask = this%At_singularity()

      ! setting gradient to zero for singularity particles
      where (mask) g = 0._r8

      if (PRESENT(H)) then
         ! because we need to invert the hessian, the singularity coordinates have to be set to zero
         ! except for the diagonal elements, which should be one
         ! this way, H is blocked and inversion is inverting each block on its own
         do i=1, SIZE(mask)
            if (mask(i)) then
               H(i,:) = 0._r8
               H(:,i) = 0._r8
               H(i,i) = this%sing_diag_H
            end if
         end do
      end if
   end subroutine singularity_particles_Fix_gradients

end module singularityParticles_m