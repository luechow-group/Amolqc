! Copyright (C) 2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module partition_m
   use kinds_m, only: r8
   use error_m, only: assert
   implicit none

   private
   public :: partition_t, assignment (=)

   type partition_t
      ! electronsPerFragment is of size MAXVAL(fragments)+1. The last element are electrons, that are not in any fragment
      integer, allocatable :: electronsPerFragment(:)
      integer :: count = 0
   contains
      procedure :: init => partition_init
      procedure :: write => partition_write
   end type partition_t

   interface assignment (=)
      module procedure partition_copy
   end interface

   interface partition_t
      module procedure partition_constructor
   end interface partition_t

contains
   subroutine partition_init(this, electronsPerFragment, count)
      class(partition_t), intent(inout) :: this
      integer, intent(in) :: electronsPerFragment(:)
      integer, intent(in) :: count

      if (allocated(this%electronsPerFragment)) deallocate(this%electronsPerFragment)

      ! using allocation on assign
      this%electronsPerFragment = electronsPerFragment

      this%count = count
   end subroutine partition_init

   function partition_constructor(electronsPerFragment, count) result(partition)
      integer, intent(in) :: electronsPerFragment(:)
      integer, intent(in) :: count
      type(partition_t) :: partition

      call partition%init(electronsPerFragment, count)
   end function partition_constructor

   subroutine partition_copy(this, other)
      class(partition_t), intent(inout) :: this
      class(partition_t), intent(in) :: other

      call this%init(other%electronsPerFragment, other%count)
   end subroutine partition_copy

   subroutine partition_write(this, iul, nSamples)
      class(partition_t), intent(in) :: this
      integer, intent(in) :: iul, nSamples

      write(iul,'(f11.8,100i4)') 1._r8 * this%count/nSamples, this%electronsPerFragment
   end subroutine partition_write

end module partition_m