! Copyright (C) 2018 Arne Luechow
! Copyright (C) 2020 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module verbosity_m
   implicit none

   private
   public :: Verbosity_t

   type :: Verbosity_t
      private
      integer :: verbose_ = 0
      integer :: verbose_unit_ = 6
   contains
      procedure :: verbose
      procedure :: set_verbose
      procedure :: verbose_unit
      procedure :: set_verbose_unit
   end type Verbosity_t

contains
   function verbose(this) result (res)
      class(Verbosity_t), intent(in) :: this
      integer                      :: res
      res = this%verbose_
   end function verbose

   subroutine set_verbose(this, val)
      class(Verbosity_t), intent(inout) :: this
      integer, intent(in)             :: val
      this%verbose_ = val
   end subroutine set_verbose

   function verbose_unit(this) result (res)
      class(Verbosity_t), intent(in) :: this
      integer                      :: res
      res = this%verbose_unit_
   end function verbose_unit

   subroutine set_verbose_unit(this, val)
      class(Verbosity_t), intent(inout) :: this
      integer, intent(in)             :: val
      this%verbose_unit_ = val
   end subroutine set_verbose_unit

end module verbosity_m