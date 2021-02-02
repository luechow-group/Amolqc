! Copyright (C) 2012, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refADT_m

   use kinds_m, only: r8
   use error_m
   implicit none

   type reference
      real(r8), pointer :: x(:) => null()
      real(r8), pointer :: y(:) => null()
      real(r8), pointer :: z(:) => null()
      real(r8)          :: f = 0
      integer         :: count = 0
   contains
      procedure :: new => ref_new
      procedure :: destroy => ref_destroy
      !!procedure :: copy => ref_assign
      procedure :: set => ref_set
      procedure :: get => ref_get
      procedure :: getCount => ref_getCount
      procedure :: setCount => ref_setCount
      procedure :: incrementCount => ref_incrCount
      procedure :: writeRef => ref_write
   end type

   interface assignment(=)
      module procedure ref_assign
   end interface

contains

   subroutine ref_new(this,n)
      class(reference)    :: this
      integer, intent(in) :: n
      call assert(.not.associated(this%x),'reference type new: already allocated')
      allocate(this%x(n),this%y(n),this%z(n))
   end subroutine ref_new
 
   subroutine ref_destroy(this)
      class(reference)    :: this
      deallocate(this%x,this%y,this%z)
   end subroutine ref_destroy

   elemental subroutine ref_assign(this,rhs)
      type(reference), intent(out) :: this
      type(reference), intent(in) :: rhs
      if (associated(this%x) .and. size(this%x)/=size(rhs%x)) deallocate(this%x,this%y,this%z)
      if (.not.associated(this%x)) allocate(this%x(size(rhs%x)),this%y(size(rhs%y)),this%z(size(rhs%z)))
      this%x = rhs%x 
      this%y = rhs%y 
      this%z = rhs%z 
      this%f = rhs%f 
      this%count = rhs%count
   end subroutine ref_assign

   subroutine ref_get(this,x,y,z,f)
      class(reference)    :: this
      real(r8), intent(inout) :: x(:)
      real(r8), intent(inout) :: y(:)
      real(r8), intent(inout) :: z(:)
      real(r8), intent(inout) :: f
      call assert(associated(this%x).and.size(this%x)==size(x),'reference type get: not allocated or size error')
      x = this%x 
      y = this%y 
      z = this%z 
      f = this%f 
   end subroutine ref_get

   subroutine ref_set(this,x,y,z,f)
      class(reference)    :: this
      real(r8), intent(in) :: x(:)
      real(r8), intent(in) :: y(:)
      real(r8), intent(in) :: z(:)
      real(r8), intent(in) :: f
      call assert(associated(this%x).and.size(this%x)==size(x),'reference type set: not allocated or size error')
      this%x = x
      this%y = y
      this%z = z
      this%f = f
   end subroutine ref_set

   integer function ref_getCount(this)
      class(reference) :: this
      ref_getcount = this%count
   end function ref_getCount

   subroutine ref_setCount(this,c)
      class(reference)    :: this
      integer, intent(in) :: c
      this%count = c
   end subroutine ref_setCount
      
   subroutine ref_incrCount(this)
      class(reference) :: this
      this%count = this%count + 1
   end subroutine ref_incrCount

   subroutine ref_write(this,iu)
      class(reference), intent(in) :: this
      integer, intent(in)          :: iu
      integer i
      write(iu,'(g15.5,i6)') this%f,this%count
      write(iu,'(i5)') size(this%x)
      do i=1,size(this%x)
         write(iu,'(i5,3g15.5)') i,this%x(i),this%y(i),this%z(i)
      end do
   end subroutine ref_write
      

end module
