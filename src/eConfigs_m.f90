! Copyright (C) 2012 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

! "class" defining array of electron configurations
! 

! cConfigsModule defines types (as "classes")
! for a electron configuration (x,y,z position of n electrons)
! and a type containing an array of electron configurations.
! Access allows pointer access to data for efficiency
! Construction only by copying for savety

! TODO: Delete unused copy routines?
!
module eConfigs_m
   use kinds_m, only: r8
   use error_m
   implicit none

   private
   public :: EConfigArray, eConfigArray_get, eConfigArray_size, eConfigArray_set, eConfigArray_new,&
           eConfigArray_destroy, eConfigArray_getPtr
   
   type EConfig
      real(r8), pointer :: x(:) => null()
      real(r8), pointer :: y(:) => null()
      real(r8), pointer :: z(:) => null()
   end type

   type EConfigArray
      type(EConfig), pointer :: eca(:) => null()
      integer :: nEConfigs = 0
   contains
      procedure :: new => eConfigArray_new
      procedure :: destroy => eConfigArray_destroy
      procedure :: get => eConfigArray_get
      procedure :: set => eConfigArray_set
      procedure :: getPtr => eConfigArray_getPtr
      procedure :: getSize => eConfigArray_size
      procedure :: getNElec => eConfigArray_nElec
   end type

   integer :: mNElecs = 0

contains

   ! methods of EConfig

   subroutine eConfig_new(self,nel)
      type(EConfig), intent(inout) :: self
      integer, intent(in)          :: nel
      integer alstat
      call assert((mNElecs==0 .and. nel>0) .or. nel==mNElecs,'eConfig_new: illegal argument')
      allocate(self%x(nel),self%y(nel),self%z(nel),stat=alstat)
      call assert(alstat==0,'eConfig_new: allocation failed')
      self%x = 0; self%y = 0; self%z = 0
      if (mNElecs==0) mNElecs = nel
   end subroutine eConfig_new

   subroutine eConfig_destroy(self)
      type(EConfig), intent(inout)  :: self
      integer alstat
      deallocate(self%x,self%y,self%z,stat=alstat)
      call assert(alstat==0,'eConfig_destroy: deallocation failed')
   end subroutine eConfig_destroy

   logical function eConfig_isAllocated(self)
      type(EConfig), intent(in) :: self
      eConfig_isAllocated = associated(self%x)
   end function eConfig_isAllocated

   integer function eConfig_size(self)
      type(EConfig), intent(in) :: self
      eConfig_size = size(self%x)
   end function eConfig_size

   subroutine eConfig_set(self,x,y,z)
      type(EConfig), intent(inout) :: self
      real(r8), intent(in)           :: x(:)
      real(r8), intent(in)           :: y(:)
      real(r8), intent(in)           :: z(:)
      call assert(size(x)>=mNElecs,'eConfig_set: illegal argument size')
      !! COPY coordinates
      self%x = x(1:mNElecs)
      self%y = y(1:mNElecs)
      self%z = z(1:mNElecs)
   end subroutine eConfig_set

   subroutine eConfig_get(self,x,y,z)
      type(EConfig), intent(in)       :: self
      real(r8), intent(inout)           :: x(:)
      real(r8), intent(inout)           :: y(:)
      real(r8), intent(inout)           :: z(:)
      call assert(size(x)>=mNElecs,'eConfig_set: illegal argument size')
      !! COPY coordinates
      x(1:mNElecs) = self%x
      y(1:mNElecs) = self%y
      z(1:mNElecs) = self%z
   end subroutine eConfig_get

   subroutine eConfig_getptr(self,x,y,z)
      type(EConfig), intent(inout)   :: self
      real(r8), pointer                :: x(:)
      real(r8), pointer                :: y(:)
      real(r8), pointer                :: z(:)
      !! return pointer to arrays
      x => self%x
      y => self%y
      z => self%z
   end subroutine eConfig_getptr

   subroutine eConfig_copy(self,copy)
      type(EConfig), intent(in)       :: self
      type(EConfig), intent(inout)    :: copy
      call assert(associated(self%x),'eConfig_copy: copy arg not associated')
      if (.not.associated(copy%x)) then
         call eConfig_new(copy,size(self%x))
      else if (size(self%x) /= size(copy%x)) then
         call eConfig_destroy(copy)
         call eConfig_new(copy,size(self%x))
      end if
      copy%x = self%x
      copy%y = self%y
      copy%z = self%z
   end subroutine eConfig_copy

   !!! assignment and addition may be added

   !!! methods for eConfigArray

   subroutine eConfigArray_new(self,nel,nec)
      class(EConfigArray),intent(inout)  :: self
      integer, intent(in)               :: nel ! # of electrons
      integer, intent(in)               :: nec ! # of electron configurations
      integer alstat,i
      call assert(nel>0 .and. nec>0,'eConfigArray_new: illegal argument values')
      allocate(self%eca(nec),stat=alstat)
      call assert(alstat==0,'eConfigArray_new: allocation failed (1)')
      do i=1,nec
         call eConfig_new(self%eca(i),nel)
      enddo
      self%nEConfigs = nec
   end subroutine eConfigArray_new

   subroutine eConfigArray_destroy(self)
      class(EConfigArray)  :: self
      integer i,alstat
      do i=1,self%nEConfigs
         call eConfig_destroy(self%eca(i))
      enddo
      deallocate(self%eca,stat=alstat)
      call assert(alstat==0,'eConfigArray_destroy: deallocation failed')
      self%nEConfigs = 0
   end subroutine eConfigArray_destroy

   subroutine eConfigArray_set(self,i,x,y,z)
      class(EConfigArray),intent(inout)  :: self
      integer, intent(in)          :: i
      real(r8), intent(in)           :: x(:)
      real(r8), intent(in)           :: y(:)
      real(r8), intent(in)           :: z(:)
      call assert(size(x)>=mNElecs .and. i>0 .and. i<= self%nEConfigs, &
                  'eConfigArray_set: illegal argument or argument size')
      call eConfig_set(self%eca(i),x,y,z)
   end subroutine eConfigArray_set

   subroutine eConfigArray_get(self,i,x,y,z)
      class(EConfigArray),intent(in)  :: self
      integer, intent(in)          :: i
      real(r8), intent(inout)        :: x(:)
      real(r8), intent(inout)        :: y(:)
      real(r8), intent(inout)        :: z(:)
      call assert(size(x)>=mNElecs .and. i>0 .and. i<= self%nEConfigs, &
                  'eConfigArray_get: illegal argument or argument size')
      call eConfig_get(self%eca(i),x,y,z)
   end subroutine eConfigArray_get

   subroutine eConfigArray_getPtr(self,i,x,y,z)
      class(EConfigArray),intent(inout)  :: self
      integer, intent(in)          :: i
      real(r8), pointer :: x(:)
      real(r8), pointer :: y(:)
      real(r8), pointer :: z(:)
      call assert(i>0 .and. i<= self%nEConfigs, &
                 'eConfigArray_getPtr: illegal argument')
      call eConfig_getPtr(self%eca(i),x,y,z)
   end subroutine eConfigArray_getPtr

   subroutine eConfigArray_copy(self,copy)
      class(EConfigArray),intent(in)    :: self
      type(EConfigArray),intent(inout) :: copy
      integer i
      call assert(self%nEConfigs>0.and.mNElecs>0,'eConfigArray_copy: copy arg not set')
      if (.not.associated(copy%eca)) then
         call eConfigArray_new(copy,mNElecs,self%nEConfigs)
      else if (copy%nEConfigs /= self%nEConfigs) then
         call eConfigArray_destroy(copy)
         call eConfigArray_new(copy,mNElecs,self%nEConfigs)
      end if
      do i=1,self%nEConfigs
         call eConfig_copy(self%eca(i),copy%eca(i))
      end do
   end subroutine eConfigArray_copy

   pure integer function eConfigArray_size(self)
      class(EConfigArray),intent(in)  :: self
      eConfigArray_size = self%nEConfigs
   end function eConfigArray_size

   pure integer function eConfigArray_nElec(self)
      class(EConfigArray),intent(in)  :: self
      eConfigArray_nElec = mNElecs
   end function eConfigArray_nElec


end module eConfigs_m
