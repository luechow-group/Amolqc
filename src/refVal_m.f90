! Copyright (C) 2013-2014 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refVal_m

   ! "value"-type references, identify references by their function value

   use kinds_m, only: r8
   use refBase_m
 
   implicit none

   type, extends(referenceContainer) :: refContainerValue
      real(r8)            :: valueThreshold = 1.d-3       ! identify function values as equal if diff is smaller
   contains
      procedure :: refcv_create
      procedure :: insert => refcv_insert
      procedure :: writeShortList => refcv_writeShortList
   end type refContainerValue


contains


   subroutine refcv_create(this,listLength,doSortFreq,valueThreshold)
      class(refContainerValue), intent(inout) :: this
      integer, intent(in)                           :: listLength
      logical, intent(in)                           :: doSortFreq
      real(r8), optional, intent(in)                  :: valueThreshold
      call refc_create(this,listLength,doSortFreq)
      if (present(valueThreshold)) this%valueThreshold = valueThreshold
   end subroutine refcv_create

   subroutine refcv_insert(this,r)
      class(refContainerValue), intent(inout) :: this
      type(reference), intent(in)              :: r
      type(reference), pointer                 :: rp
      type(refl_vlist)                         :: rl
      type(refl_vlist), pointer                :: p
      integer i
      logical found

      call r%setCount(1)
      call rl%create(this%getElemLength()+1)
      found = .false.
      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         if (abs(r%f-rp%f) < this%valueThreshold) then
            found = .true.
            ! function value is already in list, but structure might be still different!
            if (p%size() < this%getElemLength()) then
               call p%append(r)
            else
               call rp%incrementCount()    ! discard but increase count of FIRST sub list entry
            end if
            exit
         else if (r%f < rp%f) then
            ! new function value found and inserted sortedly
            found = .true.
            call rl%del()
            call rl%append(r)
            call this%insertList(i,rl)
            exit
         end if
      end do
      if (.not.found .and. this%getSize() < this%getMaxListLength()) then
         ! new last entry
         call rl%del()
         call rl%append(r)
         call this%appendList(rl)
      end if
      if (this%getSize() > this%getMaxListLength()) call this%delLastList()
      call rl%destroy()
   end subroutine refcv_insert

   subroutine refcv_writeShortList(this,iu,str)
      class(refContainerValue)             :: this
      integer, intent(in)                   :: iu     ! unit to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist), pointer             :: p
      integer i,totalCount
      real(r8) f

      write(iu,'(a/)') ' list sorted with respect to function values -ln(psi**2):'
      write(iu,'(3x,a5,a13,4x,2a10)') trim(str),'value','# found'

      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         ! sum counters of all sub list elements
         totalCount = rp%count
         f = rp%f
         write(iu,'(5x,i4,2x,f13.6,2x,i7)') i,f,totalCount
      end do
   end subroutine refcv_writeShortList


end module refVal_m
