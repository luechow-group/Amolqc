! Copyright (C) 2013-2014 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refValStr_m

   ! "value/structure"-type references, identify as equal by value 

   use kinds_m, only: r8
   use refBase_m
   use refUtils_m, only: calcRefDifference

   implicit none

   integer, parameter :: SIMPLE=0, SPINPERM=1, ALLPERM=2     ! modes for determining reference differences

   public
   private :: SIMPLE,SPINPERM,ALLPERM

   type, extends(referenceContainer) :: refContainerValStr
      real(r8)            :: valueThreshold = 1.d-3       ! identify function values as equal if diff is smaller
      integer           :: distMode = SPINPERM
      real(r8)            :: sameThreshold = 1.d-2       ! identify structures as equal if diff is smaller
   contains
      procedure :: refcvs_create
      procedure :: insert => refcvs_insert
      procedure :: writeShortList => refcvs_writeShortList
   end type refContainerValStr


contains

   subroutine refcvs_create(this,listLength,doSortFreq,el,valueThreshold,distanceMode,sameThreshold)
      class(refContainerValStr), intent(inout) :: this
      integer, intent(in)                            :: listLength
      logical, intent(in)                            :: doSortFreq
      integer, intent(in)                            :: el
      real(r8), optional, intent(in)                   :: valueThreshold
      integer, optional, intent(in)                  :: distanceMode
      real(r8), optional, intent(in)                   :: sameThreshold
      call refc_create(this,listLength,doSortFreq,el)
      if (present(valueThreshold)) this%valueThreshold = valueThreshold
      if (present(distanceMode)) this%distMode = distanceMode
      if (present(sameThreshold)) this%sameThreshold = sameThreshold
   end subroutine refcvs_create

   subroutine refcvs_insert(this,r)
      class(refContainerValStr), intent(inout) :: this
      type(reference), intent(in)              :: r
      type(reference), pointer                 :: rp,rpmin
      type(refl_vlist)                         :: rl
      type(refl_vlist), pointer                :: p
      real(r8)                                   :: dist,minDist,maxDist
      integer i,j
      logical found

      call r%setCount(1)
      call rl%create(this%getElemLength()+1)
      !!write(999,*) '*** vs insert value= ',r%f,' size=',this%getSize()
      found = .false.
      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         !!write(999,*) '   ** ',rp%f
         if (abs(r%f-rp%f) < this%valueThreshold) then
            ! function value is already in list, but structure might be still different!
            ! search in sub list for same structure
            found = .true.
            !!write(999,*)'   ** found!'
            minDist = 1.d99
            rpmin => null()
            do j=1,p%size()
               rp => p%elem(j)
               call calcRefDifference(rp,r,this%distMode,dist,maxdist)
               !!write(999,*)'      * dist =',dist,this%distMode
               if (maxdist < minDist) then
                  minDist = maxdist
                  rpmin => rp
               end if
            end do
            if (minDist > this%sameThreshold) then
               ! new structure found
               !!write(999,*) '   ** new struct'
               if (p%size() < this%getElemLength()) then
                  call p%append(r)
                  !!write(999,*) '   ** append'
               else
                  call rp%incrementCount()   ! discard increase count of LAST sub list entry
                  !!write(999,*) '   ** inc count'
               end if
            else
               ! structure already in sub list at node rpmin
               call rpmin%incrementCount()
            end if
            exit
         else if (r%f < rp%f) then
            ! new function value found. insert sorted
            found = .true.
            !!write(999,*) '   ** new value: insert here'
            call rl%del()
            call rl%append(r)
            call this%insertList(i,rl)
            exit
         end if
      end do
      if (.not.found .and. this%getSize() < this%getMaxListLength()) then
         ! new last entry
         !!write(999,*) '   ** new value: append'
         call rl%del()
         call rl%append(r)
         call this%appendList(rl)
      end if
      if (this%getSize() > this%getMaxListLength()) call this%delLastList()
      call rl%destroy()
   end subroutine refcvs_insert

   subroutine refcvs_writeShortList(this,iu,str)
      class(refContainerValStr)             :: this
      integer, intent(in)                   :: iu     ! unit to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      type(reference), pointer              :: rp
      type(refl_vlist)                      :: rl
      type(refl_vlist), pointer             :: p
      integer i,j,io,totalCount
      real(r8) f

      write(iu,'(a)') ' list sorted with respect to function value -ln(psi**2)'
      write(iu,'(a/)') ' each list element is list of different structures'

      write(iu,'(3x,a5,a13,4x,2a10)') trim(str),'value','# found','# structs'

      do i=1,this%getSize()
         call this%getElemPtr(i,p)
         rp => p%elem(1)
         ! sum counters of all sub list elements
         totalCount = rp%count
         f = rp%f
         do j=2,p%size()
            rp => p%elem(j)
            totalCount = totalCount + rp%count
         end do
         write(iu,'(5x,i4,2x,f13.6,2(2x,i7))') i,f,totalCount,p%size()
      end do
   end subroutine refcvs_writeShortList


end module refValStr_m
