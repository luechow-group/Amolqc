! Copyright (C) 2013-2014 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module refSimple_m

   ! this module stores reference electron configurations
   ! reference positions (e.g. local maxima of psi**2) can be added (append or insert)
   ! and are internally sorted according to different criteria
   ! as output, .ref-files are generated containing the sorted reference configurations
   ! this module creates an abstract base class (referenceContainer) and a "create" method
   ! that creates the actual class. Only the base class methods need to be used in applications
   ! (polymorphic objects)
   ! The basis data structure is a list of list of reference electron configurations, requiring
   ! two indices for access
   !
   ! parallelism: this container should be kept on the MASTER (only)
   !   collect data (GATHER) on the MASTER and use insert/append for all gathered data
   !   no combining/merging is done in this module! use local statistics only!

   use refBase_m

   implicit none

   public

   type, extends(referenceContainer) :: refContainerSimple
   contains
      procedure :: refcsi_create
      procedure :: insert => refcsi_insert
      procedure :: writeShortList => refcsi_writeShortList
   end type refContainerSimple


contains

   subroutine refcsi_create(this,listLength,doSortFreq)
      class(refContainerSimple), intent(inout) :: this
      integer, intent(in)                           :: listLength
      logical, intent(in)                           :: doSortFreq
      call refc_create(this,listLength,doSortFreq)
   end subroutine refcsi_create

   subroutine refcsi_insert(this,r)
      class(refContainerSimple), intent(inout) :: this
      type(reference), intent(in)              :: r
      call abortp("refcsi_insert: insert must not be called")
   end subroutine refcsi_insert


   subroutine refcsi_writeShortList(this,iu,str)
      class(refContainerSimple)             :: this
      integer, intent(in)                   :: iu     ! unit to write to
      character(len=*), intent(in)          :: str    ! output string "Ref" or "Max"
      call abortp("internal error: refcsi_writeShortList must not be called")
   end subroutine refcsi_writeShortList

end module refSimple_m
