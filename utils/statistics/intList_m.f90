! Copyright (C) ca. 2000 Arne Luechow 
!
! SPDX-License-Identifier: GPL-3.0-or-later


! attempt to construct abstract container types similar to C++ STL
! Due to the lack of templates in F90/F95 individual instantiations
! of containers can be done automatically with 'replace' using an editor
! (such as sed)
!
! Like STL containers: no access to components, but to a pointer
! to the elements (Iterator in STL)
!
! not complete: unique, reverse are missing

MODULE intList_m

  use error_m

  implicit none

  type IntNode
     private
     integer                        :: data = 0
     type(IntNode), pointer         :: next => null()
     type(IntNode), pointer         :: previous => null()
  end type IntNode

  type IntList
     private
     integer                :: listSize=0
     type(IntNode), pointer :: root    => null()  ! before first element
     type(IntNode), pointer :: lend    => null()  ! after last element
  end type IntList

  public
  private :: simpleSort

CONTAINS

!-------- constructor and destructor -----------------------------


  !-------------------------!
  subroutine createList(s)
  !-------------------------!
    type(IntList), intent(inout)     :: s
    type(IntNode), pointer           :: rootNode
    type(IntNode), pointer           :: endNode

    call assert(.not.associated(s%root),"createList: list already initialized")
    allocate(rootNode,endNode)
    s%root => rootNode        ! root element. data not to be used
    s%lend => endNode          ! end element. data not to be used
    s%lend%previous => s%root
    s%root%next => s%lend
    s%listSize = 0
  end subroutine createList

  !--------------------------!
  subroutine destroyList(s)
  !--------------------------!
    type(IntList), intent(inout)  :: s

    call clear(s)
    deallocate(s%root,s%lend)
  end subroutine destroyList

  !-------------------------------!
  logical function isCreatedList(s)
  !-------------------------------!
    type(IntList), intent(in) :: s
    isCreatedList = associated(s%root)
  end function isCreatedList


!-------- iterators (access to nodes) ---------------------------

  !------------------!
  function beginPtr(s)
  !------------------!
    type(IntNode), pointer    :: beginPtr
    type(IntList), intent(in) :: s
    beginPtr => s%root%next
  end function beginPtr

  !----------------!
  function endPtr(s)
  !----------------!
    ! note: this is *after* the end (allowing inserting at the end!)
    type(IntNode), pointer    :: endPtr
    type(IntList), intent(in) :: s
    endPtr => s%lend
  end function endPtr

  !-------------------------!
  logical function isEnd(s,p)
  !-------------------------!
    ! note: this is *after* the end
    type(IntList), intent(in) :: s
    type(IntNode), pointer    :: p
    isEnd = associated(p,s%lend)
  end function isEnd

  !---------------------------!
  logical function isBegin(s,p)
  !---------------------------!
    type(IntList), intent(in) :: s
    type(IntNode), pointer    :: p
    isBegin = associated(p,s%root%next)
  end function isBegin

  !-------------------!
  subroutine nextPtr(p)
  !-------------------!
    type(IntNode), pointer    :: p
    call assert(associated(p%next),"nextPtr: failed")
    p => p%next
  end subroutine nextPtr

  !-----------------------!
  subroutine previousPtr(p)
  !-----------------------!
    type(IntNode), pointer    :: p
    call assert(associated(p%previous),"previousPtr: failed")
    p => p%previous
  end subroutine previousPtr


!-------- properties ------------------------------------

  !-------------------------!
  logical function isEmpty(s)
  !-------------------------!
    type(IntList), intent(in) :: s
    isEmpty = (s%listSize == 0)
    ! associated(s%root%next,s%lend)
  end function isEmpty

  !-------------------------!
  integer function getSize(s)
  !-------------------------!
    type(IntList), intent(in) :: s
    getSize = s%listSize
  end function getSize

!-------- access to elements ----------------------------

  !-----------------!
  function dataPtr(p)
  !-----------------!
    ! note: returning pointer to data element (allow change of list element)
    integer, pointer    :: dataPtr
    type(IntNode), pointer :: p
    dataPtr => p%data
  end function dataPtr

  !------------------!
  function frontPtr(s)
  !------------------!
    ! note: returning pointer to data element (allow change of list element)
    integer, pointer             :: frontPtr
    type(IntList), intent(inout) :: s
    frontPtr => s%root%next%data
  end function frontPtr

  !-----------------!
  function backPtr(s)
  !-----------------!
    ! note: returning pointer to data element (allow change of list element)
    integer, pointer             :: backPtr
    type(IntList), intent(inout) :: s
    backPtr => s%lend%previous%data
  end function backPtr

!--------- list changing operations -----------------------------

  !--------------------!
  function insert(s,p,t)
  !--------------------!
    ! insert a *copy* of t into list before position p, returning pointer to new node
    type(IntNode), pointer       :: insert
    type(IntList), intent(inout) :: s
    type(IntNode), pointer       :: p
    integer, intent(in)          :: t
    type(IntNode), pointer       :: newNode
    allocate(newNode)
    newNode%data = t
    newNode%next => p
    newNode%previous => p%previous
    p%previous%next => newNode
    p%previous => newNode
    s%listSize = s%listSize + 1
    insert => newNode
  end function insert

  !-------------------------!
  subroutine insertN(s,p,n,t)
  !-------------------------!
    ! insert n *copies* of t into list before position p
    type(IntList), intent(inout) :: s
    type(IntNode), pointer       :: p
    integer, intent(in)          :: n
    integer, intent(in)          :: t
    type(IntNode), pointer       :: newNode
    integer i
    do i=1,n
       allocate(newNode)
       newNode%data = t
       newNode%next => p
       newNode%previous => p%previous
       p%previous%next => newNode
       p%previous => newNode
       newNode => null()
    enddo
    s%listSize = s%listSize + n
  end subroutine insertN

  !-----------------------------!
  subroutine insertRange(s,p,i,j)
  !-----------------------------!
    ! insert a *copies* of elements from i to j (excluding j) into list
    ! before position p
    type(IntList), intent(inout) :: s
    type(IntNode), pointer       :: p
    type(IntNode), pointer       :: i
    type(IntNode), pointer       :: j
    type(IntNode), pointer       :: newNode
    do
       if (.not.associated(i,j)) exit
       allocate(newNode)
       newNode%data = i%data
       newNode%next => p
       newNode%previous => p%previous
       p%previous%next => newNode
       p%previous => newNode
       newNode => null()
       call nextPtr(i)
       s%listSize = s%listSize + 1
    enddo
  end subroutine insertRange

  !-----------------!
  function erase(s,p)
  !-----------------!
    ! erases element p points to, returning pointer node following p (might be end)
    type(IntNode), pointer       :: erase
    type(IntList), intent(inout) :: s
    type(IntNode), pointer       :: p
    call assert(associated(p%next).and.associated(p%previous),"erase: illegal element")
    erase => p%next
    p%previous%next => p%next
    p%next%previous => p%previous
    deallocate(p)
    p => null()
    s%listSize = s%listSize - 1
  end function erase

  !-----------------!
  subroutine clear(s)
  !-----------------!
    ! erases all element p from list s
    type(IntList), intent(inout) :: s
    type(IntNode), pointer       :: p
    p => s%root%next
    do
       if (isEnd(s,p)) exit
       p => erase(s,p)
    enddo
  end subroutine clear

  !-----------------------!
  subroutine pushFront(s,t)
  !-----------------------!
    ! == insert(beginPtr(s),t), insert t before first element of list
    type(IntList), intent(inout) :: s
    integer, intent(in)          :: t
    type(IntNode), pointer       :: p
    p => beginPtr(s)
    p => insert(s,p,t)
  end subroutine pushFront

  !--------------------!
  subroutine popFront(s)
  !--------------------!
    ! == erase(beginPtr(s),t), erases first element of list
    type(IntList), intent(inout) :: s
    type(IntNode), pointer       :: p
    call assert(.not.isEmpty(s),"popFront: erasing from empty list")
    p => beginPtr(s)
    p => erase(s,p)
  end subroutine popFront

  !----------------------!
  subroutine pushBack(s,t)
  !----------------------!
    ! == insert(endPtr(s),t), insert t after last element of list
    type(IntList), intent(inout) :: s
    integer, intent(in)          :: t
    type(IntNode), pointer       :: p
    p => endPtr(s)
    p => insert(s,p,t)
  end subroutine pushBack

  !-------------------!
  subroutine popBack(s)
  !-------------------!
    ! == erase(--endPtr(s),t), erases last element of list
    type(IntList), intent(inout) :: s
    type(IntNode), pointer       :: p
    call assert(.not.isEmpty(s),"popBack: erasing from empty list")
    p => endPtr(s)
    p => erase(s,p%previous)
  end subroutine popBack

  !----------------------!
  subroutine splice(s,p,l)
  !----------------------!
    ! insert (no copy) contents of list l before element p (of s). l is empty on return
    type(IntList), intent(inout) :: s
    type(IntNode), pointer       :: p
    type(IntList), intent(inout) :: l

    call assert(.not.associated(s%root,l%root),"splice: different lists required")
    p%previous%next => l%root%next
    l%root%next%previous => p%previous
    l%lend%previous%next => p
    p%previous => l%lend%previous
    s%listSize = s%listSize + l%listSize
    l%lend%previous => l%root
    l%root%next => l%lend
    l%listSize = 0
  end subroutine splice

  !----------------------!
  subroutine simpleSort(s)
  !----------------------!
    ! sort elements
    ! this is a selection sort implementation. Not to be used for large
    ! data sets (O(n**2))
    ! note: this implementation swaps data rather than list nodes
    type(IntList), intent(inout) :: s
    type(IntNode), pointer :: p,q,qmin
    integer, pointer :: m,n
    integer tmp
    call assert(.not. isEmpty(s),'sort: empty list')
    p => beginPtr(s)
    m => dataPtr(p)
    do
       q => p
       qmin => q
       n => m
       do
          call nextPtr(q)
          if (isEnd(s,q)) exit
          if (dataPtr(q) < n) then
             qmin => q
             n => dataPtr(qmin)
          endif
       enddo
       if (n<m) then
          tmp = n
          n = m
          m = tmp
       endif
       call nextPtr(p)
       if (isEnd(s,p)) exit
       m => dataPtr(p)
    enddo
  end subroutine simpleSort

  !-----------------------------!
  recursive subroutine sortIntList(list)
  !-----------------------------!
    type(IntList), intent(inout) :: list
    type(IntList)                :: right
    type(IntNode), pointer       :: middle,middle1,p
    integer                      :: i,n,l,r
    integer, parameter           :: maxSimpleSortSize = 2

    if (getSize(list) <= maxSimpleSortSize) then
       call simpleSort(list)
    else
       ! split list into left and right part. Left remains in list!
       call createList(right)
       n = getSize(list)
       l = getSize(list)/2
       r = n - l
       p => beginPtr(list)
       do i=1,l-1
          call nextPtr(p)
       end do
       middle => p
       middle1 => middle%next
       ! splitting list into left and right list
       right%root%next => middle1
       middle1%previous => right%root
       right%lend%previous => list%lend%previous
       right%lend%previous%next => right%lend
       list%lend%previous => middle
       middle%next => list%lend
       list%listSize = l
       right%listSize = r
       call sortIntList(list)
       call sortIntList(right)
       call mergeList(list,right)
       call destroyList(right)
    end if
  end subroutine sortIntList

  !------------------------------!
  subroutine mergeList(list,list1)
  !------------------------------!
    ! merges s1 into s. Assumes both lists to be sorted
    ! on exit s1 is empty
    type(IntList), intent(inout) :: list,list1
    type(IntNode), pointer :: p,q,p1
    p => beginPtr(list)
    q => beginPtr(list1)
    do
       if (isEnd(list1,q)) exit
       if (isEnd(list,p)) then
          do
             call pushBack(list,dataPtr(q))
             q => erase(list1,q)
             if (isEnd(list1,q)) exit
          enddo
          exit
       endif
       if (dataPtr(p) < dataPtr(q)) then
          call nextPtr(p)
       else
          p1 => insert(list,p,dataPtr(q))
          q => erase(list1,q)
       endif
    enddo
  end subroutine

end MODULE intList_m


