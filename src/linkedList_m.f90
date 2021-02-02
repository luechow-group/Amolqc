! Copyright (C) 2011 Alexander Sturm
!
! SPDX-License-Identifier: GPL-3.0-or-later

module linkedList_m
  use error_m
  implicit none
  private
  public :: ListNode, List, insertNode, deallocateList

  type ListNode
    integer :: idx
    type(ListNode), pointer :: next => null()
  end type ListNode

  type List
    type(ListNode), pointer :: node => null(), last => null()
  end type List

contains

  subroutine insertNode(this, idx)
    type(List), pointer, intent(inout) :: this
    integer, intent(in) :: idx

    type(ListNode), pointer :: node
    integer :: ierror

    call assert(associated(this), "Tried to insert linked list node into null list")

    if(.not. associated(this%node)) then
      allocate(this%node, stat=ierror)
      node => this%node
    else
      allocate(this%last%next, stat=ierror)
      node => this%last%next
    endif

    if (ierror /= 0) then
      call error("linkedlist failed to allocate next node")
    endif

    this%last => node
    node%idx = idx
  end subroutine insertNode

  subroutine deallocateList(listArray)
    type(List), intent(inout) :: listArray(:, :, :)
    type(ListNode), pointer :: node, next

    integer :: i, j, k

    do i = lbound(listArray, 1), ubound(listArray, 1)
      do j = lbound(listArray, 2), ubound(listArray, 2)
        do k = lbound(listArray, 3), ubound(listArray, 3)
          node => listArray(i, j, k)%node
          do while(associated(node))
            next => node%next
            deallocate(node)
            node => next
          enddo
          !XXX set everything to null() for reuse?
        enddo
      enddo
    enddo
  end subroutine deallocateList

end module linkedList_m
