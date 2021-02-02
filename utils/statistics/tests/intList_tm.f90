! TODO UNKNOWN ORIGIN AND AUTHOR (OLDER THAN 2011?)
! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module intList_tm
  use intList_m
  use error_m
  implicit none
contains

  subroutine IntList_test
    integer i,j
    integer, pointer :: k,iPtr
    type(IntNode), pointer :: p
    type(IntList) :: l,m,l1,l2
    type(IntList), pointer :: listArray(:) => null()
    character(len=60) :: str

    !print*, "elementary tests"

    call createList(l)

    call assert(isEmpty(l),"A1")

    do i=5,1,-1
       call pushBack(l,i)
    enddo
    call assert(frontPtr(l)==5,"A2")
    call assert(backPtr(l)==1,"A3")
    call assert(.not.isEmpty(l),"A4")
    call assert(getSize(l)==5,"A5")

    k => frontPtr(l)
    call assert(k==5,"A6")

    k => backPtr(l)
    call assert(k==1,"A7")

    call listToString(l,str)
    call assert(str(1:15)=='  5  4  3  2  1',"A8")

    ! change value of 3nd element
    p => beginPtr(l)
    call nextPtr(p)
    call nextPtr(p)
    iPtr => dataPtr(p)
    iPtr = 7

    call listToString(l,str)
    call assert(str(1:15)=='  5  4  7  2  1',"A9")

    ! Einfuegen vor p
    p => beginPtr(l)
    call nextPtr(p)
    p => insert(l,p,8)

    call listToString(l,str)
    call assert(str(1:18)=='  5  8  4  7  2  1',"A10")
    !call printList(l)
    !print*, "... passed"

    !print*, "testing sort and mergeList"

    call sortIntList(l)

    call listToString(l,str)
    call assert(str(1:18)=='  1  2  4  5  7  8',"A11")

    call createList(m)
    do i=3,15,3
       call pushFront(m,i)
    enddo
    call assert(getSize(m)==5,"A12")
    call listToString(m,str)
    call assert(str(1:15)==' 15 12  9  6  3',"A12")

    !print*,"sorting decreasing list"
    call sortIntList(m)
    call listToString(m,str)
    call assert(str(1:15)=='  3  6  9 12 15',"A13")

    !print*,"mergeList test"
    !call printList(l)
    !call printList(m)
    call mergeList(l,m)
    call assert(getSize(l)==11,"A14")
    call assert(getSize(m)==0,"A15")
    call assert(isEmpty(m),"A16")
    call listToString(l,str)
    call assert(str(1:33)=='  1  2  3  4  5  6  7  8  9 12 15',"A17")

    call createList(l1)
    call assert(isEmpty(l1),"A17a")
    call createList(l2)
    do i=1,5
       call pushBack(l2,i)
    enddo
    call mergeList(l1,l2)
    call assert(getSize(l1)==5,"A17b")
    call assert(getSize(l2)==0,"A17c")
    call assert(isEmpty(l2),"A17d")

    !print* ,"... passed"

    !print*,"splice and sort test"
    p => beginPtr(l)
    do i=1,5
       p => erase(l,p)
       call nextPtr(p)
    enddo
    call assert(getSize(l)==6,"A18")
    call listToString(l,str)
    call assert(str(1:18)=='  2  4  6  8 12 15',"A19")

    do i=21,25
       call pushFront(m,i)
    enddo
    call listToString(m,str)
    call assert(str(1:15)==' 25 24 23 22 21',"A20")

    ! insert contents of l into m after first element.
    p => beginPtr(m)
    call nextPtr(p)
    call splice(m,p,l)
    call assert(getSize(l)==0,"A21")
    call assert(getSize(m)==11,"A22")

    !call printList(m)
    call listToString(m,str)
    call assert(str(1:33)==' 25  2  4  6  8 12 15 24 23 22 21',"A23")

    call sortIntList(m)
    !call printList(m)
    call listToString(m,str)
    call assert(str(1:33)=='  2  4  6  8 12 15 21 22 23 24 25',"A24")

    !print*,"... passed"
    !print*," destroying lists"
    call destroyList(l)
    call destroyList(m)
    !print*," ... done"

    !print*," testing arrays of lists"
    allocate(listArray(5))
    do i=1,5
       call createList(listArray(i))
    enddo

    do i=1,5
       do j=3,1,-1
          call pushBack(listArray(i),10*i+j)
       enddo
    enddo
    call listToString(listArray(1),str)
    call assert(str(1:9)==' 13 12 11',"A25")
    call listToString(listArray(5),str)
    call assert(str(1:9)==' 53 52 51',"A26")
    call sortIntList(listArray(3))
    call listToString(listArray(3),str)
    call assert(str(1:9)==' 31 32 33',"A27")

    !print*, "all tests passed"

  end subroutine IntList_test

  subroutine printList(list)
  type(IntList) :: list
  type(IntNode), pointer :: p
  p => beginPtr(list)
  do
     if (isEnd(list,p)) exit
     write(*,'(i3)',advance='no') dataPtr(p)
     call nextPtr(p)
  enddo
  write(*,*)
  end subroutine printList

  subroutine listToString(list,s)
  type(IntList) :: list
  character(len=*) :: s
  type(IntNode), pointer :: p
  integer i
  i=1; p => beginPtr(list)
  do
     if (isEnd(list,p)) exit
     write(s(i:i+2),'(i3)') dataPtr(p)
     i=i+3; call nextPtr(p)
  enddo
  end subroutine listToString

end module intList_tm
