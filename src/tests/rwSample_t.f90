! Copyright (C) 2009 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

program test

  use kinds_m, only: r8
  use init_m, only: initGen
  use eloc_m
  use elocData_m
  use waveFunction_m
  use randomWalker_m
  use rwSample_m

  implicit none
  
  type(RandomWalker), pointer :: rwPtr
  type(RandomWalker), pointer :: rwa(:)
  type(RandomWalker), pointer :: rwaPtr(:)
  type(RandomWalker)          :: rw
  type(RWSample)              :: sample
  integer n,i,bs
  real(r8)                      :: x(4),y(4),z(4)
  real(r8), pointer             :: xp(:),yp(:),zp(:)
  real(r8), pointer             :: xdp(:),ydp(:),zdp(:)
  character(len=120)          :: lines(1)
  type(eConfigArray)          :: ec

  ! initialize wf
  lines(1) = "verbose=0   "
  call initGen(lines,1)
  lines(1) = "read,file='test.wf'"
  call readWriteWF(lines,1)

  call assert(ne==4," ne=4 (Be) required")
  n = ne
  call setNumberOfElectrons(n)
  x(1) = 0.5; y(1) = 0  ; z(1) = 0
  x(2) = 0.0; y(2) = 1.6; z(2) = 0
  x(3) = 0.0; y(3) = 0  ; z(3) = -0.6
  x(4) = 1.5; y(4) = 0  ; z(4) = 0

  call eConfigArray_new(ec,ne,1)
  call eConfigArray_set(ec,1,x,y,z)
  call eloc(0,ec,'none')

  call rw_new(rw)
  call setTo(rw,x,y,z)
  !call display(rw,2)
  call assert(nElecs()==4,'Test rw 1')
  print*, 'Test rw 1 passed'
  call assert(sizeOfRandomWalker()==31,'Test rw 2')
  print*, 'Test rw 2 passed'
  call assert(nDWGen()==1,'Test rw 3')
  print*, 'Test rw 3 passed'

  x(1) = 0.59
  allocate(rwPtr)
  call rw_new(rwPtr)
  call setTo(rwPtr,x,y,z)
  !call display(rwPtr,2)

  allocate(rwa(5))
  z(2) = 0.1
  call rw_new(rwa(1))
  call setTo(rwa(1),x,y,z)
  z(2) = 0.2
  call rw_new(rwa(2))
  call setTo(rwa(2),x,y,z)
  z(2) = 0.3
  call rw_new(rwa(3))
  call setTo(rwa(3),x,y,z)
  call pos(rwa(2),x,y,z)
  call assertEqualAbsolute(z(2),0.2d0,1d-5,'Testa 1')
  print*, 'Test a 1 passed'
  call posPtr(rwa(1),xp,yp,zp)
  call assertEqualAbsolute(zp(2),0.1d0,1d-5,'Testa 2')
  print*, 'Test a 2 passed'

  rwaPtr => rwa(1:3)
  call assert(size(rwaPtr)==3,'Testa 3')
  print*, 'Test a 3 passed'
  call pos(rwaPtr(2),x,y,z)
  call assertEqualAbsolute(z(2),0.2d0,1d-5,'Testa 4')
  print*, 'Test a 4 passed'
  call posPtr(rwaPtr(1),xp,yp,zp)
  call assertEqualAbsolute(zp(2),0.1d0,1d-5,'Testa 5')
  print*, 'Test a 5 passed'


  call initSample(sample,5)
  call appendWalker(sample,rw)
  call appendWalker(sample,rwPtr)
  call appendWalker(sample,rwa(1))
  call appendWalker(sample,rwa(3))
  call assert(getSampleSize(sample)==4,'Test 1')
  print*, 'Test 1 passed'
  call assert(getMaxSampleSize(sample)==5,'Test 2')
  print*, 'Test 2 passed'

  !call displaySample(sample)
  call rw_destroy(rwPtr)
  deallocate(rwPtr)
  !print*, "after destroy"
  !call displaySample(sample)

  rwPtr => getFirst(sample)
  call pos(rwPtr,x,y,z)
  call assertEqualAbsolute(x(1),0.5d0,1.d-5,'Test 3')
  print*, 'Test 3 passed'
  rwPtr => getNext(sample)
  call pos(rwPtr,x,y,z)
  call assertEqualAbsolute(x(1),0.59d0,1.d-5,'Test 4')
  print*, 'Test 4 passed'
  rwPtr => getNext(sample)
  rwPtr => getNext(sample)
  call assert(isLast(sample),'Test 5')
  print*, 'Test 5 passed'

  call posPtr(rwPtr,xp,yp,zp)
  call driftPtr(rwPtr,xdp,ydp,zdp)
  x(1:n) = xp + 0.01*xdp
  y(1:n) = yp + 0.01*ydp
  z(1:n) = zp + 0.01*zdp
  rwPtr => getFirst(sample)
  call setTo(rwPtr,x,y,z)
  !print*, "after move"
  !call displaySample(sample)
  !call deleteSample(sample)

  bs = 3
  call getFirstWalkerBlock(sample,bs,rwaPtr)
  call assert(size(rwaPtr)==3,'Test block 1')
  print*, 'Test block 1 passed'
  call pos(rwaPtr(3),x,y,z)
  call assertEqualAbsolute(z(2),0.1d0,1.d-5,'Test block 2')
  print*, 'Test block 2 passed'
  call getNextWalkerBlock(sample,bs,rwaPtr)
  call assert(bs==1,'Test block 3')
  print*, 'Test block 3 passed'
  call assert(size(rwaPtr)==1,'Test block 4')
  print*, 'Test block 4 passed'
  call pos(rwaPtr(1),x,y,z)
  call assertEqualAbsolute(z(2),0.3d0,1.d-5,'Test block 5')
  print*, 'Test block 5 passed'

  print*, 'no assertation failed'

end program test
  
