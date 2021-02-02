! Copyright (C) 2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

program test
   use kinds_m, only: r8
   use sphericalIntGrids_m
   implicit none

   integer n,r,idx0,nGridPoints
   real(r8) I
   real(r8) v(3),w,sumw

   do n=1,size(RuleIdx)-1
      idx0 = RuleIdx(n) 
      nGridPoints = RuleIdx(n+1) - idx0
      print*, n,nGridPoints,' rule points'
      I = 0.d0
      sumw = 0.d0
      do r=1,nGridPoints
         v = AllGridPoints(:,idx0+r)
         w = AllWeights(idx0+r)
         I = I + w*testfunc(v)
         sumw = sumw + w
         !!!print*,w,v(1)**2 + v(2)**2 + v(3)**2
      enddo
      print*, 'rule sumw int ',n,sumw,I
   enddo

contains

   real(r8) function testfunc(v)
      real(r8), intent(in) :: v(3)
      testfunc = -v(1)**2 + v(2)**2 - v(3)**2
   end function testfunc

end program test