! Copyright (C) 2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module maximizeSampleRho_m
   use kinds_m, only: r8
   use parsing_m
   use global_m
   use rwSample_m, only: rWSample, getSampleSize, getFirst, isNext, getNext, getWalker
   use rhoMax_m, only: rhoMax_t
   use randomWalker_m, only: RandomWalker

   implicit none

   private
   public :: maximizeSampleRho, maximizeWalkerRho

contains
   subroutine maximizeSampleRho(smpl, rhoMax)
      type(RWSample), intent(inout) :: smpl
      type(rhoMax_t), intent(inout) :: rhoMax
      type(RandomWalker), pointer :: rwp

      rwp => getFirst(smpl)
      if (rhoMax%isInitialized()) then
         do
            call rhoMax%opt(rwp, update_walker=.true.)
            if (.not. isNext(smpl)) exit
            rwp => getNext(smpl)
         enddo
         call rhoMax%writeResults(iul)
      end if
   end subroutine maximizeSampleRho

   subroutine maximizeWalkerRho(lines, nl, smpl, rhoMax)
      integer, intent(in) :: nl
      character(len=120), intent(in) :: lines(nl)
      type(RWSample), intent(inout) :: smpl
      type(rhoMax_t), intent(inout) :: rhoMax
      type(RandomWalker), pointer :: rwp
      integer :: walkerIndex, iflag

      call getinta(lines,nl,'index=',walkerIndex,iflag)

      if (rhoMax%isInitialized()) then
         rwp => getWalker(smpl,walkerIndex)
         call rhoMax%opt(rwp, update_walker=.true.)
         call rhoMax%writeResults(iul)
      end if
   end subroutine maximizeWalkerRho
end module maximizeSampleRho_m