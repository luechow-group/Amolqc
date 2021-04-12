! Copyright (C) 2018 Arne Luechow
! Copyright (C) 2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module maximizeSample_m
   use kinds_m, only: r8
   use parsing_m, only: getinta, getloga
   use global_m
   use mpiInterface_m, only: myMPIReduceSumInteger
   use rWSample_m, only: rWSample, getSampleSize, getFirst, isNext, getNext, getWalker
   use psiMax_m, only: psimax
   use randomWalker_m, only: RandomWalker

   implicit none

   private
   public :: maximizeSample, maximizeWalker

contains
   subroutine maximizeSample(smpl, psimax_obj, wout)
      type(RWSample), intent(inout) :: smpl
      type(psimax), intent(inout) :: psimax_obj
      logical, intent(in) :: wout
      type(RandomWalker), pointer :: rwp
      integer :: n, nc, neval, niter, na(4), na_all(4)
      n = getSampleSize(smpl)
      nc = 0
      neval = 0
      niter = 0
      rwp => getFirst(smpl)
      do
         call psimax_obj%opt(rwp, update_walker=.true.)
         if (psimax_obj%is_converged()) then
            nc = nc + 1
            neval = neval + psimax_obj%function_evaluations()
            niter = niter + psimax_obj%iterations()
         end if
         if (.not. isNext(smpl)) exit
         rwp => getNext(smpl)
      enddo
      na(1) = n; na(2) = nc; na(3) = neval; na(4) = niter
      na_all = na
      call myMPIReduceSumInteger(na, na_all, 4)
      if (wout) then
         write(iul,'(a,2(i6,a))') " maximization of ", na_all(1), " electron configurations with", na_all(2), " converged"
         if (na_all(2) > 0) then
            write(iul,'(a,2(i6,a))') " mean # of function evaluations: ", na_all(3) / na_all(2)
            write(iul,'(a,2(i6,a))') " mean # of iterations          : ", na_all(4) / na_all(2)
         end if
      end if
      call psimax_obj%write_results()
   end subroutine maximizeSample

   subroutine maximizeWalker(lines, nl, smpl, psimax_obj, wout)
      integer, intent(in) :: nl
      character(len=120), intent(in) :: lines(nl)
      type(RWSample), intent(inout) :: smpl
      type(psimax), intent(inout) :: psimax_obj
      logical, intent(in) :: wout
      type(RandomWalker), pointer :: rwp
      integer :: walkerIndex, iflag
      logical :: update_walker

      call getinta(lines,nl,'index=',walkerIndex,iflag)

      call getloga(lines,nl,'update_walker=',update_walker,iflag)
      if (iflag /= 0) update_walker = .true.

      rwp => getWalker(smpl,walkerIndex)
      call psimax_obj%opt(rwp, update_walker=update_walker)
      if (wout) then
         if (psimax_obj%is_converged()) then
            write(iul,'(a,i6,a)') " maximization of electron configuration #", walkerIndex, " did converge"
            write(iul,'(a,2(i6,a))') "# of function evaluations: ", psimax_obj%function_evaluations()
            write(iul,'(a,2(i6,a))') "# of iterations          : ", psimax_obj%iterations()
         else
            write(iul,'(a,i6,a)') " maximization of electron configuration #", walkerIndex, " did not converge"
         end if
      end if
   end subroutine maximizeWalker

end module maximizeSample_m
