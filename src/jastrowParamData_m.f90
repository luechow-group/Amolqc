! Copyright (C) 2009 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

! Jastrow module for keeping data e.g. for parameter derivatives for ALL Jastrow functions

MODULE jastrowParamData_m

   use kinds_m, only: r8
   use error_m
   implicit none
   public

   logical                     :: JPD_enforceECusp = .true.

   ! derivatives dU / dp_k where p_k is the k-th parameter
   integer                     :: uk_params=0
   real(r8), allocatable, target :: uk(:)
   real(r8), allocatable, target :: ukgrad(:,:)
   real(r8), allocatable, target :: uklapl(:)
   real(r8), allocatable, target :: uklapli(:,:)

CONTAINS

   subroutine jastrowparamdata_create(np,ne)
      integer, intent(in) :: np   ! # of jastrow parameters to optimize
      integer, intent(in) :: ne   ! # of electrons
      call assert(np>0 .and. ne>0, 'jastrowparamdata_create: positive arguments required')
      if (allocated(uklapli)) then
        if (size(uklapli,1)/=ne .or. size(uklapli,2)/=np) then
           call jastrowparamdata_destroy()
           allocate(uk(np),ukgrad(3*ne,np),uklapl(np),uklapli(ne,np))
           !if (use_ecp) call vNlk_create(np)
        endif
      else
        allocate(uk(np),ukgrad(3*ne,np),uklapl(np),uklapli(ne,np))
        !if (use_ecp) call vNlk_create(np)
      endif
      uk_params = np
   end subroutine jastrowparamdata_create

   subroutine jastrowparamdata_destroy()
      deallocate(uk,ukgrad,uklapl,uklapli)
   end subroutine jastrowparamdata_destroy


END MODULE jastrowParamData_m
