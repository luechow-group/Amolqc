! Copyright (C) 1999, 2016 Arne Luechow
! Copyright (C) 2015 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later


! module for data required for CI optimization

module multiDetParam_m
  use kinds_m, only: r8
  use error_m
  use global_m, only: ne
  use wfData_m, only: nalpha,nbeta,ndet,ncsf
  use multiDet_m, only: cci,ndets,ccsf,deter,dgrad,dlapli
  implicit none
  public

  ! derivatives dphi / dp_k where p_k is the k-th parameter
  integer                     :: ci_params=0
  integer                     :: ci_param_mode=1
  real(r8), allocatable, target :: fk(:)
  real(r8), allocatable, target :: fkgrad(:,:)
  real(r8), allocatable, target :: fklapl(:)
  real(r8), allocatable, target :: fklapli(:,:)

CONTAINS

  !-------------------------------!
  subroutine mdetparam_create(mode)
  !-------------------------------!
  integer, intent(in) :: mode   ! mode: which ci params to optimize
  integer np

  ci_param_mode = mode
  if (ci_param_mode == 3) then
     ci_params = ncsf
  else
     ci_params = ncsf - 1
  endif
  np = ci_params

  if (allocated(fklapli)) then
     if (size(fklapli,1)/=ne .or. size(fklapli,2)/=np) then
        call mdetparam_destroy()
        allocate(fk(np),fkgrad(3*ne,np),fklapl(np),fklapli(ne,np))
     endif
  else
     allocate(fk(np),fkgrad(3*ne,np),fklapl(np),fklapli(ne,np))
  endif
  end subroutine mdetparam_create


  !----------------------------!
  subroutine mdetparam_destroy()
  !----------------------------!
     deallocate(fk,fkgrad,fklapl,fklapli)
  end subroutine mdetparam_destroy

  !-------------------------------!
  integer function getCIParamsCnt()
  !-------------------------------!
     getCIParamsCnt = ci_params
  end function getCIParamsCnt



  subroutine mdetcalcderivs(doElecDerivs)
  !-------------------------------------!
     logical, optional, intent(in) :: doElecDerivs
     integer i,j,n,na,k,offset
     real(r8) tmp
     logical calcElecDerivs

     calcElecDerivs = .true.
     if (present(doElecDerivs)) then
        if (.not.doElecDerivs) calcElecDerivs = .false.
     end if

     n = 0
     na = 3*nalpha
     if (ci_param_mode == 3) then
        offset = 0
     else
        offset = 1
        n = n + ndets(1)
     endif

     do k=1,ncsf-offset
        fk(k) = 0.d0
        fkgrad(1:3*ne,k) = 0.d0
        fklapli(1:ne,k) = 0.d0
        do j=1,ndets(k+offset)
           n = n+1
           tmp = ccsf(j,k+offset)
           fk(k) = fk(k) + tmp*deter(1,n)*deter(2,n)
           if (calcElecDerivs) then 
              do i=1,na
                 fkgrad(i,k) = fkgrad(i,k) + tmp*dgrad(i,1,n)*deter(2,n)
              enddo
              do i=1,3*nbeta
                 fkgrad(na+i,k) = fkgrad(na+i,k) + tmp*deter(1,n)*dgrad(i,2,n)
              enddo
              do i=1,nalpha
                 fklapli(i,k) = fklapli(i,k) + tmp*dlapli(i,1,n)*deter(2,n)
              enddo
              do i=1,nbeta
                 fklapli(nalpha+i,k) = fklapli(nalpha+i,k) + tmp*deter(1,n)*dlapli(i,2,n)
              enddo
           end if
        enddo
        fklapl(k) = 0.d0
        do i=1,ne
           fklapl(k) = fklapl(k) + fklapli(i,k)
        enddo
     enddo
     call assert(n==ndet,'mdetcalcderivs: internal error')
  end subroutine mdetcalcderivs

  !-------------------------------------!
  subroutine putCIParamsVector(optMode,p,normCI)
  !-------------------------------------!
     integer, intent(in)         :: optMode       ! optimization mode
     real(r8), intent(in)          :: p(:)         ! parameter vector
     logical,optional,intent(in) :: normCI   ! normalize ci coeefs after optimization
     integer k,k0

     if (optMode==3) then
        k0 = 1
     else
        k0 = 2
     endif

     cci(k0:ncsf) = p(1:ncsf-k0+1)
     if(normCI) cci=cci/sqrt(dot_product(cci,cci))
  end subroutine putCIParamsVector

  !-------------------------------------!
  subroutine getCIParamsVector(optMode,p)
  !-------------------------------------!
     integer, intent(in)   :: optMode       ! optimization mode
     real(r8), intent(inout) :: p(:)         ! parameter vector
     integer k,k0

     if (optMode==3) then
        k0 = 1
     else
        k0 = 2
     endif
     call assert(size(p) >= ncsf-k0+1,"getCIParamsVector: illegal parameter vector length")

     p(1:ncsf-k0+1) = cci(k0:ncsf)
   end subroutine getCIParamsVector

end module multiDetParam_m
