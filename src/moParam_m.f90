! Copyright (C) 2015-2016, 2018 Arne Luechow
! Copyright (C) 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

module moParam_m

   use kinds_m, only: r8
   use globalUtils_m, only: iul, iull, mytid, logmode, abortp
   use error_m
   use global_m, only: ne, MASTER
   use wfData_m, only: nalpha,nbeta,nbas,norb,ndet,ncsf,aomocomb
   use multiDet_m, only: cci,ndets,ccsf,deter,detia,detib,det1xa,det1xb,det1ya,det1yb,det1za,det1zb, &
                   det2a,det2b,mclist,lapack_inv,dgrad,dlapli,detsRepLst,detiaOne,detibOne
   use mos_m, only: cmo,mat,mat1x,mat1y,mat1z,mat2
   use parsing_m
   use linAlg_m, only: invdetcalc, invrupd

   implicit none
   private

   ! derivatives dphi / dp_k where p_k is the k-th parameter
   real(r8), allocatable, target :: mok(:)
   real(r8), allocatable, target :: mokgrad(:,:)
   real(r8), allocatable, target :: moklapl(:)
   real(r8), allocatable, target :: moklapli(:,:)
   real(r8), allocatable         :: cmo0(:,:)
   real(r8), allocatable         :: p_save(:)
   integer, allocatable        :: excitOp(:,:)      ! excitOp(1|2,p) yields rotation e1<->e2 for parameter p
   integer, allocatable        :: posInDeta(:,:),posInDetb(:,:)
   integer, allocatable        :: paramSymmetriseList(:,:)

   integer, parameter :: SIMPLE=1, UPDATE=2
   integer, parameter :: SUCCESSIVE=1, EXPKAPPA=2

   integer                     :: moParams = 0     ! # of orb rot parameters
   integer                     :: mUpdateMode = UPDATE
   integer                     :: mParamMode = EXPKAPPA
   integer                     :: mVerbose = 0     ! local verbosity
   integer                     :: mMaxPairs = 5    ! max number of orbital rotation pairs that are symmetrised
   integer                     :: mParamSymEntries = 0 ! length of list paramSymmetriseList

   public :: mok,mokgrad,moklapli,moklapl, moparam_create, moparam_destroy, moparam_calcderivs, &
             getMOParamsCnt, getMOParamsVector, putMOParamsVector, moparam_calcderivsOnlyMok, &
             getMoUpdateMode, getNParamSymEntries, getParamSymmetriseList, withSymmetriseList

contains

   subroutine moparam_create(mode,lines,nl,verbose)
   !----------------------------------------------!

   integer, intent(in)            :: mode           ! ignored ! check !
   integer, intent(in)            :: nl
   character(len=120), intent(in) :: lines(nl)
   integer, optional, intent(in) :: verbose         ! verbosity for moparam module
   integer, allocatable :: optArray(:,:)   ! array containing MOs to optimize (0 entry: length of row)
                                           ! format: optArray(0:,1) = [len, l1, l2, .. l_len]
                                           !         optArray(0:,2) = [len, m1, m2, .. m_len]
                                           !         optArray(0:,3) = |len, l'1, l'2, .. l'_len]
                                           !         optArray(0:,4) = |len, m'1, m'2, .. m'_len]
                                           ! meaning: all rotations \kappa_lm and all \kappa_l'm', ...
   integer, allocatable :: MOSymmetriseList(:,:)
   integer np,i,j1,j2,k,l,l0,m,nClass,len,p,orb,iflag, nSym, counter
   integer p1, e1, e2
   integer entryVec(0:mMaxPairs)
   logical found

   mVerbose = logmode
   if (present(verbose)) mVerbose = verbose

   mUpdateMode = UPDATE
   call getinta(lines,nl,"mo_update_mode=",mUpdateMode,iflag)

   mParamMode = EXPKAPPA
   call getinta(lines,nl,"mo_param_mode=",mParamMode,iflag)

   if (MASTER .and. logmode>=2) then
      write(iul,'(/a)') ' orbital rotation settings:'
      write(iul,'(3x,2(a15,i3)/)') ' mo_update_mode=',mUpdateMode,' mo_param_mode=',mParamMode
   endif

   do i=1,nl
      found = finda(lines(i),1,"orbital_rotation_list=")
      if (found) exit
   enddo
   if (i >= nl) call abortp("$optimize_parameters: params=mo requires 'orbital_rotation_list=' followed by formatted lines")
   l0 = i+1
   read(lines(l0),*) nClass
   if (MASTER .and. logmode >= 2) write(iul,'(a,i3,a)') ' reading ',nClass,' classes of orbital rotations'
   allocate(optArray(0:norb,2*nClass))
   optArray = 0
   m = 1
   counter = 0
   do l=1,2*nClass,2
      if (l0+l+1 > nl) call abortp("$optimize_parameters: illegal format in MO parameter description")
      m = l + counter
      read(lines(l0+m),*) len
      if (len > 30) then
         read(lines(l0+m),*) (optArray(i,l),i=0,30)
         read(lines(l0+m+1),*) (optArray(i,l),i=31,len)
         counter = counter + 1
      else
         read(lines(l0+m),*) (optArray(i,l),i=0,len)
      end if
      m = l + counter
      read(lines(l0+m+1),*) len
      if (len > 30) then
         read(lines(l0+m+1),*) (optArray(i,l+1),i=0,30)
         read(lines(l0+m+2),*) (optArray(i,l+1),i=31,len)
         counter = counter + 1
      else
         read(lines(l0+m+1),*) (optArray(i,l+1),i=0,len)
      end if
   enddo

   do i=1,nl
      found = finda(lines(i),1,"mo_symmetrise_list=")
      if (found) exit
   end do
   if (found) then
      ! read list of list of equivalent orbital rotations k <-> l
      ! format: n_pairs, k1, l1, k2, l2, ... kn_pairs, ln_pairs
      l0 = i+1
      read(lines(l0),*) nSym
      if (MASTER .and. logmode >= 2) write(iul,'(a,i3,a)') ' reading ',nSym,' lines of equivalent orbital rotations'
      allocate(MOSymmetriseList(0:mMaxPairs, nSym))
      MOSymmetriseList = 0
      do l=1,nSym
         read(lines(l0+l),*) len
         if (len > mMaxPairs) call abortp("reading mo_symmetrise_list: too many orbital rotation pairs")
         read(lines(l0+l),*) (MOSymmetriseList(i,l), i = 0, len)
      end do
   end if

   np = 0

   ! determine number of orbital rotation parameters from optArray
   do i=1,size(optArray,2),2
      do j1=1,optArray(0,i)
         do j2=1,optArray(0,i+1)
            if (optArray(j1,i) < optArray(j2,i+1)) then
               np = np + 1
            endif
         enddo
      enddo
   enddo

   moParams = np
   allocate(excitOp(2,np))

   ! generate list of orbital rotations
   if (MASTER .and. mVerbose >= 3) write(iul,'(/a)') ' list of orbital rotations:'
   np = 0
   do i=1,size(optArray,2),2
      do j1=1,optArray(0,i)
         do j2=1,optArray(0,i+1)
            if (optArray(j1,i) < optArray(j2,i+1)) then
               np = np + 1
               excitOp(1,np) = optArray(j1,i)
               excitOp(2,np) = optArray(j2,i+1)
               if (MASTER .and. mVerbose >= 3) write(iul,'(i5,3x,2i5)') np, excitOp(1,np), excitOp(2,np)
            endif
         enddo
      enddo
   enddo

   deallocate(optArray)

   ! generate list of equivalent rotation parameters that should symmetrised
   if (allocated(MOSymmetriseList)) then
      if (MASTER .and. mVerbose >= 3) write(iul,'(/a)') ' list of symmetry equivalent orbital rotation parameters:'
      allocate(paramSymmetriseList(0:mMaxPairs,np))
      mParamSymEntries = 0
      do p1 = 1, np
         e1 = excitOp(1,p1)
         e2 = excitOp(2,p1)
         !!!write(iul,*) "DBG:", p1, e1, e2
         call findEquivalentRots(MOSymmetriseList, p1, e1, e2, entryVec)
         !!!write(iul,'(a,10i3)') "DBG:", entryVec
         if (entryVec(0) > 0) then
            mParamSymEntries = mParamSymEntries + 1
            paramSymmetriseList(0:,mParamSymEntries) = entryVec
            if (MASTER .and. mVerbose >= 3) write(iul,'(i5,3x,20i5)') mParamSymEntries, entryVec(1:entryVec(0))
         end if
      end do
    deallocate(MOSymmetriseList)
   end if



   if (allocated(moklapli)) then
      if (size(moklapli,1)/=ne .or. size(moklapli,2)/=np) then
         deallocate(mok,mokgrad,moklapl,moklapli)
         allocate(mok(np),mokgrad(3*ne,np),moklapl(np),moklapli(ne,np))
      endif
   else
      allocate(mok(np),mokgrad(3*ne,np),moklapl(np),moklapli(ne,np))
   endif

   allocate(posInDeta(norb,ndet),posInDetb(norb,ndet))

   posInDeta = 0
   posInDetb = 0
   do k=1,ndet
      do i=1,nalpha
         orb = mclist(i,k)
         posInDeta(orb,k) = i
      enddo
      do i=1,nbeta
         orb = mclist(nalpha+i,k)
         posInDetb(orb,k) = i
      enddo
   enddo

   ! save cmo matrix to define
   if (aomocomb) call abortp("orbital optimization currently no possible with AOMO mode")
   allocate(cmo0(nbas,norb))
   cmo0 = cmo

   ! corresponding current parameter vector
   ! note: parameter vector is defined wrt the original MO coeff matrix cmo0
   allocate(p_save(np))
   p_save = 0.d0


   end subroutine moparam_create

   subroutine moparam_destroy()

   call assert(allocated(mok),"moparam_destroy called but array mok is not allocated")
   call assert(allocated(posInDeta),"moparam_destroy called but array posInDeta is not allocated")
   deallocate(mok,mokgrad,moklapl,moklapli)
   deallocate(excitOp,posInDeta,posInDetb)
   deallocate(cmo0)
   deallocate(p_save)
   if (allocated(paramSymmetriseList)) then
      deallocate(paramSymmetriseList)
      mParamSymEntries=0
   endif

   end subroutine moparam_destroy



   subroutine moparam_calcderivs(doElecDerivs)

   logical, optional, intent(in) :: doElecDerivs    ! if true calculate mokgrad/moklapl in addition to mok
   integer i,ii,j,jj,k,ll,mm,n,na,nb,p,lpos,mpos,orb,ierr,detssize
   real(r8) tmp,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,coeff,d,darr(size(detsRepLst,1))
   real(r8) ,allocatable ::  tmparr1(:,:),tmparr2(:,:),tmparr3(:,:),tmparr4(:,:)!keep results to avoid unnecessary calculations
   real(r8), allocatable :: excitedMO(:), tmpdet(:,:), tmpdeti(:,:), tmpdetx(:,:),tmpdety(:,:),tmpdetz(:,:),tmpdet2(:,:)
   integer, allocatable :: mclistexcit(:)
   logical calcElecDerivs

    detssize = size(detsRepLst,1)
   ! assert(nec==1)

   calcElecDerivs = .true.
   if (present(doElecDerivs)) then
      if (.not.doElecDerivs) calcElecDerivs = .false.
   end if

   mok = 0.d0
   mokgrad = 0.d0
   moklapl = 0.d0
   moklapli = 0.d0

   na = nalpha
   nb = nbeta

   ! do alpha electron part

   allocate(tmpdet(na,na),tmpdeti(na,na),tmpdetx(na,na),tmpdety(na,na),tmpdetz(na,na),tmpdet2(na,na))
   allocate(excitedMO(na),mclistexcit(na))
   allocate(tmparr1(na,detssize),tmparr2(na,detssize),tmparr3(na,detssize),tmparr4(na,detssize))
   do p=1,moParams   ! number of orb rot parameters: pairs of (l->m)
   n = 0
   ! calculate E_lm Phi
   ! loop over all CSF dets, count det number n
   do k=1,ncsf
      do j=1,ndets(k)
         n = n+1
         coeff = cci(k)*ccsf(j,k)

         ! loop over all orb rot parameters

            ! first do l -> m excitation: E_lm = a^T_m a_l
            ll   = excitOp(1,p)    ! MO to replace in current det
            mm   = excitOp(2,p)    ! MO to replace MO at l
            lpos = posInDeta(ll,n)   ! position (idx) of orb ll in current alpha det (0: not in det)
            mpos = posInDeta(mm,n)   ! position (idx) of orb mm in current alpha det (0: not in det)

            if (lpos /= 0 .and. mpos == 0) then
               select case (mUpdateMode)
               case (UPDATE)
                  call internal_ExcitOpAlphaUpdate(mm,lpos,+1)
               case (SIMPLE)
                  call internal_ExcitOpAlphaSimple(mm,lpos,+1)
               case default
                  call abortp("mocalcderivs: illegal mMode")
               end select
            endif

            ! then do (negative) m -> l excitation: -E_ml = -a^T_l a_m

            if (lpos == 0 .and. mpos /= 0) then
               select case (mUpdateMode)
               case (UPDATE)
                  call internal_ExcitOpAlphaUpdate(ll,mpos,-1)
               case (SIMPLE)
                  call internal_ExcitOpAlphaSimple(ll,mpos,-1)
               case default
                  call abortp("mocalcderivs: illegal mMode")
               end select
            endif

         enddo
      enddo
   enddo

   ! now the same with the beta determinant

   if (na /= nb) then
      deallocate(tmpdet,tmpdeti,tmpdetx,tmpdety,tmpdetz,tmpdet2)
      deallocate(excitedMO,mclistexcit)
      deallocate(tmparr1,tmparr2,tmparr3,tmparr4)
      allocate(tmpdet(nb,nb),tmpdeti(nb,nb),tmpdetx(nb,nb),tmpdety(nb,nb),tmpdetz(nb,nb),tmpdet2(nb,nb))
      allocate(tmparr1(nb,detssize),tmparr2(nb,detssize),tmparr3(nb,detssize),tmparr4(nb,detssize))
      allocate(excitedMO(nb),mclistexcit(nb))
   endif

         ! loop over all orb rot parameters
   do p=1,moParams   ! number of orb rot parameters: pairs of (l->m)
         n = 0
         ! calculate E_lm Phi
         ! loop over all CSF dets, count det number n
         do k=1,ncsf
            do j=1,ndets(k)
               n = n+1
               coeff = cci(k)*ccsf(j,k)

                  ! first do l -> m excitation: E_lm = a^T_m a_l
                  ll   = excitOp(1,p)    ! MO to replace in current det
                  mm   = excitOp(2,p)    ! MO to replace MO at l
                  lpos = posInDetb(ll,n)   ! position (idx) of orb ll in current alpha det (0: not in det)
                  mpos = posInDetb(mm,n)   ! position (idx) of orb mm in current alpha det (0: not in det)

                  if (lpos /= 0 .and. mpos == 0) then
                     select case (mUpdateMode)
                     case (UPDATE)
                        call internal_ExcitOpBetaUpdate(mm,lpos,+1)
                     case (SIMPLE)
                        call internal_ExcitOpBetaSimple(mm,lpos,+1)
                     case default
                        call abortp("mocalcderivs: illegal mMode")
                     end select
                  endif

                  ! then do (negative) m -> l excitation: -E_ml = -a^T_l a_m

                  if (lpos == 0 .and. mpos /= 0) then
                     select case (mUpdateMode)
                     case (UPDATE)
                        call internal_ExcitOpBetaUpdate(ll,mpos,-1)
                     case (SIMPLE)
                        call internal_ExcitOpBetaSimple(ll,mpos,-1)
                     case default
                        call abortp("mocalcderivs: illegal mMode")
                     end select
                  endif

               enddo
            enddo
     enddo

   ! sum laplacian
   do p=1,moParams
      moklapl(p) = sum(moklapli(:,p))
   enddo

   deallocate(tmpdet,tmpdeti,tmpdetx,tmpdety,tmpdetz,tmpdet2)
   deallocate(tmparr1,tmparr2,tmparr3,tmparr4)
   deallocate(excitedMO,mclistexcit)

   contains

      subroutine internal_ExcitOpAlphaUpdate(newmo,pos,sign)
         integer, intent(in) :: newmo,pos,sign
         integer             :: nn
         real(r8)              :: scff
         ! single Excitation operator E_lm(det), pos: position of l in det, newmo=m
         ! sign: +1/-1, for +/- E_lm
         ! calc inv for MO newmo at position of pos (replacing old mo); get det value d for E_lm(det)
         ! this version updates existing inverse for current det

         ! detia, det1xa, etc are not copied for repeated determinants. Use detsRepLst.
        scff=sign*coeff

         if (detsRepLst(n,1) > 0) then
            nn = detsRepLst(n,1)
         else
            nn = n
         endif
        if( detsRepLst(n,1) == 0 ) then !new dets
               darr(n) = deter(1,n)
               excitedMO(1:na) = mat(newmo,1:na,1)
               tmpdeti(1:na,1:na) = detia(1:na,1:na,nn)
               call invrupd(excitedMO,pos,na,na,tmpdeti,darr(n))

               mok(p) = mok(p) + scff * darr(n) * deter(2,n)

               if (calcElecDerivs) then

                     tmpdetx = det1xa(:,1:na,nn)
                     tmpdetx(pos,1:na) = mat1x(newmo,1:na,1)
                     tmpdety = det1ya(:,1:na,nn)
                     tmpdety(pos,1:na) = mat1y(newmo,1:na,1)
                     tmpdetz = det1za(:,1:na,nn)
                     tmpdetz(pos,1:na) = mat1z(newmo,1:na,1)
                     tmpdet2 = det2a(:,1:na,nn)
                     tmpdet2(pos,1:na) = mat2(newmo,1:na,1)
                     tmp0=scff*darr(n) * deter(2,n)
                     do i=1,na
                        ! calculate grad(det) and lapl(det) with determinant-matrix lemma
                        ! update with col for electron i with grad / lapl.
                        ! BUT: use correct MOs including the excited MO
                        tmparr1(i,n) = dot_product(tmpdetx(:,i),tmpdeti(i,:))
                        tmparr2(i,n) = dot_product(tmpdety(:,i),tmpdeti(i,:))
                        tmparr3(i,n) = dot_product(tmpdetz(:,i),tmpdeti(i,:))
                        tmparr4(i,n) = dot_product(tmpdet2(:,i),tmpdeti(i,:))

                        mokgrad(3*i-2,p) = mokgrad(3*i-2,p) + tmparr1(i,n)*tmp0
                        mokgrad(3*i-1,p) = mokgrad(3*i-1,p) + tmparr2(i,n)*tmp0
                        mokgrad(3*i,p)   = mokgrad(3*i,p)   + tmparr3(i,n)*tmp0
                        moklapli(i,p)    = moklapli(i,p)    + tmparr4(i,n)*tmp0
                        !!!write(iul,*) 'eoau:',p,i,tmp1*d,deter(2,n),mokgrad(3*i-2,p)
                        !!!write(iul,*) 'eoau:',p,i,tmp2*d,deter(2,n),mokgrad(3*i-1,p)
                        !!!write(iul,*) 'eoau:',p,i,tmp3*d,deter(2,n),mokgrad(3*i,p)
                        !!!write(iul,*) 'eoau:',p,i,tmp4*d,deter(2,n),moklapli(i,p)
                     enddo
                     tmp5=scff*darr(n)
                     do i=1,nb
                        mokgrad(3*na+3*i-2,p) = mokgrad(3*na+3*i-2,p) + dgrad(3*i-2,2,n) * tmp5
                        mokgrad(3*na+3*i-1,p) = mokgrad(3*na+3*i-1,p) + dgrad(3*i-1,2,n) * tmp5
                        mokgrad(3*na+3*i,p)   = mokgrad(3*na+3*i,p)   + dgrad(3*i,2,n)   * tmp5
                        moklapli(na+i,p)      = moklapli(na+i,p)      + dlapli(i,2,n)    * tmp5
                        !!!write(iul,*) 'eoau:',p,na+i,dgrad(3*i-2,2,n),d,mokgrad(3*na+3*i-2,p)
                        !!!write(iul,*) 'eoau:',p,na+i,dgrad(3*i-1,2,n),d,mokgrad(3*na+3*i-1,p)
                        !!!write(iul,*) 'eoau:',p,na+i,dgrad(3*i,2,n),d,mokgrad(3*na+3*i,p)
                        !!!write(iul,*) 'eoau:',p,na+i,dlapli(i,2,n),d,moklapli(na+i,p)
                     enddo

                     end if
         else  !already calculated dets

                  mok(p) = mok(p) + scff * darr(detsRepLst(n,1)) * deter(2,n)

                  if (calcElecDerivs) then

                  tmp0=scff*darr(detsRepLst(n,1)) * deter(2,n)
                  do i=1,na

                     mokgrad(3*i-2,p) = mokgrad(3*i-2,p) + tmparr1(i,detsRepLst(n,1))*tmp0
                     mokgrad(3*i-1,p) = mokgrad(3*i-1,p) + tmparr2(i,detsRepLst(n,1))*tmp0
                     mokgrad(3*i,p)   = mokgrad(3*i,p)   + tmparr3(i,detsRepLst(n,1))*tmp0
                     moklapli(i,p)    = moklapli(i,p)    + tmparr4(i,detsRepLst(n,1))*tmp0
                  enddo
                  tmp5=scff*darr(detsRepLst(n,1))
                  do i=1,nb
                     mokgrad(3*na+3*i-2,p) = mokgrad(3*na+3*i-2,p) + dgrad(3*i-2,2,n) * tmp5
                     mokgrad(3*na+3*i-1,p) = mokgrad(3*na+3*i-1,p) + dgrad(3*i-1,2,n) * tmp5
                     mokgrad(3*na+3*i,p)   = mokgrad(3*na+3*i,p)   + dgrad(3*i,2,n)   * tmp5
                     moklapli(na+i,p)      = moklapli(na+i,p)      + dlapli(i,2,n)    * tmp5
                  enddo
                  end if
        endif

      end subroutine internal_ExcitOpAlphaUpdate

      subroutine internal_ExcitOpBetaUpdate(newmo,pos,sign)
         integer, intent(in) :: newmo,pos,sign
         integer             :: nn
         real(r8)              :: scff

         ! single Excitation operator E_lm(det), pos: position of l in det, newmo=m
         ! sign: +1/-1, for +/- E_lm
         ! calc inv for MO newmo at position of pos (replacing old mo); get det value d for E_lm(det)
         ! this version updates existing inverse for current det

         ! detib, det1xb, etc are not copied for repeated determinants. Use detsRepLst.
         scff = sign * coeff

         if (detsRepLst(n,2) > 0) then
            nn = detsRepLst(n,2)
         else
            nn = n
         endif
         if(detsRepLst(n,2) == 0) then
               darr(n) = deter(2,n)
               excitedMO(1:nb) = mat(newmo,na+1:ne,1)
               tmpdeti(1:nb,1:nb) = detib(1:nb,1:nb,nn)
               call invrupd(excitedMO,pos,nb,nb,tmpdeti,darr(n))

               mok(p) = mok(p) + scff * darr(n) * deter(1,n)

               if (calcElecDerivs) then

               tmpdetx = det1xb(:,1:nb,nn)
               tmpdetx(pos,1:nb) = mat1x(newmo,na+1:ne,1)
               tmpdety = det1yb(:,1:nb,nn)
               tmpdety(pos,1:nb) = mat1y(newmo,na+1:ne,1)
               tmpdetz = det1zb(:,1:nb,nn)
               tmpdetz(pos,1:nb) = mat1z(newmo,na+1:ne,1)
               tmpdet2 = det2b(:,1:nb,nn)
               tmpdet2(pos,1:nb) = mat2(newmo,na+1:ne,1)
               tmp0=scff*darr(n) * deter(1,n)
               tmp5=scff*darr(n)
               do i=1,na
                  ! calculate grad(det) and lapl(det) with determinant-matrix lemma
                  ! update with col for electron i with grad / lapl.
                  ! BUT: use correct MOs including the excited MO
                  mokgrad(3*i-2,p) = mokgrad(3*i-2,p) + dgrad(3*i-2,1,n) * tmp5
                  mokgrad(3*i-1,p) = mokgrad(3*i-1,p) + dgrad(3*i-1,1,n) * tmp5
                  mokgrad(3*i,p)   = mokgrad(3*i,p)   + dgrad(3*i,1,n)   * tmp5
                  moklapli(i,p)    = moklapli(i,p)    + dlapli(i,1,n)    * tmp5
                  !!!write(iul,*) 'eobu:',p,i,dgrad(3*i-2,1,n),d,mokgrad(3*i-2,p)
                  !!!write(iul,*) 'eobu:',p,i,dgrad(3*i-1,1,n),d,mokgrad(3*i-1,p)
                  !!!write(iul,*) 'eobu:',p,i,dgrad(3*i,1,n),d,mokgrad(3*i,p)
                  !!!write(iul,*) 'eobu:',p,i,dlapli(i,1,n),d,moklapli(i,p)
               enddo
               do i=1,nb
                  tmparr1(i,n) = dot_product(tmpdetx(:,i),tmpdeti(i,:))
                  tmparr2(i,n) = dot_product(tmpdety(:,i),tmpdeti(i,:))
                  tmparr3(i,n) = dot_product(tmpdetz(:,i),tmpdeti(i,:))
                  tmparr4(i,n) = dot_product(tmpdet2(:,i),tmpdeti(i,:))
                  !!!write(iul,*) 'eobu:',p,d
                  !!!write(iul,'(10g12.4)') tmpdeti(i,:)
                  !!!write(iul,'(10g12.4)') tmpdetx(:,i)
                  !!!write(iul,'(10g12.4)') tmpdety(:,i)
                  !!!write(iul,'(10g12.4)') tmpdetz(:,i)
                  !!!write(iul,'(10g12.4)') tmpdet2(:,i)
                  mokgrad(3*na+3*i-2,p) = mokgrad(3*na+3*i-2,p) + tmparr1(i,n) * tmp0
                  mokgrad(3*na+3*i-1,p) = mokgrad(3*na+3*i-1,p) + tmparr2(i,n) * tmp0
                  mokgrad(3*na+3*i,p) = mokgrad(3*na+3*i,p) + tmparr3(i,n) * tmp0
                  moklapli(na+i,p) = moklapli(na+i,p) + tmparr4(i,n) * tmp0
                  !!!write(iul,*) 'eobu:',p,na+i,tmp1*d,deter(1,n),mokgrad(3*na+3*i-2,p)
                  !!!write(iul,*) 'eobu:',p,na+i,tmp2*d,deter(1,n),mokgrad(3*na+3*i-1,p)
                  !!!write(iul,*) 'eobu:',p,na+i,tmp3*d,deter(1,n),mokgrad(3*na+3*i,p)
                  !!!write(iul,*) 'eobu:',p,na+i,tmp4*d,deter(1,n),moklapli(na+i,p)
               enddo

               end if
         else

               mok(p) = mok(p) + scff * darr(detsRepLst(n,2)) * deter(1,n)

               if (calcElecDerivs) then
               tmp0=scff*darr(detsRepLst(n,2)) * deter(1,n)
               tmp5=scff*darr(detsRepLst(n,2))
               do i=1,na
                  mokgrad(3*i-2,p) = mokgrad(3*i-2,p) + dgrad(3*i-2,1,n) * tmp5
                  mokgrad(3*i-1,p) = mokgrad(3*i-1,p) + dgrad(3*i-1,1,n) * tmp5
                  mokgrad(3*i,p)   = mokgrad(3*i,p)   + dgrad(3*i,1,n)   * tmp5
                  moklapli(i,p)    = moklapli(i,p)    + dlapli(i,1,n)    * tmp5
               enddo
               do i=1,nb
                  mokgrad(3*na+3*i-2,p) = mokgrad(3*na+3*i-2,p) + tmparr1(i,detsRepLst(n,2)) * tmp0
                  mokgrad(3*na+3*i-1,p) = mokgrad(3*na+3*i-1,p) + tmparr2(i,detsRepLst(n,2)) * tmp0
                  mokgrad(3*na+3*i,p) = mokgrad(3*na+3*i,p) + tmparr3(i,detsRepLst(n,2)) * tmp0
                  moklapli(na+i,p) = moklapli(na+i,p) + tmparr4(i,detsRepLst(n,2)) * tmp0
               enddo
             endif
          endif

      end subroutine internal_ExcitOpBetaUpdate

      subroutine internal_ExcitOpAlphaSimple(newmo,pos,sign)
         integer, intent(in) :: newmo,pos,sign

         ! single Excitation operator E_lm(det), pos: position of l in det, newmo=m
         ! sign: +1/-1, for +/- E_lm
         ! calc inv for MO newmo at position of pos (replacing old mo); get det value d for E_lm(det)
         ! this version constructs Slater matrix for single excitation of current det from scratch

         ! note: deta contains LU decomp of Slater matrix after call to (lapack_)inv. Reconstruct matrix
         mclistexcit(1:na) = mclist(1:na,n)
         mclistexcit(pos) = newmo
         do i=1,na
            ii = i
            do jj=1,na
               orb = mclistexcit(jj)
               tmpdet(jj,i)   = mat(orb,ii,1)
               tmpdetx(jj,i) = mat1x(orb,ii,1)
               tmpdety(jj,i) = mat1y(orb,ii,1)
               tmpdetz(jj,i) = mat1z(orb,ii,1)
               tmpdet2(jj,i)  = mat2(orb,ii,1)
            enddo
         enddo

         call lapack_inv(na,tmpdet(:,:),tmpdeti(:,:),d,ierr)
         if (ierr /= 0) then
            write(iull,*) "MO rot opt, singular matrix error in process ",mytid
            call abortp("MO rot opt")
         endif

         mok(p) = mok(p) + sign * coeff * d * deter(2,n)

         if (calcElecDerivs) then

         ! construct gradients and laplacians of det as scalar prod with deti
         ! this is O(n**2)
         do i=1,na          ! inverse deti is changed on ALL positions
            tmp1 = dot_product(tmpdetx(:,i),tmpdeti(i,:))
            tmp2 = dot_product(tmpdety(:,i),tmpdeti(i,:))
            tmp3 = dot_product(tmpdetz(:,i),tmpdeti(i,:))
            tmp4 = dot_product(tmpdet2(:,i),tmpdeti(i,:))
            mokgrad(3*i-2,p) = mokgrad(3*i-2,p) + sign * coeff * tmp1*d * deter(2,n)
            mokgrad(3*i-1,p) = mokgrad(3*i-1,p) + sign * coeff * tmp2*d * deter(2,n)
            mokgrad(3*i,p) = mokgrad(3*i,p) + sign * coeff * tmp3*d * deter(2,n)
            moklapli(i,p) = moklapli(i,p) + sign * coeff * tmp4*d * deter(2,n)
            !!!write(iul,*) 'eoas:',p,i,tmp1*d,deter(2,n),mokgrad(3*i-2,p)
            !!!write(iul,*) 'eoas:',p,i,tmp2*d,deter(2,n),mokgrad(3*i-1,p)
            !!!write(iul,*) 'eoas:',p,i,tmp3*d,deter(2,n),mokgrad(3*i,p)
            !!!write(iul,*) 'eoas:',p,i,tmp4*d,deter(2,n),moklapli(i,p)
         enddo
         do i=1,nb
            mokgrad(3*na+3*i-2,p) = mokgrad(3*na+3*i-2,p) + sign * coeff * dgrad(3*i-2,2,n) * d
            mokgrad(3*na+3*i-1,p) = mokgrad(3*na+3*i-1,p) + sign * coeff * dgrad(3*i-1,2,n) * d
            mokgrad(3*na+3*i,p) = mokgrad(3*na+3*i,p) + sign * coeff * dgrad(3*i,2,n) * d
            moklapli(na+i,p) = moklapli(na+i,p) + sign * coeff * dlapli(i,2,n) * d
            !!!write(iul,*) 'eoas:',p,na+i,dgrad(3*i-2,2,n),d,mokgrad(3*na+3*i-2,p)
            !!!write(iul,*) 'eoas:',p,na+i,dgrad(3*i-1,2,n),d,mokgrad(3*na+3*i-2,p)
            !!!write(iul,*) 'eoas:',p,na+i,dgrad(3*i,2,n),d,mokgrad(3*na+3*i-2,p)
            !!!write(iul,*) 'eoas:',p,na+i,dlapli(i,2,n),d,moklapli(na+i,p)
         enddo

         end if
      end subroutine internal_ExcitOpAlphaSimple

      subroutine internal_ExcitOpBetaSimple(newmo,pos,sign)
         integer, intent(in) :: newmo,pos,sign

         ! single Excitation operator E_lm(det), pos: position of l in det, newmo=m
         ! sign: +1/-1, for +/- E_lm
         ! calc inv for MO newmo at position of pos (replacing old mo); get det value d for E_lm(det)
         ! this version constructs Slater matrix for single excitation of current det from scratch

         ! note: deta contains LU decomp of Slater matrix after call to (lapack_)inv. Reconstruct matrix
         call assert(pos <= nb, "mocalcderivs: pos parameter refers to beta determinant: 1 <= pos <= nbeta")
         mclistexcit(1:nb) = mclist(na+1:ne,n)
         mclistexcit(pos) = newmo
         do i=1,nb
            ii = na + i
            do jj=1,nb
               orb = mclistexcit(jj)
               tmpdet(jj,i)  = mat(orb,ii,1)
               tmpdetx(jj,i) = mat1x(orb,ii,1)
               tmpdety(jj,i) = mat1y(orb,ii,1)
               tmpdetz(jj,i) = mat1z(orb,ii,1)
               tmpdet2(jj,i) = mat2(orb,ii,1)
            enddo
         enddo

         call lapack_inv(nb,tmpdet(:,:),tmpdeti(:,:),d,ierr)
         if (ierr /= 0) then
            write(iull,*) "MO rot opt, singular matrix error in process ",mytid
            call abortp("MO rot opt")
         endif

         mok(p) = mok(p) + sign * coeff * deter(1,n) * d

         if (calcElecDerivs) then

         ! construct gradients and laplacians of det as scalar prod with deti
         ! this is O(n**2)
         do i=1,na
            mokgrad(3*i-2,p) = mokgrad(3*i-2,p) + sign * coeff * dgrad(3*i-2,1,n) * d
            mokgrad(3*i-1,p) = mokgrad(3*i-1,p) + sign * coeff * dgrad(3*i-1,1,n) * d
            mokgrad(3*i,p) = mokgrad(3*i,p) + sign * coeff * dgrad(3*i,1,n) * d
            moklapli(i,p) = moklapli(i,p) + sign * coeff * dlapli(i,1,n) * d
            !!!write(iul,*) 'eobs:',p,i,dgrad(3*i-2,1,n),d,mokgrad(3*i-2,p)
            !!!write(iul,*) 'eobs:',p,i,dgrad(3*i-1,1,n),d,mokgrad(3*i-1,p)
            !!!write(iul,*) 'eobs:',p,i,dgrad(3*i,1,n),d,mokgrad(3*i,p)
            !!!write(iul,*) 'eobs:',p,i,dlapli(i,1,n),d,moklapli(i,p)
         enddo
         do i=1,nb          ! inverse deti is changed on ALL positions
            tmp1 = dot_product(tmpdetx(:,i),tmpdeti(i,:))
            tmp2 = dot_product(tmpdety(:,i),tmpdeti(i,:))
            tmp3 = dot_product(tmpdetz(:,i),tmpdeti(i,:))
            tmp4 = dot_product(tmpdet2(:,i),tmpdeti(i,:))
            !!!write(iul,*) 'eobs:',p,d
            !!!write(iul,'(10g12.4)') tmpdeti(i,:)
            !!!write(iul,'(10g12.4)') tmpdetx(:,i)
            !!!write(iul,'(10g12.4)') tmpdety(:,i)
            !!!write(iul,'(10g12.4)') tmpdetz(:,i)
            !!!write(iul,'(10g12.4)') tmpdet2(:,i)
            mokgrad(3*na+3*i-2,p) = mokgrad(3*na+3*i-2,p) + sign * coeff * tmp1*d * deter(1,n)
            mokgrad(3*na+3*i-1,p) = mokgrad(3*na+3*i-1,p) + sign * coeff * tmp2*d * deter(1,n)
            mokgrad(3*na+3*i,p) = mokgrad(3*na+3*i,p) + sign * coeff * tmp3*d * deter(1,n)
            moklapli(na+i,p) = moklapli(na+i,p) + sign * coeff * tmp4*d * deter(1,n)
            !!!write(iul,*) 'eobs:',p,na+i,tmp1*d,deter(1,n),mokgrad(3*na+3*i-2,p)
            !!!write(iul,*) 'eobs:',p,na+i,tmp2*d,deter(1,n),mokgrad(3*na+3*i-1,p)
            !!!write(iul,*) 'eobs:',p,na+i,tmp3*d,deter(1,n),mokgrad(3*na+3*i,p)
            !!!write(iul,*) 'eobs:',p,na+i,tmp4*d,deter(1,n),moklapli(na+i,p)
         enddo

         end if
      end subroutine internal_ExcitOpBetaSimple

   end subroutine moparam_calcderivs


  subroutine moparam_calcderivsOnlyMok(ppMok,ie)

   integer i,ii,j,jj,k,ll,mm,n,na,nb,p,lpos,mpos,orb,ierr
   real(r8)               :: coeff,darr(size(detsRepLst,1)),d
   real(r8), allocatable  :: excitedMO(:), tmpdet(:,:), tmpdeti(:,:)
   integer, allocatable :: mclistexcit(:)
   real(r8), intent(out)  :: ppMok (:)
   integer, intent(in)  :: ie !needed for to know we did one electron update for alpha or beta

   ! assert(nec==1)

   ppMok = 0.d0

   na = nalpha
   nb = nbeta

   ! do alpha electron part

   allocate(tmpdet(na,na),tmpdeti(na,na))
   allocate(excitedMO(na),mclistexcit(na))

   ! calculate E_lm Phi
   ! loop over all CSF dets, count det number n
   do p=1,moParams
          n=0
         do k=1,ncsf
            do j=1,ndets(k)
               n = n+1
               coeff = cci(k)*ccsf(j,k)

               ! loop over all orb rot parameters
                 ! number of orb rot parameters: pairs of (l->m)

                  ! first do l -> m excitation: E_lm = a^T_m a_l
                  ll   = excitOp(1,p)    ! MO to replace in current det
                  mm   = excitOp(2,p)    ! MO to replace MO at l
                  lpos = posInDeta(ll,n)   ! position (idx) of orb ll in current alpha det (0: not in det)
                  mpos = posInDeta(mm,n)   ! position (idx) of orb mm in current alpha det (0: not in det)

                  if (lpos /= 0 .and. mpos == 0) then
                     select case (mUpdateMode)
                     case (UPDATE)
                        call internal_ExcitOpAlphaUpdateOnlyMok(mm,lpos,+1)
                     case (SIMPLE)
                        call internal_ExcitOpAlphaSimpleOnlyMok(mm,lpos,+1)
                     case default
                        call abortp("mocalcderivs: illegal mMode")
                     end select
                  endif

                  ! then do (negative) m -> l excitation: -E_ml = -a^T_l a_m

                  if (lpos == 0 .and. mpos /= 0) then
                     select case (mUpdateMode)
                     case (UPDATE)
                        call internal_ExcitOpAlphaUpdateOnlyMok(ll,mpos,-1)
                     case (SIMPLE)
                        call internal_ExcitOpAlphaSimpleOnlyMok(ll,mpos,-1)
                     case default
                        call abortp("mocalcderivs: illegal mMode")
                     end select
                  endif

               enddo !dets
            enddo    !CSFs

   enddo       !MOparams

   ! now the same with the beta determinant

   if (na /= nb) then
      deallocate(tmpdet,tmpdeti)
      deallocate(excitedMO,mclistexcit)
      allocate(tmpdet(nb,nb),tmpdeti(nb,nb))
      allocate(excitedMO(nb),mclistexcit(nb))
   endif

   n = 0
   ! calculate E_lm Phi
   ! loop over all CSF dets, count det number n
do p=1,moParams
   n=0
   do k=1,ncsf
      do j=1,ndets(k)
         n = n+1
         coeff = cci(k)*ccsf(j,k)

         ! loop over all orb rot parameters
         ! number of orb rot parameters: pairs of (l->m)

            ! first do l -> m excitation: E_lm = a^T_m a_l
            ll   = excitOp(1,p)    ! MO to replace in current det
            mm   = excitOp(2,p)    ! MO to replace MO at l
            lpos = posInDetb(ll,n)   ! position (idx) of orb ll in current alpha det (0: not in det)
            mpos = posInDetb(mm,n)   ! position (idx) of orb mm in current alpha det (0: not in det)

            if (lpos /= 0 .and. mpos == 0) then
               select case (mUpdateMode)
               case (UPDATE)
                  call internal_ExcitOpBetaUpdateOnlyMok(mm,lpos,+1)
               case (SIMPLE)
                  call internal_ExcitOpBetaSimpleOnlyMok(mm,lpos,+1)
               case default
                  call abortp("mocalcderivs: illegal mMode")
               end select
            endif

            ! then do (negative) m -> l excitation: -E_ml = -a^T_l a_m

            if (lpos == 0 .and. mpos /= 0) then

               select case (mUpdateMode)
               case (UPDATE)
                  call internal_ExcitOpBetaUpdateOnlyMok(ll,mpos,-1)
               case (SIMPLE)
                  call internal_ExcitOpBetaSimpleOnlyMok(ll,mpos,-1)
               case default
                  call abortp("mocalcderivs: illegal mMode")
               end select
            endif

         enddo
      enddo
   enddo

   ! sum laplacian
   deallocate(tmpdet,tmpdeti)
   deallocate(excitedMO,mclistexcit)

   contains

      subroutine internal_ExcitOpAlphaUpdateOnlyMok(newmo,pos,sign)
         integer, intent(in) :: newmo,pos,sign
         integer nn

         ! single Excitation operator E_lm(det), pos: position of l in det, newmo=m
         ! sign: +1/-1, for +/- E_lm
         ! calc inv for MO newmo at position of pos (replacing old mo); get det value d for E_lm(det)
         ! this version updates existing inverse for current det

         ! detia, det1xa, etc are not copied for repeated determinants. Use detsRepLst.
         if (detsRepLst(n,1) > 0) then
            nn = detsRepLst(n,1)
         else
            nn = n
         endif
         if (detsRepLst(n,1)==0) then
            darr(n) = deter(1,n)
            excitedMO(1:na) = mat(newmo,1:na,1)
            if(ie .le. nalpha) then
               call invdetcalc(excitedMO,na,na,detiaOne(1:na,pos,nn),darr(n))
                else
               call invdetcalc(excitedMO,na,na,detia(1:na,pos,nn),darr(n))
            endif
          else
             darr(n)=darr(detsRepLst(n,1))
         endif

         ppMok(p) = ppMok(p) + sign * coeff * darr(n) * deter(2,n)


      end subroutine internal_ExcitOpAlphaUpdateOnlyMok

      subroutine internal_ExcitOpBetaUpdateOnlyMok(newmo,pos,sign)
         integer, intent(in) :: newmo,pos,sign
         integer nn

         ! single Excitation operator E_lm(det), pos: position of l in det, newmo=m
         ! sign: +1/-1, for +/- E_lm
         ! calc inv for MO newmo at position of pos (replacing old mo); get det value d for E_lm(det)
         ! this version updates existing inverse for current det

         ! detib, det1xb, etc are not copied for repeated determinants. Use detsRepLst.
         if (detsRepLst(n,2) > 0) then
            nn = detsRepLst(n,2)
         else
            nn = n
         endif
         if (detsRepLst(n,2)==0) then
            darr(n) = deter(2,n)
            excitedMO(1:nb) = mat(newmo,na+1:ne,1)
            if (ie > nalpha) then
                call invdetcalc(excitedMO,nb,nb,detibOne(1:nb,pos,nn),darr(n))
             else
                call invdetcalc(excitedMO,nb,nb,detib(1:nb,pos,nn),darr(n))
            endif
             else
              darr(n)=darr(detsRepLst(n,2))
         endif

         ppMok(p) = ppMok(p) + sign * coeff * darr(n) * deter(1,n)

      end subroutine internal_ExcitOpBetaUpdateOnlyMok

      subroutine internal_ExcitOpAlphaSimpleOnlyMok(newmo,pos,sign)
         integer, intent(in) :: newmo,pos,sign

         ! single Excitation operator E_lm(det), pos: position of l in det, newmo=m
         ! sign: +1/-1, for +/- E_lm
         ! calc inv for MO newmo at position of pos (replacing old mo); get det value d for E_lm(det)
         ! this version constructs Slater matrix for single excitation of current det from scratch

         ! note: deta contains LU decomp of Slater matrix after call to (lapack_)inv. Reconstruct matrix
         mclistexcit(1:na) = mclist(1:na,n)
         mclistexcit(pos) = newmo
         do i=1,na
            ii = i
            do jj=1,na
               orb = mclistexcit(jj)
               tmpdet(jj,i)   = mat(orb,ii,1)
            enddo
         enddo

         call lapack_inv(na,tmpdet(:,:),tmpdeti(:,:),d,ierr)
         if (ierr /= 0) then
            write(iull,*) "MO rot opt, singular matrix error in process ",mytid
            call abortp("MO rot opt")
         endif

         ppMok(p) = ppMok(p) + sign * coeff * d * deter(2,n)

      end subroutine internal_ExcitOpAlphaSimpleOnlyMok

      subroutine internal_ExcitOpBetaSimpleOnlyMok(newmo,pos,sign)
         integer, intent(in) :: newmo,pos,sign

         ! single Excitation operator E_lm(det), pos: position of l in det, newmo=m
         ! sign: +1/-1, for +/- E_lm
         ! calc inv for MO newmo at position of pos (replacing old mo); get det value d for E_lm(det)
         ! this version constructs Slater matrix for single excitation of current det from scratch

         ! note: deta contains LU decomp of Slater matrix after call to (lapack_)inv. Reconstruct matrix
         call assert(pos <= nb, "mocalcderivs: pos parameter refers to beta determinant: 1 <= pos <= nbeta")
         mclistexcit(1:nb) = mclist(na+1:ne,n)
         mclistexcit(pos) = newmo
         do i=1,nb
            ii = na + i
            do jj=1,nb
               orb = mclistexcit(jj)
               tmpdet(jj,i)  = mat(orb,ii,1)
            enddo
         enddo

         call lapack_inv(nb,tmpdet(:,:),tmpdeti(:,:),d,ierr)
         if (ierr /= 0) then
            write(iull,*) "MO rot opt, singular matrix error in process ",mytid
            call abortp("MO rot opt")
         endif

         ppMok(p) = ppMok(p) + sign * coeff * deter(1,n) * d

      end subroutine internal_ExcitOpBetaSimpleOnlyMok

   end subroutine moparam_calcderivsOnlyMok

   integer function getMOParamsCnt()
      getMOParamsCnt = moParams
   end function getMOParamsCnt


   subroutine putMOParamsVector(optMode,p)
   !-------------------------------------!
      integer, intent(in)   :: optMode       ! optimization mode ! ignored
      real(r8), intent(in)    :: p(:)         ! parameter vector
      integer k, l, m, ierr, i, j, mu, lwork
      real(r8) X(norb,norb), X2(norb,norb), W(norb,norb), R(norb,norb), W1(norb,norb)
      real(r8) lambda(norb), tau(norb), work(4*norb)
      real(r8) zero, one
      character*1 uplo1,uplo2
      real(r8), parameter :: EPS = 1.d-8

      ! overwrite current MO matrix
      cmo = cmo0

      if (mParamMode==SUCCESSIVE) then
         ! successive 2x2 rotations for each l,m pair
         do k=1,moParams   ! number of orb rot parameters: pairs of (l->m)
            l   = excitOp(1,k)    ! MO to replace in current det
            m   = excitOp(2,k)    ! MO to replace MO at l
            if (l>norb .or. m>norb) call abortp("orbital optimization: inconsistent parameters (1) in putMOParamsVector")
            call rotateMO(l,m,p(k))
         enddo
      else if (mParamMode==EXPKAPPA) then
         ! Eigenvector method for squared kappa matrix, see Helgaker book eq. 3.1.31
         X = 0
         do k=1,moParams   ! number of orb rot parameters: pairs of (l->m)
            l   = excitOp(1,k)    ! MO to replace in current det
            m   = excitOp(2,k)    ! MO to replace MO at l
            if (l>m) call abortp("orbital optimization: illegal parameter in putMOParamsVector")
            if (l>norb .or. m>norb) call abortp("orbital optimization: inconsistent parameters (2) in putMOParamsVector")
            X(l,m) = p(k)
            X(m,l) = -p(k)
         enddo

         ! X2 = X*X
         uplo1 = 'N'
         uplo2 = 'N'
         one = 1.D0
         zero = 0.D0
         call dgemm(uplo1, uplo2, norb, norb, norb, one, X, norb, X, norb, zero, X2, norb)
         ! get eigenvalues lambda and eigenvectors W for X2
         lwork = 4*norb
         W = X2   
         call dsyev('V', 'U', norb, W, norb, lambda, work, lwork, ierr)
         if (ierr /= 0) call abortp('(putMOParamsVector): diagonalization failed')

         do i=1,norb
            if (lambda(i) > EPS ) then
               write(iull,*) "WARNING: positive eigenvalue:",i,lambda(i)
            endif
            if (lambda(i) > zero .and. lambda(i) < EPS ) then
                lambda(i)= zero
                if (MASTER)  then
                   write(iull,*) " WARNING: very small positive eigenvalue is set to zero",i
                end if
             end if
            tau(i) = sqrt(-lambda(i))
         enddo


         if (MASTER .and. logmode>=4) then
            write(iul,'(/a)') " eigenvalues: lambda, tau:"
            do i=1,norb
               write(iul,'(2g20.8)') lambda(i), tau(i)
            enddo
            write(iul,'(/a)') " eigenvectors: W:"
            do i=1,norb
               write(iul,'(12g10.4)') (W(i,j),j=1,norb)
            enddo
         endif

         ! W * cos(tau), cos(tau) as diagonal matrix
         do i=1,norb
            W1(:,i) = W(:,i) * cos(tau(i))
         enddo
         ! R = W * cos(tau) * W^T
         uplo1 = 'N'
         uplo2 = 'T'
         call dgemm(uplo1,uplo2,norb,norb,norb,one,W1,norb,W,norb,zero,R,norb)

         ! W * tau^-1 * sin(tau), again though of as diagonal matrices
         ! note: sin(tau)/tau = 1 for tau -> 0
         do i=1,norb
            if (abs(tau(i)) > EPS) then
               W1(:,i) = W(:,i) * sin(tau(i)) / tau(i)
            else
               W1(:,i) = W(:,i)
            endif
         enddo

         ! overwrite X2 = W1*W^T, R = R + X2*X
         uplo1 = 'N'
         uplo2 = 'T'
         call dgemm(uplo1,uplo2,norb,norb,norb,one,W1,norb,W,norb,zero,X2,norb)
         uplo1 = 'N'
         uplo2 = 'N'
         call dgemm(uplo1,uplo2,norb,norb,norb,one,X2,norb,X,norb,one,R,norb)

         if (MASTER .and. logmode >=4) then
            write(iul,'(/a)') ' rotation matrix R:'
            do i=1,norb
               do j=1,norb
                  write(iul,'(2i5,g20.8)') i,j,R(i,j)
               enddo
            enddo
         endif

         ! now multiply C_MO with the rotation matrix (from rhs!!)
         ! Multiplying with transposed matrix for consistency with calculated derivatives CHECK!
         uplo1 = 'N'
         uplo2 = 'T'
         call dgemm(uplo1,uplo2,nbas,norb,norb,one,cmo0,nbas,R,norb,zero,cmo,nbas)

      else
         call abortp("orbital optimization: illegal optMode in putMOParamsVector")
      endif

      p_save = p

      if (MASTER .and. mVerbose >= 4) then
         write(iul,*) ' putMOParamsVector: display of all changed values of MO coeff matrix '
         write(iul,*) '    mu      j        c_mu,j (orig)  ->   c_mu,j (new)'
         do j=1,norb
            do mu=1,nbas
               if (abs(cmo(mu,j) - cmo0(mu,j)) > 1.d-10) then
                  write(iul,'(2i5,2f20.10)') mu,j,cmo0(mu,j),cmo(mu,j)
               endif
            enddo
         enddo
      endif

   end subroutine putMOParamsVector


   subroutine rotateMO(i,j,kappa)
      integer, intent(in) :: i,j      ! indices of MOs to rotate
      real(r8), intent(in)  :: kappa    ! rotation angle
      integer mu
      real(r8) ci,cj

      do mu=1,nbas
         ci =  cos(kappa)*cmo(mu,i) + sin(kappa)*cmo(mu,j)
         cj = -sin(kappa)*cmo(mu,i) + cos(kappa)*cmo(mu,j)
         cmo(mu,i) = ci
         cmo(mu,j) = cj
      enddo
   end subroutine rotateMO




  subroutine getMOParamsVector(optMode,p)
  !-------------------------------------!
     integer, intent(in)   :: optMode       ! optimization mode  ! ignored !
     real(r8), intent(inout) :: p(:)         ! parameter vector
     integer k,k0

     call assert(size(p) == size(p_save),"orbital optimization: illegal parameter vector length")
     p = p_save
   end subroutine getMOParamsVector

   pure integer function getMoUpdateMode()
  !-------------------------------------!
  ! return mo update mode
      getMoUpdateMode = mUpdateMode
   end function getMoUpdateMode


   function withSymmetriseList() result(res)
      logical res
      res = mParamSymEntries > 0
   end function withSymmetriseList


   function getNParamSymEntries() result(res)
      integer res
      res = mParamSymEntries
   end function getNParamSymEntries


   subroutine getParamSymmetriseList(idx,list)
      integer, intent(in)                 :: idx
      integer, allocatable, intent(inout) :: list(:)
      if (allocated(list)) deallocate(list)
      allocate(list(paramSymmetriseList(0,idx)))
      list = paramSymmetriseList(1:paramSymmetriseList(0,idx),idx)
   end subroutine


   subroutine findEquivalentRots(MOSymmetriseList, p0, e1, e2, entryVec)
   !-------------------------------------------------------------------!
      ! find all parameters that are symmetry equivalent to the e1->e2 rotation
      ! using the excitOp list and MOSymmetriseList
      integer, intent(in)    :: MOSymmetriseList(0:,:)
      integer, intent(in)    :: p0          ! current parameter index
      integer, intent(in)    :: e1, e2      ! corresponding orbital rotation
      integer, intent(inout) :: entryVec(0:)   ! list of equivalent parameters (0: len of list)
      integer idx, mo1, mo2, mosl1, mosl2, i, j, p

      entryVec = 0
      mosl1 = findInList(e1)
      mosl2 = findInList(e2)

      !!!write(iul,*) "DBG:fer1:", e1, e2, mosl1, mosl2

      if (mosl1 > 0 .and. mosl2 > 0) then
         idx = 1
         entryVec(idx) = p0
         do i = 1, MOSymmetriseList(0,mosl1)
            mo1 = MOSymmetriseList(i,mosl1)
            do j = 1, MOSymmetriseList(0,mosl2)
               mo2 = MOSymmetriseList(j,mosl2)
               !!!write(iul,*) "DBG:fer:", idx, i, j, mo1, mo2
               p = getParamNumber(mo1,mo2)
               !!!write(iul,*) "DBG:fer:", p
               if (p > p0) then ! add p to list
                  idx = idx + 1
                  if (idx > size(entryVec)) call abortp("findEquivalentRots: too many entries!")
                  entryVec(idx) = p
                  !!!write(iul,*) "DBG:fer:add:", idx, p
               end if
            end do
         end do
         if (idx > 1) then
            entryVec(0) = idx
         else
            entryVec = 0
         end if
      end if

   contains

      function findInList(e) result(res)
         integer, intent(in) :: e
         integer             :: res
         integer idx, i

         res = 0
         OUTER: do idx = 1, size(MOSymmetriseList,2)
            INNER: do i = 1, MOSymmetriseList(0,idx)
               if (e == MOSymmetriseList(i,idx)) then
                  res = idx
                  exit OUTER
               end if
            end do INNER
         end do OUTER
      end function findInList

      function getParamNumber(mo1, mo2) result(res)
         integer, intent(in) :: mo1, mo2
         integer             :: res
         integer i, np

         np = size(excitOp,2)
         res = 0
         do i = 1, np
            if (mo1 == excitOp(1,i) .and. mo2 == excitOp(2,i)) then
               res = i
               exit
            end if
         end do
      end function getParamNumber

   end subroutine

end module moParam_m
