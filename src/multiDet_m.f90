! Copyright (C) 1996, 2006-2007, 2011-2012, 2015, 2018 Arne Luechow
! Copyright (C) 2000 Sebastian Manten
! Copyright (C) 2005 Tony C. Scott
! Copyright (C) 2013 Alexander Sturm
! Copyright (C) 2014-2017 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later


! F90 module 'mdet' for handling multi determinant part of wavefunction

! QMC program amolqc. rewritten in 2006. AL
!
! free form f90 version. AL, 2011

!---------!
MODULE multiDet_m
!---------!

use kinds_m, only: r8
use error_m
use wfData_m
use mos_m, only: mat,mat1x,mat1y,mat1z,mat2,mMOElecConfigs
use utils_m, only: tokenize
use linAlg_m, only: invcupd, invdetcalc, invrupd, inv1
use, intrinsic :: ieee_arithmetic, only: ieee_is_normal

implicit none

#ifdef MAXDETSPERCSF
   integer, parameter :: ndetcsfmax = MAXDETSPERCSF  ! max # of dets per csf
#else
   integer, parameter :: ndetcsfmax = 120  ! max # of dets per csf
#endif
   integer, parameter :: MDET_NONE=0, MDET_LU_ERR=1, MDET_INV_ERR=2, MDET_LU_ZERO_DET=3

   integer, allocatable :: mclist(:,:) ! configuration list
   integer, allocatable :: ndets(:)    ! # dets in k-th CSF
   integer, allocatable :: detsRepLst(:,:)  ! list of repeated determinats in CSFs
   real(r8), allocatable  :: cci(:)      ! CI (multi CSF!) coeff.
   real(r8), allocatable  :: ccsf(:,:)   ! CSF coefficients
   ! determinants and its derivatives
   ! deta,b(j,i,nci)  : j-th MO, i-th el., nci-th CI-det, for a(lpha) and b(eta) spin
   ! det1x,y,z|a,b           : as det, but containing x,y,z derivatives
   ! det2a,b                : as det, but containing laplacians
   real(r8), allocatable  :: deta(:,:,:),detb(:,:,:), det1xa(:,:,:),det1xb(:,:,:), &
                           det1ya(:,:,:),det1yb(:,:,:), det1za(:,:,:),det1zb(:,:,:), &
                           det2a(:,:,:), det2b(:,:,:)
   ! corresponding inverse matrices
   real(r8), allocatable :: detia(:,:,:),detib(:,:,:), detiaOne(:,:,:),detibOne(:,:,:)
   real(r8), allocatable :: det1xia(:,:),det1xib(:,:), &
                          det1yia(:,:),det1yib(:,:),det1zia(:,:),det1zib(:,:), &
                          det2ia(:,:), det2ib(:,:)
   ! deter(ispin,nci)    : respective determinants of det
   ! dgrad(3*ne,ispin,nci),dlapli(ne,ispin,nci) : gradient and laplacian per electron
   real(r8), allocatable :: deter(:,:),dgrad(:,:,:),dlapli(:,:,:), olddet(:)
   logical :: dets_reordered = .false.

   integer :: nexalpha, nexbeta        ! number of virtual alpha/beta orbitals
   integer, allocatable :: virt(:,:)                ! indices of virtual orbitals
   integer, allocatable :: nexc(:,:)   ! number of excitations for each determinant
   integer, allocatable :: exidx(:,:,:,:)! indices of excited orbitals for each det
   real(r8), allocatable :: Ta(:,:), Tb(:,:)
   real(r8), allocatable :: T1xa(:,:), T1ya(:,:), T1za(:,:), T2a(:, :)
   real(r8), allocatable :: T1xb(:,:), T1yb(:,:), T1zb(:,:), T2b(:, :)
   real(r8), allocatable :: excma(:,:), excwa(:,:), excm1xa(:,:), excm1ya(:,:), &
                          excm1za(:,:), excm2a(:,:), excmb(:,:), excwb(:,:), &
                          excm1xb(:,:), excm1yb(:,:), excm1zb(:,:), excm2b(:,:)

CONTAINS

   !----------------------------!
   subroutine mdetinput(lines,nl)
   !----------------------------!

! initialize module from input
! mdetinput read # of determinants and the corresponding determinants
! from 'lines' and bring them in maximal coincidence
! one determinant CSFs are generated here!

   character(len=*), intent(in) :: lines(:)! lines array
   integer, intent(in)          :: nl      ! actual # of lines
   integer i,j,m,n,k,alstat,ios,mmm,nnn
   integer ii1,ii2,nw,na
   integer spin,tmp
   character(len=40) :: words(5)

   call assert(.not. allocated(cci) .and. .not. allocated(deta),'mdetinput: data already allocated')

   call tokenize(lines(2),words,nw)
   if (nw >= 2) then
      ! new format
      if (words(1)=='single') then  ! single determinant
         ncsf = 1; ndet = 1
         allocate(cci(ndet),mclist(ne,ndet),ndets(ndet),ccsf(1,ndet),stat=alstat)
         call assert(alstat==0, 'mdetinput: allocation failed')
         cci(1) = 1.d0
         ndets(1) = 1
         ccsf(1,1) = 1.d0
         if (words(2)=='restricted') then
            na = getNAlpha()
            do i=1,na
               mclist(i,1) = i
            end do
            do i=1,getNBeta()
               mclist(na+i,1) = i
            end do
         else if (words(2)=='unrestricted') then
            do i=1,getNElec()
               mclist(i,1) = i
            end do
         else
            call abortp('$dets: single det input: "single [un]restricted"')
         end if
      else if (words(1)=='multi') then
         call abortp('$dets: multi det input not yet implemented')
      else
         call abortp('$dets: format: "single|multi [un]restricted"')
      end if
   else  ! old format
      read(lines(2),*,iostat=ios) ncsf
      call assert(ios==0 .and. ncsf>0, 'mdetinput: reading # dets failed')

      ndet = ncsf
      allocate(cci(ndet),mclist(ne,ndet),ndets(ndet),ccsf(1,ndet),stat=alstat)
      call assert(alstat==0, 'mdetinput: allocation failed')


      do n=1,ncsf
         ! Give each determinant pair as orbital occupancy. First alpha
         ! electrons then beta electrons
         read(lines(n+2),*,iostat=ios) cci(n),(mclist(i,n),i=1,ne)
         call assert(ios==0,'mdetinput: reading determinants in det mode failed')
         ndets(n) = 1
         ccsf(1,n) = 1.0d0
      enddo
      ! bring excited determinants in maximal coincidence with reference determinant
      ! (required for proper work of mdetcalc) !
      call reorder_dets()

   end if

   allocate(detsRepLst(ndet,2),stat=alstat)
   call assert(alstat==0,'mdetinput: allocate2 failed')
   detsRepLst=0
   if (repeatedDetsOpt) call build_replst()
   ! allocate determinant arrays
   call mdet_alloc()

   end subroutine mdetinput


   !-------------------------------!
   subroutine mdetcsfinput(lines,nl)
   !-------------------------------!

   ! mdetinput read # of csfs and the corresponding determinants
   ! from 'lines' and bring them in maximal coincidence
   character(len=*), intent(in) :: lines(:)! lines array
   integer, intent(in)          :: nl      ! actual # of lines
   integer i,j,m,n,k,ios,alstat,idx
   integer ii1,ii2,nw,na,nnn,mmm
   integer spin,tmp
   integer, allocatable :: mclist0(:,:)
   character(len=40) :: words(5)

   call assert(.not. fastmdet, "fastdet not yet implemented for CSFs")

   call tokenize(lines(2),words,nw)
   if (nw >= 2) then
      ! new format
      if (words(1)=='single') then  ! single determinant
         ncsf = 1; ndet = 1
         allocate(cci(ncsf),mclist0(ne,ncsf*ndetcsfmax),ndets(ncsf),ccsf(ndetcsfmax,ncsf),stat=alstat)
         call assert(alstat==0, 'mcsfinput: allocation failed')
         cci(1) = 1.d0
         ndets(1) = 1
         ccsf(1,1) = 1.d0
         if (words(2)=='restricted') then
            na = getNAlpha()
            do i=1,na
               mclist0(i,1) = i
            end do
            do i=1,getNBeta()
               mclist0(na+i,1) = i
            end do
         else if (words(2)=='unrestricted') then
            do i=1,getNElec()
               mclist0(i,1) = i
            end do
         else
            call abortp('$dets: single det input: "single [un]restricted"')
         end if
      else if (words(1)=='multi') then
         call abortp('$dets: multi det input not yet implemented')
      else
         call abortp('$dets: format: "single|multi [un]restricted"')
      end if
   else  ! old format
      read(lines(2),*,iostat=ios) ncsf
      call assert(ios==0 .and. ncsf>0, 'mdetcsfinput: reading # csfs failed')

      allocate(cci(ncsf),mclist0(ne,ncsf*ndetcsfmax),ndets(ncsf),ccsf(ndetcsfmax,ncsf),stat=alstat)
      call assert(alstat==0, 'mdetinput: allocation failed')


      n = 0    ! count for dets
      ndet = 0 ! ditto
      idx=3
      do k=1,ncsf
         read(lines(idx),*,iostat=ios) cci(k),ndets(k)
         call assert(ios==0,'mdetcsfinput: reading ci coeff and ndet per csf failed')
         call assert(ndets(k)<=ndetcsfmax,'mdetcsfinput: ndet per csf too large. Compile with -DMAXDETSPERCSF=___')
         ndet = ndet + ndets(k)
         ! Give each determinant pair as orbital occupancy. First alpha
         ! electrons then beta electrons
         do j=1,ndets(k)
            n = n+1
            read(lines(idx+j),*,iostat=ios) ccsf(j,k),(mclist0(i,n),i=1,ne)
            call assert(ios==0,'mdetcsfinput: reading csf coeff and det failed')
         enddo
         idx=idx+ndets(k)+1
      enddo
      call assert(n==ndet,"mdetcsfinput: something strange is wrong!")
   end if

   allocate(mclist(ne,ndet),stat=alstat)
   call assert(alstat==0,'mdetcsfinput: allocate2 failed')
   allocate(detsRepLst(ndet,2),stat=alstat)
   call assert(alstat==0,'mdetcsfinput: allocate3 failed')
   detsRepLst = 0
   mclist(:,:) = mclist0(:,1:ndet)
   deallocate(mclist0)
   call reorder_dets()
   if (repeatedDetsOpt) call build_replst()

   ! allocate determinant arrays
   call mdet_alloc()

   end subroutine mdetcsfinput



   subroutine build_replst()
   !-----------------------!
   !build list of repeated determinants in csfs
   integer :: detscnt,elcnt ! determinants and electrons counter
   integer :: jj,uniqalpha=0,uniqbeta=0
   logical :: repeated
   detsRepLst(1,:)=(0)
   !first for alpha dets
   do detscnt=2,ndet
      repeated=.FALSE.
      jj=1
      do while (detscnt>jj .and. ( repeated .eqv. .FALSE.) )
         if (all(mclist(1:nalpha,detscnt)==mclist(1:nalpha,jj))) then
            repeated=.TRUE.
            detsRepLst(detscnt,1)=(jj)
         endif
         jj=jj+1
      end do
      if (repeated .eqv. .FALSE.) detsRepLst(detscnt,1)=0
   end do
   !first for alpha dets
   do detscnt=2,ndet
      repeated=.FALSE.
      jj=1
      do while (detscnt>jj .and. ( repeated .eqv. .FALSE.) )
         if (all(mclist(nalpha+1:,detscnt)==mclist(nalpha+1:,jj))) then
            repeated=.TRUE.
            detsRepLst(detscnt,2)=(jj)
         endif
         jj=jj+1
      end do
      if (repeated .eqv. .FALSE.) detsRepLst(detscnt,2)=0
   end do

   do jj=1,ndet
      if (detsRepLst(jj,1) == 0) uniqalpha = uniqalpha + 1
      if (detsRepLst(jj,2) == 0) uniqbeta = uniqbeta + 1
   end do
   if (logmode>=2  .and. MASTER) then
           write(iul,*)""
           write(iul,'(a,i8,a,i8)') '    NCSFs = ', ncsf,' Ndets = ',ndet
           write(iul,'(i5,a,i5,a)') uniqalpha, " unique alpha and", uniqbeta, " unique beta determinants"
           write(iul,*) "                         are going to be calculated."
           write(iul,*) ""
   endif
   end subroutine build_replst


   !---------------------!
   subroutine mdet_alloc()
   !---------------------!
   integer :: alstat
   ! matrices
   allocate(deta(nalpha,nalpha,ndet),detb(nbeta,nbeta,ndet),     &
            det1xa(nalpha,nalpha,ndet),det1xb(nbeta,nbeta,ndet), &
            det1ya(nalpha,nalpha,ndet),det1yb(nbeta,nbeta,ndet), &
            det1za(nalpha,nalpha,ndet),det1zb(nbeta,nbeta,ndet), &
            det2a(nalpha,nalpha,ndet),det2b(nbeta,nbeta,ndet),olddet(ndet), &
            stat=alstat)
   call assert(alstat==0,'mdet_alloc: allocation1 failed')
   ! inverse matrices
   allocate(detia(nalpha,nalpha,ndet),detib(nbeta,nbeta,ndet), stat=alstat)
   call assert(alstat==0,'mdet_alloc: allocation2 failed')
   allocate(deter(2,ndet),dgrad(3*ne,2,ndet),dlapli(ne,2,ndet),stat=alstat)
   call assert(alstat==0,'mdet_alloc: allocation3 failed')
   allocate(detiaOne(nalpha,nalpha,ndet),detibOne(nbeta,nbeta,ndet), stat=alstat)
   call assert(alstat==0,'mdet_alloc: allocation2 failed')
   if (fastmdet) then
      allocate(det1xia(nalpha,nalpha),det1xib(nbeta,nbeta), &
               det1yia(nalpha,nalpha),det1yib(nbeta,nbeta), &
               det1zia(nalpha,nalpha),det1zib(nbeta,nbeta), &
               det2ia(nalpha,nalpha),det2ib(nbeta,nbeta), &
               stat=alstat)
      call assert(alstat==0,'mdet_alloc: allocation4 failed')
      allocate(Ta(norb,norb), Tb(norb,norb), &
               T1xa(norb,norb), T1ya(norb,norb), T1za(norb,norb), T2a(norb, norb), &
               T1xb(norb,norb), T1yb(norb,norb), T1zb(norb,norb), T2b(norb, norb), &
               stat=alstat)
      call assert(alstat==0,'mdet_alloc: allocation5 failed')
      allocate(excma(nexalpha,nexalpha), excwa(nexalpha,nexalpha), &
               excm1xa(nexalpha,nexalpha), excm1ya(nexalpha,nexalpha), &
               excm1za(nexalpha,nexalpha), excm2a(nexalpha,nexalpha), &
               excmb(nexbeta,nexbeta), excwb(nexbeta,nexbeta), &
               excm1xb(nexbeta,nexbeta), excm1yb(nexbeta,nexbeta), &
               excm1zb(nexbeta,nexbeta), excm2b(nexbeta,nexbeta), stat=alstat)
      call assert(alstat==0,'mdet_alloc: allocation6 failed')
   endif
   end subroutine


   subroutine reorder_dets()
   !-----------------------!

   ! reorder the mos in all determinants such that maximal coincidence with
   ! the 1st determinant (as reference determinant) is obtained
   ! also generate the list of used virtual alpha and beta orbitals for each
   ! determinant
   integer spin,j,ii1,ii2,m,n,exc,alstat,i,k,l,mmm,nnn,jj
   real(r8) tmp

   do j=1,ndet
      jj=j
      mmm=0
      do while (jj>0)
         mmm=mmm+1
         jj=jj-ndets(mmm)
      enddo
      nnn=ndets(mmm)+jj
      do spin=1,2
         if (spin==1) then
            ii1=1
            ii2=nalpha
         else
            ii1=nalpha+1
            ii2=ne
         endif

         !!  ensure that each MO appears only once in each determinant
         do m=ii1, ii2
           do n=ii1, ii2
             if (mclist(m,j)==mclist(n,j) .and. (n/=m)) then
               call abortp('mdetcalc: inconsistent mdet input')
             endif
           enddo
         enddo

         !! bring in maximal coincidence
         do m=ii1, ii2
            if (mclist(m,j)/=mclist(m,1)) then
              do n=ii1, ii2
                if (mclist(n,j)==mclist(m,1) .and. (m/=n)) then
                  tmp = mclist(n,j)
                  mclist(n,j) = mclist(m,j)
                  mclist(m,j) = tmp
                  ccsf(nnn,mmm) = -ccsf(nnn,mmm)
                  if (.not.(dets_reordered)) dets_reordered = .true.
                endif
              enddo
            endif
         enddo
      enddo !spin loop
   enddo

   if (logmode>=2 .and. (dets_reordered) .and. MASTER) then
      write(iul,*)
      write(iul,*) '------------------------------------------'
      write(iul,*) 'Multi determinant input has been modified '
      write(iul,*) 'in order to get maximal coincidence of the'
      write(iul,*) 'excited determinants and the reference.   '
      write(iul,*) '------------------------------------------'
      write(iul,*)
      if (logmode>=3) then
         write(iul,*) 'new det/csf list:'
         n=0
         do k=1,ncsf
            do l=1,ndets(k)
               n = n + 1
               write(iul,'(i4,i3,i4,2f8.5,100i3)') k,l,n,cci(k),ccsf(l,k),(mclist(i,n),i=1,ne)
            enddo
         enddo
         write(iul,*)
      endif
   endif

   if (fastmdet) then
      allocate(virt(norb,2),stat=alstat)
      call assert(alstat==0, 'reoder_dets: allocation failed')
      nexalpha = 0
      nexbeta = 0
      do n = 2, ndet
         ELOOPA: do i = 1, nalpha
            if(mclist(i,n) > nalpha) then
               do k = 1, nexalpha
                  if(virt(k,1) == mclist(i,n)) cycle ELOOPA
               enddo
               nexalpha = nexalpha + 1
               virt(nexalpha,1) = mclist(i,n)
            endif
         enddo ELOOPA
         ELOOPB: do i = 1, nbeta
            if(mclist(i+nalpha,n) > nbeta) then
               do k = 1, nexbeta
                  if(virt(k,2) == mclist(i+nalpha,n)) cycle ELOOPB
               enddo
               nexbeta = nexbeta + 1
               virt(nexbeta,2) = mclist(i+nalpha,n)
            endif
         enddo ELOOPB
      enddo

      allocate(exidx(max(nexalpha,nexbeta),2,ndet,2),nexc(ndet,2), stat=alstat)
      call assert(alstat==0, 'mdetinput: allocation3 failed')
      do n = 2, ndet
         exc = 0
         do j = 1, nalpha
            if(mclist(j, n) /= mclist(j, 1)) then
               exc = exc + 1
               exidx(exc,1,n,1) = mclist(j, 1)
               exidx(exc,2,n,1) = mclist(j, n)
            endif
         enddo
         nexc(n,1) = exc

         exc = 0
         do j = nalpha + 1, ne
            if(mclist(j, n) /= mclist(j,1)) then
               exc = exc + 1
               exidx(exc,1,n,2) = mclist(j, 1)
               exidx(exc,2,n,2) = mclist(j, n)
            endif
         enddo
         nexc(n,2) = exc
      enddo
   endif

   end subroutine reorder_dets



   subroutine mdetoutput(iu)
   !-----------------------!

! mdetoutput write CSFs to file unit 'iu'

   integer iu
   integer i,j,k,n

   if (ncsf==1 .and. ndets(1)==1 .and. mclist(nalpha,1)==nalpha &
      .and. mclist(ne,1)==ne) then
      write(iu,*) 'single unrestricted'
   else if (ncsf==1 .and. ndets(1)==1 .and. mclist(nalpha,1)==nalpha &
      .and. mclist(ne,1)==nbeta) then
      write(iu,*) 'single restricted'
   else
      write(iu,'(I5)') ncsf
      n = 0
      if (norb > 99) then
         do k=1,ncsf
            write(iu,'(g12.6,i5)') cci(k),ndets(k)
            do j=1,ndets(k)
               n = n+1
               write(iu,'(g12.6,1000(i4))') ccsf(j,k),(mclist(i,n),i=1,ne)
            enddo
         enddo
      else
         do k=1,ncsf
            write(iu,'(g12.6,i5)') cci(k),ndets(k)
            do j=1,ndets(k)
               n = n+1
               write(iu,'(g12.6,1000(i3))') ccsf(j,k),(mclist(i,n),i=1,ne)
            enddo
         enddo
      endif
      call assert(n==ndet,"mdetoutput: wrong det count")
   end if

   end subroutine mdetoutput


   subroutine mdetcalc(ie, nec, phi, fgrad, flapli, flapl, error_code)
   !-----------------------------------------------------------------!

! MDETCALC calculates the (multi-)determinantal part of the wavefunction,
! including the derivatives for a given configuration (x,y,z).
! All determinants with derivatives are kept and SAVEd for updating
! and faster evaluation of excited determinants.
! NEW: nec == 1 required. more than one electron config deprecated

   integer, intent(in)          :: ie       ! == 0 for all electron, > 0 for one electron calcn
   integer, intent(in)          :: nec      ! # of electron configurations
   real(r8), intent(inout)  :: phi(:)   ! Slater determinant part (for all elec configs)
   real(r8), intent(inout)  :: flapl(:),fgrad(:,:),flapli(:,:)  ! derivatives of phi: laplacian and gradient
   integer, intent(inout)         :: error_code
   integer w 

   call assert(size(fgrad,1)>=3*ne .and. size(flapli,1)>=ne,'mdetcalc: wrong size of argument (# elecs)')
   call assert(nec == 1 .and. size(phi) >= nec .and. nec <= mMOElecConfigs,  &
       'mdetcalc: incorrect size of argument (# elec configs)')

   error_code = MDET_NONE
   phi = 0
   flapl = 0
   fgrad = 0
   flapli = 0
   w = 1

   
   if (fastmdet .and. ie == 0) then
      call mdetcalcfast(phi, fgrad, flapli, flapl, w)
   else if (directDet .and. ie == 0) then
      call mdetDirectDet(phi, fgrad, flapli, flapl, w, error_code)
   else
      call mdetInvUpdateDet(ie, phi, fgrad, flapli, flapl, w, error_code)
   endif

   end subroutine mdetcalc




   subroutine mdetcalcfast(phi,fgrad,flapli,flapl,w)
   !-----------------------------------------------!
   ! uses fast matrix operations to calculate determinants
   ! only supports all electron move so far
   real(r8), intent(inout) :: phi(:)
   real(r8), intent(inout) :: flapl(:),fgrad(:,:),flapli(:,:)
   integer, intent(in) :: w

   integer :: i, ii, j, jj, k, n, na, nci, orb, ierr
   real(r8) :: d, tmp, f, g
   real(r8) :: dxa(nalpha), dya(nalpha), dza(nalpha), d2a(nalpha)
   real(r8) :: dxb(nbeta), dyb(nbeta), dzb(nbeta), d2b(nbeta)

   !!!!! not possible with ECP (why?). New assert!
   !!call assert(.not. use_ecp, "Cannot use fast mdet calculation with ECPs")

   ! Construct ground state matrix
   do i=1, nalpha
      do j=1, nalpha
         orb = mclist(j,1)
         deta(j,i,1) = mat(orb,i,w)
      enddo
   enddo

   do i = 1, nbeta
      ii = i + nalpha
      do j = 1, nbeta
         orb = mclist(nalpha + j,1)
         detb(j,i,1) = mat(orb,ii,w)
      enddo
   enddo

   call inv1(deta(1,1,1),nalpha,nalpha,detia(1,1,1),deter(1,1))
   call inv1(detb(1,1,1),nbeta ,nbeta ,detib(1,1,1),deter(2,1))

   do j = 1, nexalpha
      Ta(1:nalpha,virt(j,1)) = matmul(mat(virt(j,1),1:nalpha,w),detia(1:nalpha,1:nalpha,1))
   enddo
   do j = 1, nexbeta
      Tb(1:nbeta,virt(j,2)) = matmul(mat(virt(j,2),nalpha+1:ne,w),detib(1:nbeta,1:nbeta,1))
   enddo

   do nci = 2, ndet
      ii = nexc(nci, 1)
      do i = 1, ii
         do j = 1, ii
            excma(i, j) = Ta(exidx(i,1,nci,1), exidx(j,2,nci,1))
         enddo
      enddo
      call inv1(excma(:,:),ii,nexalpha,excwa,d)
      deter(1,nci) = deter(1,1) * d

      ii = nexc(nci, 2)
      do i = 1, ii
         do j = 1, ii
            excmb(i, j) = Tb(exidx(i,1,nci,2),exidx(j,2,nci,2))
         enddo
      enddo
      call inv1(excmb(:,:),ii,nexbeta,excwb,d)
      deter(2,nci) = deter(2,1) * d
   enddo

   dgrad(:,1,1)  = deter(1,1)
   dlapli(:,1,1) = deter(1,1)
   do i = 1, nalpha
      do j = 1, nalpha
         orb = mclist(j,1)
         dxa(j) = mat1x(orb,i,w)
         dya(j) = mat1y(orb,i,w)
         dza(j) = mat1z(orb,i,w)
         d2a(j) = mat2 (orb,i,w)
      enddo

      det1xia(:,:) = detia(:,:,1)
      det1yia(:,:) = detia(:,:,1)
      det1zia(:,:) = detia(:,:,1)
      det2ia(:,:)  = detia(:,:,1)
      call invcupd(dxa(:), i, nalpha, nalpha, det1xia(:,:), dgrad(3*i-2,1,1))
      call invcupd(dya(:), i, nalpha, nalpha, det1yia(:,:), dgrad(3*i-1,1,1))
      call invcupd(dza(:), i, nalpha, nalpha, det1zia(:,:), dgrad(3*i  ,1,1))
      call invcupd(d2a(:), i, nalpha, nalpha, det2ia(:,:),  dlapli(i   ,1,1))

      do j = 1, nexalpha
         jj = virt(j,1)

         d = mat(jj,i,w)
         mat(jj,i,w) = mat1x(jj,i,w)
         T1xa(1:nalpha,jj) = matmul(mat(jj,1:nalpha,w),det1xia(1:nalpha,1:nalpha))
         mat(jj,i,w) = mat1y(jj,i,w)
         T1ya(1:nalpha,jj) = matmul(mat(jj,1:nalpha,w),det1yia(1:nalpha,1:nalpha))
         mat(jj,i,w) = mat1z(jj,i,w)
         T1za(1:nalpha,jj) = matmul(mat(jj,1:nalpha,w),det1zia(1:nalpha,1:nalpha))
         mat(jj,i,w) = mat2(jj,i,w)
         T2a(1:nalpha,jj)  = matmul(mat(jj,1:nalpha,w),det2ia(1:nalpha,1:nalpha))
         mat(jj,i,w) = d
      enddo

      do nci = 2, ndet
         ii = nexc(nci,1)
         do k = 1, ii
            do n = 1, ii
               excm1xa(k, n) = T1xa(exidx(k,1,nci,1), exidx(n,2,nci,1))
               excm1ya(k, n) = T1ya(exidx(k,1,nci,1), exidx(n,2,nci,1))
               excm1za(k, n) = T1za(exidx(k,1,nci,1), exidx(n,2,nci,1))
               excm2a (k, n) = T2a (exidx(k,1,nci,1), exidx(n,2,nci,1))
            enddo
         enddo

         dgrad(3*i-2,1,nci) = fast_det(excm1xa,ii,nexalpha,excwa,dgrad(3*i-2,1,1))
         dgrad(3*i-1,1,nci) = fast_det(excm1ya,ii,nexalpha,excwa,dgrad(3*i-1,1,1))
         dgrad(3*i  ,1,nci) = fast_det(excm1za,ii,nexalpha,excwa,dgrad(3*i  ,1,1))
         dlapli(i   ,1,nci) = fast_det(excm2a, ii,nexalpha,excwa,dlapli(i   ,1,1))
      enddo
   enddo

   dgrad(:,2,1)  = deter(2,1)
   dlapli(:,2,1) = deter(2,1)
   do i = 1, nbeta
      ii = i + nalpha
      do j = 1, nbeta
         orb = mclist(nalpha + j,1)
         dxb(j) = mat1x(orb,ii,w)
         dyb(j) = mat1y(orb,ii,w)
         dzb(j) = mat1z(orb,ii,w)
         d2b(j)  = mat2(orb,ii,w)
      enddo

      det1xib(:,:) = detib(:,:,1)
      det1yib(:,:) = detib(:,:,1)
      det1zib(:,:) = detib(:,:,1)
      det2ib(:,:)  = detib(:,:,1)
      call invcupd(dxb(:), i, nbeta, nbeta, det1xib(:,:), dgrad(3*i-2,2,1))
      call invcupd(dyb(:), i, nbeta, nbeta, det1yib(:,:), dgrad(3*i-1,2,1))
      call invcupd(dzb(:), i, nbeta, nbeta, det1zib(:,:), dgrad(3*i  ,2,1))
      call invcupd(d2b(:), i, nbeta, nbeta, det2ib(:,:),  dlapli(i   ,2,1))

      do j = 1, nexbeta
         d = mat(virt(j,2),ii,w)
         mat(virt(j,2),ii,w) = mat1x(virt(j,2),ii,w)
         T1xb(1:nbeta,virt(j,2)) = matmul(mat(virt(j,2),nalpha+1:ne,w),det1xib(1:nbeta,1:nbeta))
         mat(virt(j,2),ii,w) = mat1y(virt(j,2),ii,w)
         T1yb(1:nbeta,virt(j,2)) = matmul(mat(virt(j,2),nalpha+1:ne,w),det1yib(1:nbeta,1:nbeta))
         mat(virt(j,2),ii,w) = mat1z(virt(j,2),ii,w)
         T1zb(1:nbeta,virt(j,2)) = matmul(mat(virt(j,2),nalpha+1:ne,w),det1zib(1:nbeta,1:nbeta))
         mat(virt(j,2),ii,w) = mat2(virt(j,2),ii,w)
         T2b(1:nbeta,virt(j,2))  = matmul(mat(virt(j,2),nalpha+1:ne,w),det2ib(1:nbeta,1:nbeta))
         mat(virt(j,2),ii,w) = d
      enddo

      do nci = 2, ndet
         ii = nexc(nci,2)
         do k = 1, ii
            do n = 1, ii
               excm1xb(k, n) = T1xb(exidx(k,1,nci,2), exidx(n,2,nci,2))
               excm1yb(k, n) = T1yb(exidx(k,1,nci,2), exidx(n,2,nci,2))
               excm1zb(k, n) = T1zb(exidx(k,1,nci,2), exidx(n,2,nci,2))
               excm2b (k, n) = T2b (exidx(k,1,nci,2), exidx(n,2,nci,2))
            enddo
         enddo

         dgrad(3*i-2,2,nci) = fast_det(excm1xb,ii,nexbeta,excwb,dgrad(3*i-2,2,1))
         dgrad(3*i-1,2,nci) = fast_det(excm1yb,ii,nexbeta,excwb,dgrad(3*i-1,2,1))
         dgrad(3*i  ,2,nci) = fast_det(excm1zb,ii,nexbeta,excwb,dgrad(3*i  ,2,1))
         dlapli(i   ,2,nci) = fast_det(excm2b, ii,nexbeta,excwb,dlapli(i   ,2,1))
      enddo
   enddo

   phi(w) = 0d0
   fgrad(1:3*ne,w) = 0d0
   flapli(1:ne,w) = 0d0

   n = 0
   na = 3*nalpha
   do k=1,ncsf
      do j=1,ndets(k)
         n = n+1
         tmp = cci(k)*ccsf(j,k)
         phi(w) = phi(w) + tmp*deter(1,n)*deter(2,n)
         do i=1,na
            fgrad(i,w) = fgrad(i,w) + tmp*dgrad(i,1,n)*deter(2,n)
         enddo
         do i=1,3*nbeta
            fgrad(na+i,w) = fgrad(na+i,w) + tmp*deter(1,n)*dgrad(i,2,n)
         enddo
         do i=1,nalpha
            flapli(i,w) = flapli(i,w) + tmp*dlapli(i,1,n)*deter(2,n)
         enddo
         do i=1,nbeta
            flapli(nalpha+i,w) = flapli(nalpha+i,w) + tmp*deter(1,n)*dlapli(i,2,n)
         enddo
      enddo
   enddo

   ! Sum individual laplacians to total laplacian
   flapl(w) = sum(flapli(1:ne,w))

   end subroutine mdetcalcfast

   !-----------------------------------!
   !DEC$ ATTRIBUTES INLINE :: fast_det
   real(r8) function fast_det(a,n,nmax,a1,det)
   !-----------------------------------!
      implicit none
      integer, intent(in) :: n,nmax
      real(r8), intent(in) :: a(nmax,nmax), det
      real(r8), intent(inout) :: a1(nmax,nmax)
      real(r8) :: d

      if(n == 1) then
         fast_det = det * a(1,1)
      else if(n == 2) then
         fast_det = det * (a(1,1) * a(2,2) - a(1,2) * a(2,1))
      else
         call inv1(a,n,nmax,a1,d)
         fast_det = det * d
      endif
   end function fast_det



   subroutine mdetDirectDet(phi, fgrad, flapli, flapl, w, error_code)
   !----------------------------------------------------------------!
   ! uses LU decomp to calculate all determinants incl derivatives
   ! for debugging only !
   real(r8), intent(inout) :: phi(:)
   real(r8), intent(inout) :: flapl(:),fgrad(:,:),flapli(:,:)
   integer, intent(in)     :: w
   integer, intent(inout)  :: error_code

   real(r8) :: cola(size(deta,1)), colb(size(detb,1))
   real(r8) :: d, tmp
   integer :: i, ii, j, k, n, na, nci, orb, error_this_det, offset


   offset = 0

   ! CI loop over products of determinants
   CILOOPA: do nci=1,ndet

   error_this_det = MDET_NONE

   if (detsRepLst(nci, 1) == 0) then !  only for new dets
      ! Construct matrices for determinants
      do i = 1, nalpha
         ii = i + offset
         do j = 1, nalpha
            orb = mclist(j + offset,nci)
            deta(j,i,nci)   = mat(orb,ii,w)
            det1xa(j,i,nci) = mat1x(orb,ii,w)
            det1ya(j,i,nci) = mat1y(orb,ii,w)
            det1za(j,i,nci) = mat1z(orb,ii,w)
            det2a(j,i,nci)  = mat2(orb,ii,w)
         end do
      end do

      cola = 0
      call lapack_det(nalpha, deta(:,:,nci), cola, 0, d, error_this_det)
      deter(1, nci) = d
      do i = 1, nalpha
         cola = det1xa(:, i, nci)
         call lapack_det(nalpha, deta(:,:,nci), cola, i, d, error_this_det)
         dgrad(3*i-2, 1, nci) = d
         cola = det1ya(:, i, nci)
         call lapack_det(nalpha, deta(:,:,nci), cola, i, d, error_this_det)
         dgrad(3*i-1, 1, nci) = d
         cola = det1za(:, i, nci)
         call lapack_det(nalpha, deta(:,:,nci), cola, i, d, error_this_det)
         dgrad(3*i, 1, nci) = d
         cola = det2a(:, i, nci)
         call lapack_det(nalpha, deta(:,:,nci), cola, i, d, error_this_det)
         dlapli(i, 1, nci) = d               
      end do

   else !detsRepLst already calculated

      do i = 1, nalpha            ! inverse deti is changed on ALL positions
         dgrad(3*i-2,1,nci) = dgrad(3*i-2,1,detsRepLst(nci,1))
         dgrad(3*i-1,1,nci) = dgrad(3*i-1,1,detsRepLst(nci,1))
         dgrad(3*i,1,nci)   = dgrad(3*i,1,detsRepLst(nci,1))
         dlapli(i,1,nci)    = dlapli(i,1,detsRepLst(nci,1))
      enddo
      deter(1,nci) = deter(1,detsRepLst(nci,1))

   end if !detsRepLst

   end do CILOOPA

   ! Now the same with beta!
   offset = nalpha

   ! CI loop over products of determinants
   CILOOPB: do nci=1,ndet

   error_this_det = MDET_NONE

   if (detsRepLst(nci, 2) == 0) then ! only for new dets
      ! Construct matrices for determinants
      do i = 1, nbeta
         ii = i + offset
         do j = 1, nbeta
            orb = mclist(j + offset,nci)
            detb(j,i,nci)   = mat(orb,ii,w)
            det1xb(j,i,nci) = mat1x(orb,ii,w)
            det1yb(j,i,nci) = mat1y(orb,ii,w)
            det1zb(j,i,nci) = mat1z(orb,ii,w)
            det2b(j,i,nci)  = mat2(orb,ii,w)
         enddo
      enddo

      colb = 0
      call lapack_det(nbeta, detb(:,:,nci), colb, 0, d, error_this_det)
      deter(2, nci) = d
      do i = 1, nbeta
         colb = det1xb(:, i, nci)
         call lapack_det(nbeta, detb(:,:,nci), colb, i, d, error_this_det)
         dgrad(3*i-2, 2, nci) = d
         colb = det1yb(:, i, nci)
         call lapack_det(nbeta, detb(:,:,nci), colb, i, d, error_this_det)
         dgrad(3*i-1, 2, nci) = d
         colb = det1zb(:, i, nci)
         call lapack_det(nbeta, detb(:,:,nci), colb, i, d, error_this_det)
         dgrad(3*i, 2, nci) = d
         colb = det2b(:, i, nci)
         call lapack_det(nbeta, detb(:,:,nci), colb, i, d, error_this_det)
         dlapli(i, 2, nci) = d
      end do

   else   !detsRepLst. already calculated

      do i=1,nbeta            ! inverse deti is changed on ALL positions

         dgrad(3*i-2,2,nci) = dgrad(3*i-2,2,detsRepLst(nci,2))
         dgrad(3*i-1,2,nci) = dgrad(3*i-1,2,detsRepLst(nci,2))
         dgrad(3*i,2,nci)   = dgrad(3*i,2,detsRepLst(nci,2))
         dlapli(i,2,nci)    = dlapli(i,2,detsRepLst(nci,2))
      enddo
      deter(2,nci) = deter(2,detsRepLst(nci,2))

   endif  !detsRepLst

   enddo CILOOPB

   ! Get orbital part phi and its derivatives
   phi(w) = 0.0_r8
   fgrad(1:3*ne,w) = 0.0_r8
   flapli(1:ne,w) = 0.0_r8

   n = 0
   na = 3*nalpha
   if (nbeta > 0) then
      do k=1,ncsf
         do j=1,ndets(k)
            n = n+1
            tmp = cci(k)*ccsf(j,k)
            phi(w) = phi(w) + tmp*deter(1,n)*deter(2,n)
            do i=1,na
               fgrad(i,w) = fgrad(i,w) + tmp*dgrad(i,1,n)*deter(2,n)
            enddo
            do i=1,3*nbeta
               fgrad(na+i,w) = fgrad(na+i,w) + tmp*deter(1,n)*dgrad(i,2,n)
            enddo
            do i=1,nalpha
               flapli(i,w) = flapli(i,w) + tmp*dlapli(i,1,n)*deter(2,n)
            enddo
            do i=1,nbeta
               flapli(nalpha+i,w) = flapli(nalpha+i,w) + tmp*deter(1,n)*dlapli(i,2,n)
            enddo
         enddo
      enddo
   else    ! no beta electrons => no beta dets
      do k=1,ncsf
         do j=1,ndets(k)
            n = n+1
            tmp = cci(k)*ccsf(j,k)
            phi(w) = phi(w) + tmp*deter(1,n)
            do i=1,na
               fgrad(i,w) = fgrad(i,w) + tmp*dgrad(i,1,n)
            enddo
            do i=1,nalpha
               flapli(i,w) = flapli(i,w) + tmp*dlapli(i,1,n)
            enddo
         enddo
      enddo
   endif

   ! Sum individual laplacians to total laplacian
   flapl(w) = sum(flapli(1:ne,w))

   end subroutine mdetDirectDet




   subroutine mdetInvUpdateDet(ie, phi, fgrad, flapli, flapl, w, error_code)
   !-------------------------------------------------------------------------!

! MDETCALC calculates the (multi-)determinantal part of the wavefunction,
! including the derivatives for a given configuration (x,y,z).
! All determinants with derivatives are kept and SAVEd for updating
! and faster evaluation of excited determinants.

   integer, intent(in)          :: ie       ! == 0 for all electron, > 0 for one electron calcn
   real(r8), intent(inout)  :: phi(:)   ! Slater determinant part (for all elec configs)
   real(r8), intent(inout)  :: flapl(:),fgrad(:,:),flapli(:,:)  ! derivatives of phi: laplacian and gradient
   integer, intent(in)          :: w     ! deprecated electron config index
   integer, intent(out)         :: error_code 
   integer i,j,k,l,orb,ii,ii1,ii2
   integer nci,nnci,n,na,offset,ierr, error_this_det
   real(r8) tmp,tmp1,tmp2,tmp3,tmp4, rerr
   real(r8) :: cola(size(deta,1)), colb(size(detb,1))
   real(r8) d
   logical :: update 
#ifdef CHKNANUP
   real(r8) :: switchDirectThreshold, absDetFirst, detInverseThreshold
   logical :: checkMDetError
#endif

   error_code = MDET_NONE
   phi = 0
   flapl = 0
   fgrad = 0
   flapli = 0

#ifdef CHKNANUP
   switchDirectThreshold = getSwitchDirectThreshold()
   detInverseThreshold = getDetInverseThreshold()
#endif

   ! Note on algorithm:
   ! First ALL alpha dets and then ALL beta dets are calculated
   ! (but only alpha or beta in one electron update case: ie > 0)
   ! This allows updating the inverse matrix from det to det
   ! alpha and beta arrays are separate to allow for spin polarized systems
   ! while keeping memory contiguous

   offset = 0
   ii1 = 0
   ii2 = 0
   if (ie == 0) then
      ii1  = 1
      ii2  = nalpha
   else if(ie <= nalpha) then
      ii1 = ie
      ii2 = ie
   endif

   ! if this is a OEM call, and the electron is a beta electron, skip over the
   ! alpha loop
   if(ii1 /= 0) then

   ! CI loop over products of determinants
   CILOOPA: do nci = 1, ndet

   error_this_det = MDET_NONE

   if (detsRepLst(nci,1)==0) then !  only for new dets
      ! Construct matrices for determinants
      do i = ii1, ii2
         ii = i + offset
         do j = 1, nalpha
            orb = mclist(j + offset,nci)
            deta(j,i,nci)   = mat(orb,ii,w)
            det1xa(j,i,nci) = mat1x(orb,ii,w)
            det1ya(j,i,nci) = mat1y(orb,ii,w)
            det1za(j,i,nci) = mat1z(orb,ii,w)
            det2a(j,i,nci)  = mat2(orb,ii,w)
         enddo
      enddo

      ! Evaluation of the determinants
      if (nalpha == 1) then

         deter(1,nci)  =deta(1,1,nci)
         detia(1,1,nci) = 1 / deta(1,1,nci)
         dgrad(1,1,nci)=det1xa(1,1,nci)
         dgrad(2,1,nci)=det1ya(1,1,nci)
         dgrad(3,1,nci)=det1za(1,1,nci)
         dlapli(1,1,nci)=det2a(1,1,nci)

      else if (nalpha > 1) then

         if (ie == 0) then               ! all electrons new

            if (nci == 1) then

               ! calculate inverse matrix of det (LU decomposition) (is O(n**3) )
               ! returns inverse matrix in deti
               ! Note: ierr==MDET_LU_ZERO_DET means det==0, singular matrix
               ! det, grad, and lapl exist, but require direct calculation of grad and lapl
               ! TODO, but see excited dets below 
               call lapack_inv(nalpha, deta(:,:,nci), detia(:,:,nci), d, ierr)
               if (ierr > 0) then
                  error_code = ierr
                  if (logmode > 3) write(iull,*) "mdetInvUpdateDet: LU error 1st det (alpha)", ierr
                  goto 999
               end if
               do nnci = 2, ndet
                  detia(:nalpha,:nalpha,nnci) = detia(:nalpha,:nalpha,nci)
                  deter(1, nnci) = d
               end do
#ifdef CHKNANUP
               absDetFirst = abs(d)
#endif
            else      ! nci > 1, excited dets
               ! update inverse matrix for all rows (MOs) that are "excited"
               d = deter(1, nci)
#ifdef CHKNANUP
               update = .true.
#endif
               do j = 1, nalpha
                  if (mclist(j+offset,nci) /= mclist(j+offset,1)) then
                     call invrupd(deta(j,:,nci), j, nalpha, nalpha, detia(:,:,nci), d)
#ifdef CHKNANUP
                     if (.not. (ieee_is_normal(d) .and. ALL(ieee_is_normal(detia(:,:,nci)))) &
                        .or. abs(d)/absDetFirst < switchDirectThreshold) then
                        update = .false.
                        exit
                     end if
#endif
                  end if
               end do
               !!! if a matrix during update becomes singular calculate without updates
#ifdef CHKNANUP
               if (.not.update) then

                  call lapack_inv(nalpha, deta(:,:,nci), detia(:,:,nci), d, ierr)

                  if (ierr > 0) then
                     error_this_det = ierr
                     if (error_this_det /= MDET_LU_ZERO_DET) then
                        error_code = error_this_det
                        if (logmode > 3) write(iull,*) "mdetInvUpdateDet: LU error excited det (alpha)"
                        goto 999 ! leave mdet
                     end if
                  end if
               end if
#endif
            end if

         else ! ie > 0, one electron update

            ! update inverse matrix deti of det
            d = deter(1,nci)
            ! first arg is "pointer" to ii1-th column, doesn't work for rows
            call invcupd(deta(1,ii1,nci), ii1, nalpha, nalpha, detia(1,1,nci), d)

            !!! CAREFUL: CHECK SINGULARITIES HERE LIKE IN ALL ELECTRON MOVES

         end if

         ! construct gradients and laplacians of det as scalar prod with deti
         ! this is O(n**2)
#ifdef CHKNANUP
         if (error_this_det == MDET_LU_ZERO_DET .or. abs(d)/absDetFirst < detInverseThreshold) then
    
            deter(1, nci) = d  
            ! direct LU decomposition necessary
            do i = 1, nalpha
               cola = det1xa(:, i, nci)
               call lapack_det(nalpha, deta(:,:,nci), cola, i, d, ierr)
               dgrad(3*i-2, 1, nci) = d
               cola = det1ya(:, i, nci)
               call lapack_det(nalpha, deta(:,:,nci), cola, i, d, ierr)
               dgrad(3*i-1, 1, nci) = d
               cola = det1za(:, i, nci)
               call lapack_det(nalpha, deta(:,:,nci), cola, i, d, ierr)
               dgrad(3*i, 1, nci) = d
               cola = det2a(:, i, nci)
               call lapack_det(nalpha, deta(:,:,nci), cola, i, d, ierr)
               dlapli(i, 1, nci) = d               
            end do

         else
#endif

            do i = 1, nalpha            ! inverse deti is changed on ALL positions
               tmp1 = dot_product(det1xa(:,i,nci),detia(i,:,nci))
               tmp2 = dot_product(det1ya(:,i,nci),detia(i,:,nci))
               tmp3 = dot_product(det1za(:,i,nci),detia(i,:,nci))
               tmp4 = dot_product(det2a(:,i,nci),detia(i,:,nci))
               dgrad(3*i-2,1,nci) = tmp1*d
               dgrad(3*i-1,1,nci) = tmp2*d
               dgrad(3*i,1,nci) = tmp3*d
               dlapli(i,1,nci) =  tmp4*d
            end do
            deter(1, nci) = d

#ifdef CHKNANUP
         end if
#endif

      endif ! nalpha==1

   else !detsRepLst already calculated

      do i = 1, nalpha            ! inverse deti is changed on ALL positions
         dgrad(3*i-2,1,nci) = dgrad(3*i-2,1,detsRepLst(nci,1))
         dgrad(3*i-1,1,nci) = dgrad(3*i-1,1,detsRepLst(nci,1))
         dgrad(3*i,1,nci)   = dgrad(3*i,1,detsRepLst(nci,1))
         dlapli(i,1,nci)    = dlapli(i,1,detsRepLst(nci,1))
      end do
      deter(1, nci) = deter(1, detsRepLst(nci,1))

   end if !detsRepLst

   end do CILOOPA
   
   end if ! ii1 /= 0


   ! Now the same with beta!
   offset = nalpha
   ii1 = 0
   ii2 = 0
   if (ie == 0) then
      ii1  = 1
      ii2  = nbeta
   else if(ie > nalpha) then
      ii1 = ie - offset
      ii2 = ie - offset
   endif

   if(ii1 /= 0) then

   ! CI loop over products of determinants
   CILOOPB: do nci=1,ndet

   error_this_det = MDET_NONE

   if (detsRepLst(nci,2)==0) then ! only for new dets
      ! Construct matrices for determinants
      do i = ii1, ii2
         ii = i + offset
         do j = 1, nbeta
            orb = mclist(j + offset,nci)
            detb(j,i,nci)   = mat(orb,ii,w)
            det1xb(j,i,nci) = mat1x(orb,ii,w)
            det1yb(j,i,nci) = mat1y(orb,ii,w)
            det1zb(j,i,nci) = mat1z(orb,ii,w)
            det2b(j,i,nci)  = mat2(orb,ii,w)
         enddo
      enddo

      ! Evaluation of the determinants
      if (nbeta==1) then

         deter(2,nci)  =detb(1,1,nci)
         detib(1,1,nci) = 1 / detb(1,1,nci)
         dgrad(1,2,nci)=det1xb(1,1,nci)
         dgrad(2,2,nci)=det1yb(1,1,nci)
         dgrad(3,2,nci)=det1zb(1,1,nci)
         dlapli(1,2,nci)=det2b(1,1,nci)

      else if (nbeta > 1) then

         if (ie == 0) then               ! all electrons new

            if (nci == 1) then

               ! calculate inverse matrix of det (LU decomposition) (is O(n**3) )
               ! returns inverse matrix in deti
               call lapack_inv(nbeta, detb(:,:,nci), detib(:,:,nci), d, ierr)
               if (ierr > 0) then                        
                  error_code = ierr
                  if (logmode > 3) write(iull,*) "mdetInvUpdateDet: LU error 1st det (beta)", ierr
                  goto 999
               end if
               do nnci = 2, ndet
                  detib(:nbeta,:nbeta,nnci) = detib(:nbeta,:nbeta,nci)
                  deter(2,nnci) = d
               end do
#ifdef CHKNANUP
               absDetFirst = abs(d)
#endif
            else      ! nci > 1, excited dets
               ! update inverse matrix for all rows (MOs) that are "excited"
               d = deter(2, nci)
#ifdef CHKNANUP
               update = .true.
#endif
               do j = 1, nbeta
                  if (mclist(j+offset,nci) /= mclist(j+offset,1)) then
                     call invrupd(detb(j,:,nci), j, nbeta, nbeta, detib(:,:,nci), d)
#ifdef CHKNANUP
                     if (.not. (ieee_is_normal(d) .and. ALL(ieee_is_normal(detia(:,:,nci)))) &
                        .or. abs(d)/absDetFirst < switchDirectThreshold) then
                        update = .false.
                        exit
                     end if
#endif
                  end if
               end do
               !!! if a matrix during update is singular calculate without updates
#ifdef CHKNANUP
               if (.not.update) then

                  call lapack_inv(nbeta, detb(:,:,nci), detib(:,:,nci), d, ierr)

                  if (ierr > 0) then
                     error_this_det = ierr
                     if (error_this_det /= MDET_LU_ZERO_DET) then
                        error_code = error_this_det
                        if (logmode > 3) write(iull,*) "mdetInvUpdateDet: LU error excited det (beta)"
                        goto 999 ! leave mdet
                     end if
                  end if
               end if
#endif
            endif

         else ! ie > 0, one electron update

            ! update inverse matrix deti of det
            d = deter(2,nci)
            ! first arg is "pointer" to ii1-th column, doesn't work for rows
            call invcupd(detb(1,ii1,nci), ii1, nbeta, nbeta, detib(1,1,nci), d)

            !!! CAREFUL: CHECK SINGULARITIES HERE LIKE IN ALL ELECTRON MOVES

         endif

         ! construct gradients and laplacians of det as scalar prod with deti
         ! this is O(n**2)
#ifdef CHKNANUP
         if (error_this_det == MDET_LU_ZERO_DET .or. abs(d)/absDetFirst < detInverseThreshold) then
    
            deter(2, nci) = d
            ! direct LU decomposition necessary            
            do i = 1, nbeta
               colb = det1xb(:, i, nci)
               call lapack_det(nbeta, detb(:,:,nci), colb, i, d, ierr)
               dgrad(3*i-2, 2, nci) = d
               colb = det1yb(:, i, nci)
               call lapack_det(nbeta, detb(:,:,nci), colb, i, d, ierr)
               dgrad(3*i-1, 2, nci) = d
               colb = det1zb(:, i, nci)
               call lapack_det(nbeta, detb(:,:,nci), colb, i, d, ierr)
               dgrad(3*i, 2, nci) = d
               colb = det2b(:, i, nci)
               call lapack_det(nbeta, detb(:,:,nci), colb, i, d, ierr)
               dlapli(i, 2, nci) = d
            end do

         else
#endif
            do i = 1, nbeta            ! inverse deti is changed on ALL positions
               tmp1 = dot_product(det1xb(:,i,nci),detib(i,:,nci))
               tmp2 = dot_product(det1yb(:,i,nci),detib(i,:,nci))
               tmp3 = dot_product(det1zb(:,i,nci),detib(i,:,nci))
               tmp4 = dot_product(det2b(:,i,nci),detib(i,:,nci))
               dgrad(3*i-2, 2, nci) = tmp1*d
               dgrad(3*i-1, 2, nci) = tmp2*d
               dgrad(3*i, 2, nci) = tmp3*d
               dlapli(i, 2, nci) =  tmp4*d
            enddo
            deter(2, nci) = d
#ifdef CHKNANUP
         end if
#endif

      end if ! nbeta == 1

   else   !detsRepLst. already calculated

      do i = 1, nbeta            ! inverse deti is changed on ALL positions
         dgrad(3*i-2,2,nci) = dgrad(3*i-2,2,detsRepLst(nci,2))
         dgrad(3*i-1,2,nci) = dgrad(3*i-1,2,detsRepLst(nci,2))
         dgrad(3*i,2,nci)   = dgrad(3*i,2,detsRepLst(nci,2))
         dlapli(i,2,nci)    = dlapli(i,2,detsRepLst(nci,2))
      end do
      deter(2,nci) = deter(2,detsRepLst(nci,2))

   end if  !detsRepLst

   end do CILOOPB

   end if ! ii1 /= 0



   ! Get orbital part phi and its derivatives
   phi(w) = 0d0
   fgrad(1:3*ne,w) = 0d0
   flapli(1:ne,w) = 0d0

   n = 0
   na = 3*nalpha
   if (nbeta > 0) then
      do k=1,ncsf
         do j=1,ndets(k)
            n = n+1
            tmp = cci(k)*ccsf(j,k)
            phi(w) = phi(w) + tmp*deter(1,n)*deter(2,n)
            do i=1,na
               fgrad(i,w) = fgrad(i,w) + tmp*dgrad(i,1,n)*deter(2,n)
            enddo
            do i=1,3*nbeta
               fgrad(na+i,w) = fgrad(na+i,w) + tmp*deter(1,n)*dgrad(i,2,n)
            enddo
            do i=1,nalpha
               flapli(i,w) = flapli(i,w) + tmp*dlapli(i,1,n)*deter(2,n)
            enddo
            do i=1,nbeta
               flapli(nalpha+i,w) = flapli(nalpha+i,w) + tmp*deter(1,n)*dlapli(i,2,n)
            enddo
         enddo
      enddo
   else    ! no beta electrons => no beta dets
      do k=1,ncsf
         do j=1,ndets(k)
            n = n+1
            tmp = cci(k)*ccsf(j,k)
            phi(w) = phi(w) + tmp*deter(1,n)
            do i=1,na
               fgrad(i,w) = fgrad(i,w) + tmp*dgrad(i,1,n)
            enddo
            do i=1,nalpha
               flapli(i,w) = flapli(i,w) + tmp*dlapli(i,1,n)
            enddo
         enddo
      enddo
   endif

   ! Sum individual laplacians to total laplacian
   flapl(w) = sum(flapli(1:ne,w))

   999 continue

   end subroutine mdetInvUpdateDet



   !--------------------------------!
   subroutine mdet1calc(ie,phi,iflag)
   !--------------------------------!

! MDET1CALC calculates the (multi-)determinantal part of the wavefunction,
! WITHOUT the derivatives for a given configuration (x,y,z).
! All determinants with derivatives are kept and SAVEd for updating
! and faster evaluation of excited determinants.

   integer, intent(in)          :: ie       ! == 0 for all electron, > 0 for one electron calcn
   real(r8), intent(inout)        :: phi      ! Slater determinant part
   integer, intent(out)         :: iflag    ! 0: OK, 1: sing matrix, 2: error in dgetri

   integer offset
   integer i,j,k,orb,ii,ii1,ii2
   integer nci,nnci,n,na
   integer tmp,ierr
   real(r8) d

   ! Note on algorithm:
   ! First ALL alpha dets and then ALL beta dets are calculated
   ! (but only alpha or beta in one electron update case: ie > 0)
   ! This allows updating the inverse matrix from det to det
   ! alpha and beta arrays are separate to allow for spin polarized systems
   ! while keeping memory contiguous

   iflag = 0
   phi = 0
   offset = 0
   ii1 = 0
   ii2 = 0
   if (ie == 0) then
      ii1  = 1
      ii2  = nalpha
   else if(ie <= nalpha) then
      ii1 = ie
      ii2 = ie
   endif

   if(ii1 /= 0) then
   ! CI loop over products of determinants
   CILOOPA: do nci=1,ndet

   ! Construct matrices for determinants
   do i=ii1,ii2
      ii = i + offset
      do j=1,nalpha
         orb = mclist(j + offset,nci)
         deta(j,i,nci)   = mat(orb,ii,1)
      enddo
   enddo

   ! Evaluation of the determinants
   if (nalpha==1) then
      deter(1,nci)  =deta(1,1,nci)
   else
      if (ie == 0) then               ! all electrons new
         if (nci == 1) then

            ! calculate inverse matrix of det (LU decomposition) (is O(n**3) )
            ! returns inverse matrix in deti
            call lapack_inv(nalpha,deta(:,:,nci),detia(:,:,nci),d,ierr)
            if (ierr > 0) iflag = ierr
            goto 999
            do nnci=nci+1,ndet
               detia(:nalpha,:nalpha,nnci) = detia(:nalpha,:nalpha,nci)
               deter(1,nnci) = d
            enddo
         else      ! nci > 1, excited dets
            ! update inverse matrix for all rows (MOs) that are "excited"
            d = deter(1,nci)
            do j=1,nalpha
               if (mclist(j+offset,nci) /= mclist(j+offset,1)) then
                  call invrupd(deta(j,:,nci),j,nalpha,nalpha,detia(1,1,nci),d)
               endif
            enddo
         endif

      else ! ie > 0, one electron update

         ! update inverse matrix deti of det
         d = deter(1,nci)
         ! first arg is "pointer" to ii1-th column, doesn't work for rows
         call invcupd(deta(1,ii1,nci),ii1,nalpha,nalpha,detia(1,1,nci),d)
      endif

      deter(1,nci) = d

   endif ! nalpha==1

   enddo CILOOPA
   endif ! ii1 /= 0

   ! Now the same with beta!
   offset = nalpha
   ii1 = 0
   ii2 = 0
   if (ie == 0) then
      ii1  = 1
      ii2  = nbeta
   else if(ie > nalpha) then
      ii1 = ie - offset
      ii2 = ie - offset
   endif

   if(ii1 /= 0) then
   ! CI loop over products of determinants
   CILOOPB: do nci=1,ndet

   ! Construct matrices for determinants
   do i=ii1,ii2
      ii = i + offset
      do j=1,nbeta
         orb = mclist(j + offset,nci)
         detb(j,i,nci)   = mat(orb,ii,1)
      enddo
   enddo

   ! Evaluation of the determinants
   if (nbeta==1) then
      deter(2,nci)  =detb(1,1,nci)
   else
      if (ie == 0) then               ! all electrons new
         if (nci == 1) then

            ! calculate inverse matrix of det (LU decomposition) (is O(n**3) )
            ! returns inverse matrix in deti
            call lapack_inv(nbeta,detb(:,:,nci),detib(:,:,nci),d,ierr)
            if (ierr > 0) iflag = ierr
            goto 999
            do nnci=nci+1,ndet
               detib(:nbeta,:nbeta,nnci) = detib(:nbeta,:nbeta,nci)
               deter(2,nnci) = d
            enddo
         else      ! nci > 1, excited dets
            ! update inverse matrix for all rows (MOs) that are "excited"
            d = deter(2,nci)
            do j=1,nbeta
               if (mclist(j+offset,nci) /= mclist(j+offset,1)) then
                  call invrupd(detb(j,:,nci),j,nbeta,nbeta,detib(1,1,nci),d)
               endif
            enddo
         endif

      else ! ie > 0, one electron update

         ! update inverse matrix deti of det
         d = deter(2,nci)
         ! first arg is "pointer" to ii1-th column, doesn't work for rows
         call invcupd(detb(1,ii1,nci),ii1,nbeta,nbeta,detib(1,1,nci),d)
      endif

      deter(2,nci) = d

   endif ! nbeta == 1

   enddo CILOOPB
   endif ! ii1 /= 0


   ! Get orbital part phi
   phi = 0d0

   n = 0
   na = 3*nalpha
   do k=1,ncsf
      do j=1,ndets(k)
         n = n+1
         tmp = cci(k)*ccsf(j,k)
         phi = phi + tmp*deter(1,n)*deter(2,n)
      enddo
   enddo

   999 continue

   end subroutine mdet1calc

   subroutine mdet2calc(ie,phi,fk,ci_param_mode,resetDet,moUpdateMode)
   !----------------------------------------------------------!
   ! calculates the (multi)determinant part of the wave function for a potential
   ! one electron move. note that no matrices are updated by this function, thus
   ! this should only be used when the proposed electron move is not actually
   ! (permanently) executed. this function is mainly used in the Lebedev integration
   ! of the pseudopotential localization. only supports one electron moves.

   integer, intent(in)          :: ie    ! electron number, required to be != 0
   real(r8), intent(out)          :: phi   ! Slater determinant part
   integer, intent(in),optional :: ci_param_mode !if we need fk we need this
   real(r8), intent(out),optional :: fk(:) ! divative are needed only for vNlk
   logical,intent(in),optional  :: resetDet ! necessary for mo optimization
   integer,intent(in),optional  :: moUpdateMode ! for mok
   integer :: i,ii,j,k,n,offset
   integer :: nci
   real(r8) :: tmp
   real(r8) :: olda(1:nalpha), oldb(1:nbeta)
   integer, parameter :: SIMPLE=1, UPDATE=2 !for mok

   if(present(moUpdateMode)) then
      if(moUpdateMode==UPDATE) then
         call internal_update_det_and_inv()
       else
          call internal_update_det()
       endif
    else
      call internal_update_det()
   endif


   ! Get orbital part phi
   phi = 0d0

   n = 0
   do k=1,ncsf
      do j=1,ndets(k)
         n = n+1
         tmp = cci(k)*ccsf(j,k)
         phi = phi + tmp*deter(1,n)*deter(2,n)
      enddo
   enddo
   if(present(fk)) then ! calculate drivatives for vNlk
           n = 0
           if (ci_param_mode == 3) then
              offset = 0
           else
              offset = 1
              n = n + ndets(1)
           endif
           do k=1,ncsf-offset
              fk(k) = 0.d0
              do j=1,ndets(k+offset)
                 n = n+1
                 tmp = ccsf(j,k+offset)
                 fk(k) = fk(k) + tmp*deter(1,n)*deter(2,n)
              enddo
           enddo
           call assert(n==ndet,'mdet2cal: internal error')
  endif

   if( .not. present(resetDet)) then
          call resetToOld(ie)
   elseif ( resetDet) then
          call resetToOld(ie)
   endif
    contains
      subroutine internal_update_det() ! only updates determinant
           if(ie <= nalpha) then
               olddet(1:ndet) = deter(1,:)
               ! CI loop over products of determinants
               CILOOPA: do nci=1,ndet
                if (detsRepLst(nci,1)==0) then  ! calculation only for new dets
                  i = ie
                  ii = ie
                  !olda(1:nalpha) = deta(:,i,nci) ! please take a look this line might be unnecessary!
                  ! Construct matrices for determinants
                  deta(1:nalpha, i, nci) = mat(mclist(1:nalpha,nci), ii, 1)

                  if (nalpha==1) then
                     deter(1,nci) = deta(1,1,nci)
                  else
                     ! the following call is equivalent to
                     ! deter(1,nci) = deter(1,nci) * dot_product(deta(1:nalpha,i,nci),detia(i,1:nalpha,nci))
                     ! but depending on the user settings uses the built in dot_product
                     ! or the MKL ddot() function
                     call invdetcalc(deta(1:nalpha,i,nci),nalpha,nalpha,detia(i,1:nalpha,nci),deter(1,nci))
                  endif
                  !deta(1:nalpha,i,nci) = olda(:) ! please take a look this line, It might be unnecessary!
                 else  ! detsRepLst. already calculated
                  deter(1,nci)=deter(1,detsRepLst(nci,1))
                 endif
               enddo CILOOPA
            else
               olddet(1:ndet) = deter(2,:)
               ! CI loop over products of determinants
               CILOOPB: do nci=1,ndet
                 if (detsRepLst(nci,2)==0) then ! calculation only for new dets
                  i = ie - nalpha
                  ii = ie
                  !oldb(1:nbeta) = detb(:,i,nci) ! please take a look this line. It might be unnecessary!
                  ! Construct matrices for determinants
                  detb(1:nbeta, i, nci) = mat(mclist(nalpha + 1:nalpha + nbeta,nci), ii, 1)

                  if (nbeta==1) then
                     deter(2,nci) = detb(1,1,nci)
                  else
                     ! equivalent to
                     ! deter(2,nci) = deter(2,nci) * dot_product(detb(1:nbeta,i,nci),detib(i,1:nbeta,nci))

                     call invdetcalc(detb(1:nbeta,i,nci),nbeta,nbeta,detib(i,1:nbeta,nci),deter(2,nci))

                  endif
                  !detb(1:nbeta,i,nci) = oldb(:) ! please take a look this line, It might be unnecessary!
                 else ! detsRepLst already calculated
                  deter(2,nci)=deter(2,detsRepLst(nci,2))
                 endif !detsRepLst
               enddo CILOOPB
            endif
     end subroutine internal_update_det

     subroutine internal_update_det_and_inv() !updates determinant and inv matrix

           if(ie <= nalpha) then
               olddet(1:ndet) = deter(1,:)
               detiaOne=detia
               !detibOne=detib
               ! CI loop over products of determinants
               CILOOPA: do nci=1,ndet
                if (detsRepLst(nci,1)==0) then  ! calculation only for new dets
                  i = ie
                  ii = ie
                  ! Construct matrices for determinants
                  deta(1:nalpha, i, nci) = mat(mclist(1:nalpha,nci), ii, 1)

                  if (nalpha==1) then
                     deter(1,nci) = deta(1,1,nci)
                  else
                     call invcupd(deta(1:nalpha,i,nci),i,nalpha,nalpha,detiaOne(1:nalpha,1:nalpha,nci),deter(1,nci))

                  endif
                 else  ! detsRepLst. already calculated
                  deter(1,nci)=deter(1,detsRepLst(nci,1))
                  detiaOne(1:nalpha,1:nalpha,nci)=detiaOne(1:nalpha,1:nalpha,detsRepLst(nci,1))
                 endif
               enddo CILOOPA
            else
               olddet(1:ndet) = deter(2,:)
               !detiaOne=detia
               detibOne=detib

               ! CI loop over products of determinants
               CILOOPB: do nci=1,ndet
                 if (detsRepLst(nci,2)==0) then ! calculation only for new dets
                  i = ie - nalpha
                  ii = ie

                  ! Construct matrices for determinants
                  detb(1:nbeta, i, nci) = mat(mclist(nalpha + 1:nalpha + nbeta,nci), ii, 1)

                  if (nbeta==1) then
                     deter(2,nci) = detb(1,1,nci)
                  else
                     call invcupd(detb(1:nbeta,i,nci),i,nbeta,nbeta,detibOne(1:nbeta,1:nbeta,nci),deter(2,nci))
                  endif
                 else ! detsRepLst already calculated
                  deter(2,nci)=deter(2,detsRepLst(nci,2))
                  detibOne(1:nbeta,1:nbeta,nci)=detibOne(1:nbeta,1:nbeta,detsRepLst(nci,2))
                 endif !detsRepLst
               enddo CILOOPB
            endif
     end subroutine internal_update_det_and_inv

   end subroutine mdet2calc

   subroutine resetToOld(ie)
      integer, intent(in) :: ie
      if(ie <= nalpha) then
         deter(1,1:ndet) = olddet(:)
      else
         deter(2,1:ndet) = olddet(:)
      endif
   end subroutine resetToOld


   subroutine lapack_inv(N, A, deti, det, error_code)

      integer, parameter :: nb=64      ! Note:  nb=ilaenv(1,'dgetri','u',N,N,-1,-1)

      integer, intent(in)         :: N                     ! real dimension of matrix A
      real(r8), intent(inout) :: A(:,:)                ! A and deti (A_LU) = P*L*U
      real(r8), intent(inout) :: deti(:,:)             ! LU decomp. and inverse of A
      real(r8), intent(out)   :: det                   ! determinant of A
      integer, intent(out)        :: error_code                 

      integer lwork                 ! work space dimensions
      real(r8)  p                     ! permutation P (see ipiv)
      real(r8)  work(N*nb)            ! workspace!
      integer ipiv(N)               ! pivot indices
      integer info1,info2
      integer i,j

      error_code = MDET_NONE

      !  call dcopy(ndmax*ndmax,A,1,deti,1) for a continuous block
      deti = A

      !  LU decomposition using Lapack/NAG routine
      call dgetrf(N,N,deti,size(deti,1),ipiv(1),info1)

      if (info1 == 0) then

         ! calculate determinant
         det = 1.0_r8

         do i=1,N
            if (ipiv(i) /= i) then
               p = -1.0_r8
            else
               p = 1.0_r8
            end if
            det = det * p * deti(i,i)
         end do

         lwork = size(deti,1)*nb
         call dgetri(N,deti(1,1),size(deti,1),ipiv(1),work(1),lwork,info2)

         if (info2 /= 0) then
            if (logmode >= 2) write(iull,*) " lapack_inv: LAPACK (dgetri) error in inverse INFO =",info2
            error_code = MDET_INV_ERR
         end if

      else

         det = 0.0_r8

         if (info1 > 0) then
            error_code = MDET_LU_ZERO_DET
            if (logmode >= 4) write(iull,*) " lapack_inv: LAPACK (dgetrf) singular matrix in LU decomp. INFO =", info1
         else
            error_code = MDET_LU_ERR
            if (logmode >= 2) write(iull,*) " lapack_inv: LAPACK (dgetrf) error in LU decomp. INFO =", info1
         end if
      
      end if

   end subroutine lapack_inv


   subroutine lapack_det(N, A, col, k, det, error_code)
      
      ! calculate determinant of A with k-th column replaced by col with LU decomposition
      ! calculate determinant of A for k==0
      integer, intent(in)  :: N                     ! real dimension of matrix A
      integer, intent(in)  :: k                     ! replace k-th column
      real(r8), intent(in) :: A(:,:)                ! A 
      real(r8), intent(in) :: col(:)                ! column
      real(r8), intent(inout):: det                   ! determinant of A
      integer, intent(inout) :: error_code 

      real(r8) :: A1(N, N), p
      integer ipiv(N)               ! pivot indices
      integer info, i

      error_code = MDET_NONE

      call assert(size(col)==N .and. size(A,1)==N, "lapack_det: size mismatch")

      A1 = A
      if (k > 0 .and. k <= N) then
         A1(:, k) = col 
      end if

      call dgetrf(N, N, A1, N, ipiv(1), info)

      if (info == 0) then

         ! calculate determinant
         det = 1.0_r8

         do i = 1, N
            if (ipiv(i) /= i) then
               p = -1.0_r8
            else
               p = 1.0_r8
            end if
            det = det * p * A1(i,i)
         end do

      else
         det = 0.0_r8
         if (info > 0) then
            error_code = MDET_LU_ZERO_DET
            if (logmode >= 4) write(iull,*) " lapack_det: LAPACK (dgetrf) singular matrix in LU decomp. INFO =", info
         else
            error_code = MDET_LU_ERR
            if (logmode >= 2) write(iull,*) " lapack_det: LAPACK (dgetrf) error in LU decomp. INFO =", info
         end if
      end if

   end subroutine lapack_det

   function relative_error(value1, value2)
      real(r8), intent(in) :: value1, value2 
      real(r8) :: relative_error
      real(r8), parameter :: min_value = 1.0e-15_r8
      relative_error = abs(value1 - value2) / max(abs(value1), min_value) 
      if (max(abs(value1), abs(value2)) < min_value) relative_error = min_value
   end function relative_error

END MODULE multiDet_m
