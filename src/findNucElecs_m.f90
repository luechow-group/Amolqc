! Copyright (C) 2012-2015, 2017-2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module findNucElecs_m

!  note that there is also Rene's findcore module (findcore_m.f90)
!
!  this module uses simple array data structure to identify
!  and possibly swap electrons close to the nuclei
!
   use kinds_m, only: r8
   use error_m
   use global_m
   use wfData_m
   use sorting_m, only: inverse_permutation

   implicit none

   private
   public :: createNucList, find1sCoreElecs, findNewNucElecs, findNucElecs, findNucElecsSlist, findAndSetNucElecs
   public :: writeElecPositions


contains

   subroutine createNucList(nucList, heavyAtomsOnly)
      ! creates a list of nucleus numbers taken from geometry data (=atoms(a)) but w/o pseudoatoms
      ! and w/o H and He atoms when heavyAtomsOnly=.true.
      integer, allocatable, intent(inout) :: nucList(:)
      logical, intent(in), optional       :: heavyAtomsOnly
      integer a, j, n
      logical heavyOnly

      heavyOnly = .true.
      if (present(heavyAtomsOnly)) heavyOnly = heavyAtomsOnly

      j = 0
      do a = 1, getNNuc()
         if (heavyOnly) then
            if ( atoms(a)%za > 2 .and. (.not.atoms(a)%ecp) ) then
               j = j + 1
            end if
         else
            if (.not.atoms(a)%ecp) then
               j = j + 1
            end if
         end if            
      end do
      n = j
      if (allocated(nucList)) deallocate(nucList)
      allocate(nucList(n))
      j = 0
      do a = 1, getNNuc()
         if (heavyOnly) then
            if ( atoms(a)%za > 2 .and. (.not.atoms(a)%ecp) ) then
               j = j + 1
               nucList(j) = a
            end if
         else
            if (.not.atoms(a)%ecp) then
               j = j + 1
               nucList(j) = a
            end if
         end if            
      end do
   end subroutine createNucList


   subroutine find1sCoreElecs(nucList, x, y, z, idx)
      ! identifies K shell core electrons as closest alpha/beta electrons
      ! swap electrons such that core electrons are the first in alpha and
      ! beta list and put core electrons at nucleus position
      integer, intent(in)   :: nucList(:)         ! list of (heavy atom) nuclei w/o pseudo potential
      real(r8), intent(inout) :: x(:), y(:), z(:)   ! coords
      integer, intent(inout), optional :: idx(:)  ! j=idx(i) means x(i) was before x(j). idx is updated
                                                  ! meaning x(i) was x(j) prior to all swaps
      integer i, ia, ib, a, iamin, ibmin, ica, icb, na, itmp
      real(r8) dmin
      real(r8) rai(getNNuc(), getNElec())
      real(r8), parameter :: EPS = 0.d0   ! put electrons EPS outside nucleus to avoid singularities in Vpot!

      call calcNucElecDists(x, y, z, rai)

      na = getNAlpha()
      do ia = 1, size(nucList)
         a = nucList(ia)
         iamin = minloc(rai(a,ia:na), dim=1) + ia - 1
         ! copy electron ia to iamin (deleting iamin) and set Nuc electron ia to nucleus pos (slightly outside by EPS)
         ! copy also rai entry
         ica = ia
         x(iamin) = x(ica); y(iamin) = y(ica); z(iamin) = z(ica)
         rai(:,iamin) = rai(:,ica)
         x(ica) = atoms(a)%cx; y(ica) = atoms(a)%cy; z(ica) = atoms(a)%cz+EPS
         if (present(idx)) then
            itmp = idx(ica); idx(ica) = idx(iamin); idx(iamin) = itmp
         end if
      end do
      do ib = 1, size(nucList)
         a = nucList(ib)
         ibmin = minloc(rai(a,na+ib:getNElec()), dim=1) + na + ib - 1
         icb = na + ib
         x(ibmin) = x(icb); y(ibmin) = y(icb); z(ibmin) = z(icb)
         rai(:,ibmin) = rai(:,icb)
         x(icb) = atoms(a)%cx; y(icb) = atoms(a)%cy; z(icb) = atoms(a)%cz-EPS
         if (present(idx)) then
            itmp = idx(icb); idx(icb) = idx(ibmin); idx(ibmin) = itmp
         end if
      end do
   end subroutine find1sCoreElecs

   subroutine findNewNucElecs(nucList, thresh, naNucElecs, aNucElecList, nbNucElecs, bNucElecList, &
                              x, y, z, found, idx)
      ! identifies new electrons with distance to a nucleus < thresh/Z (nuclei in nucList only)
      ! assumes sorted electrons such that the first n[a|b]NucElecs are on entry at nuclei
      ! [a|b]NucElecList contains the corresponding nuclei
      ! swaps these electrons such that they are at the end of list of nuclear electrons
      ! electrons identified as electrons at nucleus are put on nucleus position
      ! idx keeps track of the electron permutations
      integer, intent(in)     :: nucList(:)      ! list of possible nuclei
      real(r8), intent(in)      :: thresh          ! thresh / Z is threshold to identify electron at nucleus
      integer, intent(inout)  :: naNucElecs      ! # of alpha elecs at nucleus (inout!)
      integer, intent(inout)  :: aNucElecList(:) ! list of nuclei of the alpha elecs at nucleus
      integer, intent(inout)  :: nbNucElecs      ! # of beta elecs at nucleus (inout!)
      integer, intent(inout)  :: bNucElecList(:) ! list of nuclei of the beta elecs at nucleus
      real(r8), intent(inout)   :: x(:),y(:),z(:)  ! electron coords
      logical, intent(out)    :: found
      integer, intent(inout), optional :: idx(:) ! j=idx(i) means x(i) was originally x(j)

      integer i, ia, ib, a, iamin, ibmin, ica, icb, na, itmp, k
      real(r8) dmin
      real(r8) rai(getNNuc(), getNElec())
      real(r8), parameter :: EPS = 0.d0       ! put electrons EPS outside nucleus to avoid singularities!

      found = .false.
      call calcNucElecDists(x, y, z, rai)   ! not all distances are necessary!

      ! alpha electrons
      na = getNAlpha()
      do k = 1, size(nucList)
         ia = naNucElecs + 1
         a = nucList(k)
         if ( any(aNucElecList(1 : naNucElecs) == a) ) cycle    ! never two alpha elecs at one nuc
         dmin = minval(rai(a, ia : na))
         iamin = minloc(rai(a, ia : na), dim = 1) + ia - 1
         if (dmin < thresh / atoms(a)%za) then
            ! copy electron ia to iamin (deleting iamin) and set NucElec electron ia to nucleus pos (slightly outside by EPS)
            ! copy also rai entry
            ica = ia
            x(iamin) = x(ica); y(iamin) = y(ica); z(iamin) = z(ica)
            rai(:, iamin) = rai(:, ica)
            x(ica) = atoms(a)%cx; y(ica) = atoms(a)%cy; z(ica) = atoms(a)%cz + EPS
            naNucElecs = naNucElecs + 1
            aNucElecList(naNucElecs) = a
            found = .true.
            if (present(idx)) then
               itmp = idx(ica); idx(ica) = idx(iamin); idx(iamin) = itmp
            end if
        end if
      end do

      ! beta electrons
      do k = 1, size(nucList)
         ib = nbNucElecs + 1
         a = nucList(k)
         if (any(bNucElecList(1 : nbNucElecs) == a)) cycle    ! never two beta elecs at one nuc
         dmin = minval(rai(a, na + ib : getNElec()))
         ibmin = minloc(rai(a, na + ib : getNElec()), dim = 1) + na + ib - 1
         if (dmin < thresh / atoms(a)%za) then
            icb = na + ib
            x(ibmin) = x(icb); y(ibmin) = y(icb); z(ibmin) = z(icb)
            rai(:, ibmin) = rai(:, icb)
            x(icb) = atoms(a)%cx; y(icb) = atoms(a)%cy; z(icb) = atoms(a)%cz-EPS
            nbNucElecs = nbNucElecs + 1
            bNucElecList(nbNucElecs) = a
            found = .true.
            if (present(idx)) then
               itmp = idx(icb); idx(icb) = idx(ibmin); idx(ibmin) = itmp
            end if
         end if
      enddo

   end subroutine findNewNucElecs


   subroutine findNucElecs(thresh, x, y, z, naNuc, nbNuc, aNucElec, bNucElec, nucList)
      ! identifies elecs close to a nucleus core (if dist < thresh / Z)
      ! build a/bNucArray(na/bNuc)
      ! a/bNucElec differ from a/bNucElecList (findNewNucElecs). Here it is an array
      ! containing for each nucleus the electron number if there is an elec with dist < thresh / Z 
      ! do NOT swap these electrons (use findNewNucElecs for this)
      ! do NOT check if electrons are sorted
      real(r8), intent(in)             :: thresh            ! threshold to identify H core electrons
      real(r8), intent(in)             :: x(:),y(:),z(:)    ! coords
      integer, intent(out)           :: naNuc, nbNuc      ! alpha, beta elecs at nucleus (not varied)
      integer, intent(out), optional :: aNucElec(:), bNucElec(:) ! nuc a -> elec i (at nuc)
      integer, intent(in), optional  :: nucList(:)        ! check only nuclei in list (e.g. no ecp nuclei)
      integer i, ia, ib, a, iamin, ibmin, ica, icb, na, k
      real(r8) dmin, tmp
      real(r8) rai(getNNuc(), getNElec())

      call assert(size(x)==getNElec(), '(findNucElecs): illegal sizes on entry')

      call calcNucElecDists(x, y, z, rai)
      na = getNAlpha()

      if (present(nucList)) then

         ! alpha electrons
         naNuc = 0
         if (present(aNucElec)) aNucElec = 0
         do k = 1, size(nucList)
            a = nucList(k)
            dmin = minval(rai(a, 1 : na))
            iamin = minloc(rai(a, 1 : na), dim = 1)
            if (dmin < thresh / atoms(a)%za) then
               if (present(aNucElec)) aNucElec(a) = iamin
               naNuc = naNuc + 1
            end if
         end do

         ! beta electrons
         nbNuc = 0
         if (present(bNucElec)) bNucElec = 0
         do k = 1, size(nucList)
            a = nucList(k)
            dmin = minval(rai(a, na + 1 : getNElec()))
            ibmin = minloc(rai(a, na + 1 : getNElec()), dim = 1) + na
            if (dmin < thresh / atoms(a)%za) then
               if (present(bNucElec)) bNucElec(a) = ibmin
               nbNuc = nbNuc + 1
           end if
           !!print'(i3,2f10.5,l4,10i3)',a,dmin,thresh,swap,ibmin,nbNuc,bNucElec
         end do

      else

         ! alpha electrons
         naNuc = 0
         if (present(aNucElec)) aNucElec = 0
         do a = 1, getNNuc()
            dmin = minval(rai(a, 1 : na))
            iamin = minloc(rai(a, 1 : na), dim = 1)
            if (dmin < thresh / atoms(a)%za) then
               if (present(aNucElec)) aNucElec(a) = iamin
               naNuc = naNuc + 1
            end if
         end do

         ! beta electrons
         nbNuc = 0
         if (present(bNucElec)) bNucElec = 0
         do a = 1, getNNuc()
            dmin = minval(rai(a, na + 1 : getNElec()))
            ibmin = minloc(rai(a, na + 1 : getNElec()), dim = 1) + na
            if (dmin < thresh / atoms(a)%za) then
               if (present(bNucElec)) bNucElec(a) = ibmin
               nbNuc = nbNuc + 1
           end if
           !!print'(i3,2f10.5,l4,10i3)',a,dmin,thresh,swap,ibmin,nbNuc,bNucElec
         end do
      end if
   end subroutine findNucElecs

   subroutine findNucElecsSlist(thresh, x, y, z, slist, nnuc)
      ! identifies elecs close to a nucleus core (if dist < thresh / Z)
      ! build a/bNucArray(na/bNuc)
      ! a/bNucElec differ from a/bNucElecList (findNewNucElecs). Here it is an array
      ! containing for each nucleus the electron number if there is an elec with dist < thresh / Z 
      ! do NOT swap these electrons (use findNewNucElecs for this)
      ! do NOT check if electrons are sorted
      real(r8), intent(in)             :: thresh            ! threshold to identify H core electrons
      real(r8), intent(in)             :: x(:),y(:),z(:)    ! coords
      integer, intent(inout)         :: slist(:)          ! slist(i)==a: electron i at nuc a; ==0: not at a nuc
      integer, intent(inout)         :: nnuc              ! 
      integer i, k, n, iamin
      real(r8) dmin
      real(r8) rai(getNNuc(), getNElec())

      call assert(size(x)==getNElec(), '(findNucElecsSlist): illegal xsize on entry')
      call assert(size(slist)==size(x), "(findNucElecsSlist): illegal slist size on entry")

      call calcNucElecDists(x, y, z, rai)
      n = getNElec()
      slist = 0
      nnuc = 0

      do i = 1, n
         dmin = minval(rai(1 : getNNuc(), i))
         iamin = minloc(rai(1 : getNNuc(), i), dim = 1)
         if (dmin < thresh) then
            nnuc = nnuc + 1
            slist(i) = iamin
         end if
      end do
   end subroutine findNucElecsSlist


   subroutine findAndSetNucElecs(thresh, x, y, z, naNuc, nbNuc, isSet, aNucElec, bNucElec, nucList)
      ! identifies elecs close to a nucleus (if dist < thresh / Z)
      ! sets them at nuc position
      ! build a/bNucArray(na/bNuc)
      ! a/bNucElec differ from a/bNucElecList (findNewNucElecs). Here it is an array
      ! containing for each nucleus the electron number if there is an elec with dist < thresh / Z 
      ! do NOT swap these electrons (use findNewNucElecs for this)
      ! do NOT check if electrons are sorted
      real(r8), intent(in)             :: thresh            ! threshold to identify electrons at nucleus
      real(r8), intent(inout)          :: x(:),y(:),z(:)    ! coords
      integer, intent(out)           :: naNuc, nbNuc      ! alpha, beta elecs at nucleus
      logical, intent(out)           :: isSet             ! .true. if coords are changed on exit
      integer, intent(inout), optional :: aNucElec(:), bNucElec(:) ! nuc a -> elec i (at nuc)
      integer, intent(in), optional  :: nucList(:)        ! check only nuclei in list (e.g. no ecp nuclei)
      integer i, ia, ib, a, iamin, ibmin, ica, icb, na, k
      real(r8) dmin, tmp
      real(r8) rai(getNNuc(), getNElec())

      call assert(size(x)==getNElec() .and. size(aNucElec)==getNNuc(), &
                 '(findNucElecs): illegal sizes on entry')

      call calcNucElecDists(x, y, z, rai)
      na = getNAlpha()
      isSet = .false.

      if (present(nucList)) then

         ! alpha electrons
         naNuc = 0
         do k = 1, size(nucList)
            a = nucList(k)
            dmin = minval(rai(a, 1 : na))
            iamin = minloc(rai(a, 1 : na), dim = 1)
            if (dmin < thresh / atoms(a)%za) then
               if (present(aNucElec)) then
                  if (aNucElec(a) == 0) then
                     aNucElec(a) = iamin                     
                     x(iamin) = atoms(a)%cx
                     y(iamin) = atoms(a)%cy
                     z(iamin) = atoms(a)%cz
                     isSet = .true.
                  else if (aNucElec(a) /= iamin) then
                     call abortp("FINDANDSETNUCELECS: this should not happen")
                  end if
               else
                  x(iamin) = atoms(a)%cx
                  y(iamin) = atoms(a)%cy
                  z(iamin) = atoms(a)%cz
                  isSet = .true.
               end if
               naNuc = naNuc + 1
            end if
         end do

         ! beta electrons
         nbNuc = 0
         do k = 1, size(nucList)
            a = nucList(k)
            dmin = minval(rai(a, na + 1 : getNElec()))
            ibmin = minloc(rai(a, na + 1 : getNElec()), dim = 1) + na
            if (dmin < thresh / atoms(a)%za) then
               if (present(bNucElec)) then
                  if (bNucElec(a) == 0) then
                     bNucElec(a) = ibmin                     
                     x(ibmin) = atoms(a)%cx
                     y(ibmin) = atoms(a)%cy
                     z(ibmin) = atoms(a)%cz
                     isSet = .true.
                  else if (bNucElec(a) /= ibmin) then
                     call abortp("FINDANDSETNUCELECS: this should not happen")
                  end if
               else
                  x(ibmin) = atoms(a)%cx
                  y(ibmin) = atoms(a)%cy
                  z(ibmin) = atoms(a)%cz
                  isSet = .true.
               end if
               nbNuc = nbNuc + 1
            end if
           !!print'(i3,2f10.5,l4,10i3)',a,dmin,thresh,swap,ibmin,nbNuc,bNucElec
         end do

      else

         ! alpha electrons
         naNuc = 0
         do a = 1, getNNuc()
            dmin = minval(rai(a, 1 : na))
            iamin = minloc(rai(a, 1 : na), dim = 1)
            if (dmin < thresh / atoms(a)%za) then
               if (present(aNucElec)) then
                  if (aNucElec(a) == 0) then
                     aNucElec(a) = iamin                     
                     x(iamin) = atoms(a)%cx
                     y(iamin) = atoms(a)%cy
                     z(iamin) = atoms(a)%cz
                     isSet = .true.
                  else if (aNucElec(a) /= iamin) then
                     call abortp("FINDANDSETNUCELECS: this should not happen")
                  end if
               else
                  x(iamin) = atoms(a)%cx
                  y(iamin) = atoms(a)%cy
                  z(iamin) = atoms(a)%cz
                  isSet = .true.
               end if
               naNuc = naNuc + 1
            end if
         end do

         ! beta electrons
         nbNuc = 0
         do a = 1, getNNuc()
            dmin = minval(rai(a, na + 1 : getNElec()))
            ibmin = minloc(rai(a, na + 1 : getNElec()), dim = 1) + na
            if (dmin < thresh / atoms(a)%za) then
               if (present(bNucElec)) then
                  if (bNucElec(a) == 0) then
                     bNucElec(a) = ibmin                     
                     x(ibmin) = atoms(a)%cx
                     y(ibmin) = atoms(a)%cy
                     z(ibmin) = atoms(a)%cz
                     isSet = .true.
                  else if (bNucElec(a) /= ibmin) then
                     call abortp("FINDANDSETNUCELECS: this should not happen")
                  end if
               else
                  x(ibmin) = atoms(a)%cx
                  y(ibmin) = atoms(a)%cy
                  z(ibmin) = atoms(a)%cz
                  isSet = .true.
               end if
               nbNuc = nbNuc + 1
            end if
           !!print'(i3,2f10.5,l4,10i3)',a,dmin,thresh,swap,ibmin,nbNuc,bNucElec
         end do
      end if
   end subroutine findAndSetNucElecs


   subroutine writeElecPositions(iunit, x, y, z, m, n, ff, tag, idx)
      ! write electron positions in "ref" format to unit "iunit"
      integer, intent(in) :: iunit                ! file unit to write to
      real(r8), intent(in)  :: x(:), y(:), z(:)     ! electron coordinates
      integer, intent(in) :: m, n                 ! counter for 1st line
      real(r8), intent(in)  :: ff                   ! function value
      character(len=*)    :: tag                  ! 1st line tag
      integer, optional, intent(in) :: idx(:)     ! permutation, print original order
      integer, allocatable :: p(:), invp(:)
      integer i

      write(iunit,'(a, 2i6, a, f14.5, a)') trim(tag)//":", m, n, "   F(Max):",  ff, "  found: 0 0 0 0 0"
      write(iunit,'(i5)') getNElec()
      if (present(idx)) then
         allocate(p(size(idx)), invp(size(idx)))
         p = idx
         invp = [ (i, i=1, size(idx)) ]
         call inverse_permutation(p, invp)
         do i = 1, getNElec()
            write(iunit,'(3f14.7)') x(invp(i)), y(invp(i)), z(invp(i))
         enddo
      else
         do i = 1, getNElec()
            write(iunit,'(3f14.7)') x(i), y(i), z(i)
         enddo
      end if
   end subroutine writeElecPositions


end module findNucElecs_m
