! Copyright (C) ca. 2000, 2012-2014, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module atom_m
   use kinds_m, only: r8
   use error_m
   use sorting_m, only: quickSortIndex
   implicit none

   type atom
      real(r8)          :: cx = 0   ! x,y,z coordinates (in bohr)
      real(r8)          :: cy = 0
      real(r8)          :: cz = 0
      real(r8)          :: za = 0   ! core charge. With ECPs, this may be unequal to elemIdx
      integer           :: elemIdx = 0
      integer           :: pa = 0   ! charge/separated electrons for sampling
      integer           :: sa = 0   ! index for atoms that are "the same"
      character(len=2)  :: elem = '  '
      character(len=20) :: ba = ''
      logical           :: ecp=.false.
   contains
      procedure :: Get_position => atom_getPosition
      procedure :: Get_atomic_number => atom_getAtomicNumber
   end type

   integer, parameter          :: PSEMAX=118
   character(len=2), parameter :: pse(0:PSEMAX) =    &
     ['X ','H ','He',&
      'Li','Be','B ','C ','N ','O ','F ','Ne',&
      'Na','Mg','Al','Si','P ','S ','Cl','Ar',&
      'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe',&
      'Co','Ni','Cu','Zn','Ga','Ge','As','Se',&
      'Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo',&
      'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',&
      'Sb','Te','I ','Xe','Cs','Ba','La','Ce',&
      'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',&
      'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ',&
      'Re','Os','Ir','Pt','Au','Hg','Tl','Pb',&
      'Bi','Po','At','Rn','Fr','Ra','Ac','Th',&
      'Pa','U ','Np','Pu','Am','Cm','Bk','Cf',&
      'Es','Fm','Md','No','Lr','Rf','Db','Sg',&
      'Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl',&
      'Mc','Lv','Ts','Og']
   real(r8),parameter  :: Bohr2angs = 0.529177249d0


   ! in and output
   logical :: mWithSA = .false.
   logical :: mWithCharge = .false.
   logical :: mWithBasis = .false.
   logical :: mAngstrom = .true.
   ! same atom properties
   logical :: mIgnoreHAtoms = .false.

   ! atomic distance array
   integer :: mNAtoms = 0
   real(r8), allocatable :: mAtomDist(:,:)

   ! thresholds for atoms_whatPosition
   real(r8)  :: mNucThresh = 0.01
   real(r8)  :: mCoreThresh = 2.0
   real(r8)  :: mBondDiff = 0.4


   private
   public :: atom,atoms_configure,atom_init,atoms_read,atoms_readFromFile,atoms_write,atoms_withSA, &
             atoms_withBasis,atoms_withCharge,atom_getSA,atoms_ignoreHAtoms,atoms_doNotIgnoreHAtoms, &
             atoms_getNSCenter,atoms_getNCLast,atoms_countElectrons,atoms_checkAtomOrder,atoms_getBox, &
             getPSEIdx,getPSESymbol,atoms_initPositionThresholds,atoms_whatPosition,atoms_whatType, &
             atoms_spaceWarpTransform, atom_getPosition, atom_getDistance, atoms_getPositionMatrix, &
             atom_getAtomicNumber

contains



   subroutine atoms_configure(angstrom,withSA,withCharge,withBasis)
   !-------------------------------------------------------------!
   logical, intent(in)           :: angstrom
   logical, intent(in)           :: withSA
   logical, intent(in)           :: withCharge
   logical, intent(in)           :: withBasis
   mAngstrom = angstrom
   mWithSA = withSA
   mWithCharge = withCharge
   mWithBasis = withBasis
   end subroutine atoms_configure

   elemental function atom_getAtomicNumber(self) result(atomic_number)
      class(atom),intent(in) :: self
      integer :: atomic_number

      atomic_number = self%elemIdx
   end function atom_getAtomicNumber

   function atom_getPosition(self) result(position)
      class(atom),intent(in) :: self
      real(r8) :: position(3)

      position(1) = self%cx
      position(2) = self%cy
      position(3) = self%cz
   end function atom_getPosition

   function atom_getDistance(self, other) result(distance)
      type(atom),intent(in) :: self, other
      real(r8) :: distance

      distance = NORM2(atom_getPosition(self) - atom_getPosition(other))
   end function atom_getDistance

   subroutine atom_init(self,x,y,z,atomIdx,symbol,charge,satom,basis)
   !----------------------------------------------------------------!
   type(atom),intent(inout)      :: self
   real(r8), intent(in)            :: x,y,z
   integer, intent(in), optional :: atomIdx
   integer, intent(in), optional :: charge
   integer, intent(in), optional :: satom
   character(len=2), optional    :: symbol
   character(len=20), optional   :: basis

   call assert(present(atomIdx).or.present(symbol),'atom: either element number or symbol required')
   call assert((present(charge).eqv.mWithCharge) .and. (present(satom).eqv.mWithSA) .and. (present(basis).eqv.mWithBasis), &
               'atom: init inconsistent with configuration')

   if (mAngstrom) then
      self%cx = x / bohr2angs
      self%cy = y / bohr2angs
      self%cz = z / bohr2angs
   else
      self%cx = x
      self%cy = y
      self%cz = z
   end if

   if (present(atomIdx)) then
      self%elemIdx = atomIdx
      self%elem = getPSESymbol(atomIdx)
   else if (present(symbol)) then
      self%elemIdx = getPSEIdx(symbol)
      self%elem = symbol
   else
      call error('atom: either element number or symbol required')
   end if
   self%za = self%elemIdx

   if (present(charge)) then
      self%pa = charge
   end if

   if (present(basis)) then
      self%ba = basis
   end if

   if (present(satom)) then
      self%sa = satom
   end if
   end subroutine atom_init

   function atoms_getPositionMatrix(atoms) result(matrix)
      type(atom),intent(in) :: atoms(:)
      real(r8) :: matrix(3, SIZE(atoms))
      integer :: i

      do i = 1, SIZE(atoms)
         matrix(:, i) = atom_getPosition(atoms(i))
      end do
   end function atoms_getPositionMatrix

   subroutine atoms_readFromFile(self,iu,withSymbol)
   !-----------------------------------------------!
   type(atom),intent(inout)      :: self(:)
   integer, intent(in)           :: iu
   logical, intent(in)           :: withSymbol
   real(r8)   :: x,y,z
   integer  :: aidx,ch,sa,i
   character(len=2)  :: sym
   character(len=20) :: bas
   logical ordered

   if (withSymbol) then
      if (.not.mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(iu,*) sym,x,y,z
            call atom_init(self(i),x,y,z,symbol=sym)
         end do
      else if (.not.mWithSA .and. .not.mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(iu,*) sym,x,y,z,bas
            call atom_init(self(i),x,y,z,symbol=sym,basis=bas)
         end do
      else if (.not.mWithSA .and. mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(iu,*) sym,x,y,z,ch
            call atom_init(self(i),x,y,z,symbol=sym,charge=ch)
         end do
      else if (.not.mWithSA .and. mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(iu,*) sym,x,y,z,ch,bas
            call atom_init(self(i),x,y,z,symbol=sym,charge=ch,basis=bas)
         end do
      else if (mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(iu,*) sym,sa,x,y,z
            call atom_init(self(i),x,y,z,symbol=sym,satom=sa)
         end do
      else if (mWithSA .and. .not.mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(iu,*) sym,sa,x,y,z,bas
            call atom_init(self(i),x,y,z,symbol=sym,satom=sa,basis=bas)
         end do
      else if (mWithSA .and. mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(iu,*) sym,sa,x,y,z,ch
            call atom_init(self(i),x,y,z,symbol=sym,charge=ch,satom=sa)
         end do
      else if (mWithSA .and. mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(iu,*) sym,sa,x,y,z,ch,bas
            call atom_init(self(i),x,y,z,symbol=sym,charge=ch,satom=sa,basis=bas)
         end do
      end if
   else ! read atom number instead of element symbol
      if (.not.mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(iu,*) aidx,x,y,z
            call atom_init(self(i),x,y,z,atomIdx=aidx)
         end do
      else if (.not.mWithSA .and. .not.mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(iu,*) aidx,x,y,z,bas
            call atom_init(self(i),x,y,z,atomIdx=aidx,basis=bas)
         end do
      else if (.not.mWithSA .and. mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(iu,*) aidx,x,y,z,ch
            call atom_init(self(i),x,y,z,atomIdx=aidx,charge=ch)
         end do
      else if (.not.mWithSA .and. mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(iu,*) aidx,x,y,z,ch,bas
            call atom_init(self(i),x,y,z,atomIdx=aidx,charge=ch,basis=bas)
         end do
      else if (mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(iu,*) aidx,sa,x,y,z
            call atom_init(self(i),x,y,z,atomIdx=aidx,satom=sa)
         end do
      else if (mWithSA .and. .not.mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(iu,*) aidx,sa,x,y,z,bas
            call atom_init(self(i),x,y,z,atomIdx=aidx,satom=sa,basis=bas)
         end do
      else if (mWithSA .and. mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(iu,*) aidx,sa,x,y,z,ch
            call atom_init(self(i),x,y,z,atomIdx=aidx,charge=ch,satom=sa)
         end do
      else if (mWithSA .and. mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(iu,*) aidx,sa,x,y,z,ch,bas
            call atom_init(self(i),x,y,z,atomIdx=aidx,charge=ch,satom=sa,basis=bas)
         end do
      end if
   end if

   if (.not.mWithSA) then
      call calcSameAtomIndices(self)
   end if

   ordered = atoms_checkAtomOrder(self)
   call assert(ordered,'atom:atoms_read: atoms must be ordered')

   call atoms_initDistances(self)

   end subroutine atoms_readFromFile



   subroutine atoms_read(self,lines,nl,withSymbol)
   !---------------------------------------------!
   type(atom),intent(inout)      :: self(:)
   character(len=*), intent(in)  :: lines(:)       ! string array containing geometry
   integer, intent(in)           :: nl             ! actual # of lines in string array
   logical, intent(in)           :: withSymbol
   real(r8)   :: x,y,z
   integer  :: aidx,ch,sa,i,of
   character(len=2)  :: sym
   character(len=20) :: bas
   logical ordered

   call assert(nl>=size(self)+2,'(atoms_read): illegal format in atoms input')
   of = 2  ! start reading atoms from line 3
   if (withSymbol) then
      if (.not.mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) sym,x,y,z
            call atom_init(self(i),x,y,z,symbol=sym)
         end do
      else if (.not.mWithSA .and. .not.mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) sym,x,y,z,bas
            call atom_init(self(i),x,y,z,symbol=sym,basis=bas)
         end do
      else if (.not.mWithSA .and. mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) sym,x,y,z,ch
            call atom_init(self(i),x,y,z,symbol=sym,charge=ch)
         end do
      else if (.not.mWithSA .and. mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) sym,x,y,z,ch,bas
            call atom_init(self(i),x,y,z,symbol=sym,charge=ch,basis=bas)
         end do
      else if (mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) sym,sa,x,y,z
            call atom_init(self(i),x,y,z,symbol=sym,satom=sa)
         end do
      else if (mWithSA .and. .not.mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) sym,sa,x,y,z,bas
            call atom_init(self(i),x,y,z,symbol=sym,satom=sa,basis=bas)
         end do
      else if (mWithSA .and. mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) sym,sa,x,y,z,ch
            call atom_init(self(i),x,y,z,symbol=sym,charge=ch,satom=sa)
         end do
      else if (mWithSA .and. mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) sym,sa,x,y,z,ch,bas
            call atom_init(self(i),x,y,z,symbol=sym,charge=ch,satom=sa,basis=bas)
         end do
      end if
   else ! read atom number instead of element symbol
      if (.not.mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) aidx,x,y,z
            call atom_init(self(i),x,y,z,atomIdx=aidx)
         end do
      else if (.not.mWithSA .and. .not.mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) aidx,x,y,z,bas
            call atom_init(self(i),x,y,z,atomIdx=aidx,basis=bas)
         end do
      else if (.not.mWithSA .and. mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) aidx,x,y,z,ch
            call atom_init(self(i),x,y,z,atomIdx=aidx,charge=ch)
         end do
      else if (.not.mWithSA .and. mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) aidx,x,y,z,ch,bas
            call atom_init(self(i),x,y,z,atomIdx=aidx,charge=ch,basis=bas)
         end do
      else if (mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) aidx,sa,x,y,z
            call atom_init(self(i),x,y,z,atomIdx=aidx,satom=sa)
         end do
      else if (mWithSA .and. .not.mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) aidx,sa,x,y,z,bas
            call atom_init(self(i),x,y,z,atomIdx=aidx,satom=sa,basis=bas)
         end do
      else if (mWithSA .and. mWithCharge .and. .not.mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) aidx,sa,x,y,z,ch
            call atom_init(self(i),x,y,z,atomIdx=aidx,charge=ch,satom=sa)
         end do
      else if (mWithSA .and. mWithCharge .and. mWithBasis) then
         do i=1,size(self)
            read(lines(of+i),*) aidx,sa,x,y,z,ch,bas
            call atom_init(self(i),x,y,z,atomIdx=aidx,charge=ch,satom=sa,basis=bas)
         end do
      end if
   end if

   if (.not.mWithSA) then
      call calcSameAtomIndices(self)
   end if

   ordered = atoms_checkAtomOrder(self)
   call assert(ordered,'atom:atoms_read: atoms must be ordered')

   call atoms_initDistances(self)

   end subroutine atoms_read



   subroutine atoms_write(self,iu)
   !-----------------------------!
   type(atom),intent(inout)      :: self(:)
   integer, intent(in)           :: iu
   integer :: i
   real(r8)  :: x(size(self)),y(size(self)),z(size(self))

   x = self(:)%cx * bohr2angs
   y = self(:)%cy * bohr2angs
   z = self(:)%cz * bohr2angs

   if (.not.mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
      do i=1,size(self)
         write(iu,'(a2,3f10.5)') self(i)%elem,x(i),y(i),z(i)
      end do
   else if (.not.mWithSA .and. .not.mWithCharge .and. mWithBasis) then
      do i=1,size(self)
         write(iu,'(a2,3f10.5,1x,a20)') self(i)%elem,x(i),y(i),z(i),self(i)%ba
      end do
   else if (.not.mWithSA .and. mWithCharge .and. .not.mWithBasis) then
      do i=1,size(self)
         write(iu,'(a2,3f10.5,i3)') self(i)%elem,x(i),y(i),z(i),self(i)%pa
      end do
   else if (.not.mWithSA .and. mWithCharge .and. mWithBasis) then
      do i=1,size(self)
         write(iu,'(a2,3f10.5,i3,1x,a20)') self(i)%elem,x(i),y(i),z(i),self(i)%pa,self(i)%ba
      end do
   else if (mWithSA .and. .not.mWithCharge .and. .not.mWithBasis) then
      do i=1,size(self)
         write(iu,'(a2,i5,3f10.5)') self(i)%elem,self(i)%sa,x(i),y(i),z(i)
      end do
   else if (mWithSA .and. .not.mWithCharge .and. mWithBasis) then
      do i=1,size(self)
         write(iu,'(a2,i5,3f10.5,1x,a20)') self(i)%elem,self(i)%sa,x(i),y(i),z(i),self(i)%ba
      end do
   else if (mWithSA .and. mWithCharge .and. .not.mWithBasis) then
      do i=1,size(self)
         write(iu,'(a2,i5,3f10.5,i3)') self(i)%elem,self(i)%sa,x(i),y(i),z(i),self(i)%pa
      end do
   else if (mWithSA .and. mWithCharge .and. mWithBasis) then
      do i=1,size(self)
         write(iu,'(a2,i5,3f10.5,i3,1x,a20)') self(i)%elem,self(i)%sa,x(i),y(i),z(i),self(i)%pa,self(i)%ba
      end do
   end if
   end subroutine atoms_write

   integer function atoms_countElectrons(self)
      type(atom), intent(in) :: self(:)
      atoms_countElectrons = sum(self(:)%elemIdx)
   end function atoms_countElectrons


   logical function atoms_withSA()
      atoms_withSA = mWithSA
   end function atoms_withSA

   logical function atoms_withCharge()
      atoms_withCharge = mWithCharge
   end function atoms_withCharge

   logical function atoms_withBasis()
      atoms_withBasis = mWithBasis
   end function atoms_withBasis

   integer function atom_getSA(self)
      type(atom), intent(in) :: self
      atom_getSA = self%sa
   end function atom_getSA

   subroutine atoms_ignoreHAtoms()
      mIgnoreHAtoms = .true.
   end subroutine atoms_ignoreHAtoms

   subroutine atoms_doNotIgnoreHAtoms()
      mIgnoreHAtoms = .false.
   end subroutine atoms_doNotIgnoreHAtoms


   integer function atoms_getNSCenter(self)
      type(atom), intent(inout) :: self(:)
      logical ordered

      call assert(size(self)>0,'atom: getNSCenter: empty atom array')
      if (any(self(:)%sa==0)) then
         call calcSameAtomIndices(self)
      end if
      ordered = atoms_checkAtomOrder(self)
      call assert(ordered,'atom: getNSCenter: requires ordered atom input')
      atoms_getNSCenter = self(size(self))%sa
   end function atoms_getNSCenter


   integer function atoms_getNCLast(self)
   !------------------------------------!
      type(atom), intent(in) :: self(:)
      logical ordered
      integer n

      if (mIgnoreHAtoms) then
         ordered = atoms_checkAtomOrder(self)
         call assert(ordered,'atoms_getNCLast: requires ordered atoms')
         call assert(self(size(self))%elemIdx == 1,  &
                     'atoms_getNCLast: ignoring H atoms only when H atoms are last')
         do n=size(self),1,-1
            if (self(n)%elemIdx /= 1) exit
         end do
         atoms_getNCLast = n
      else
         atoms_getNCLast = size(self)
      end if
   end function atoms_getNCLast



   logical function atoms_checkAtomOrder(self)
   !-----------------------------------------!
      type(atom), intent(in) :: self(:)
      integer i,idx
      logical failed

      call assert(size(self)>0,'atom:atoms_checkAtomOrder empty atoms array')
      failed = .false.
      if (self(1)%sa /= 1) then
         failed = .true.
      else
         idx = 1
         do i=2,size(self)
            if (self(i)%sa /= idx) then
               if (self(i)%sa == idx+1) then
                  idx = idx+1
               else
                  failed = .true.
                  exit
               end if
            end if
         end do
      end if
      atoms_checkAtomOrder = .not. failed
   end function atoms_checkAtomOrder


   subroutine atoms_getBox(self,dist,ax,bx,ay,by,az,bz)
   !--------------------------------------------------!
      type(atom), intent(in) :: self(:)
      real(r8), intent(in)     :: dist
      real(r8), intent(inout)  :: ax,bx,ay,by,az,bz
      ! atoms_getBox calculates a tight box around the atoms
      ! with a distance of at least 'dist' along all three axis
      ! on return: intervals [ax,bx], [ay,by], [az,bz] denoting the box.
      ax = minval(self(:)%cx) - dist
      bx = maxval(self(:)%cx) + dist
      ay = minval(self(:)%cy) - dist
      by = maxval(self(:)%cy) + dist
      az = minval(self(:)%cz) - dist
      bz = maxval(self(:)%cz) + dist
   end subroutine atoms_getBox



   subroutine atoms_initDistances(self)
   !----------------------------------!
      type(atom), intent(in) :: self(:)
      integer                :: i,j
      mNAtoms = size(self)
      if (allocated(mAtomDist)) deallocate(mAtomDist)
      allocate(mAtomDist(mNAtoms,mNAtoms))
      mAtomDist = 0.d0
      do i=1,mNAtoms
         do j=i+1,mNAtoms
            mAtomDist(i,j) = sqrt( (self(i)%cx-self(j)%cx)**2 + (self(i)%cy-self(j)%cy)**2 + &
                                    (self(i)%cz-self(j)%cz)**2 )
            mAtomDist(j,i) = mAtomDist(i,j)
         end do
      end do
   end subroutine atoms_initDistances



   subroutine  atoms_initPositionThresholds(nucThresh,coreThresh,bondDiff)
   !---------------------------------------------------------------------!
      real(r8), intent(in)  :: nucThresh,coreThresh,bondDiff
      mNucThresh = nucThresh
      mCoreThresh = coreThresh
      mBondDiff = bondDiff
   end subroutine



   subroutine atoms_whatPosition(self,r,resultString,resultTypeCode,resultIndexCode)
   !-------------------------------------------------------------------------------!
      type(atom), intent(in)         :: self(:)
      real(r8), intent(in)             :: r(3)
      character(len=*),intent(inout),optional :: resultString
      integer, intent(inout), optional        :: resultTypeCode
      integer, intent(inout), optional        :: resultIndexCode
      integer i
      real(r8), allocatable            :: rDist(:)
      integer, allocatable           :: idx(:)

      allocate(rDist(mNAtoms),idx(mNAtoms))
      do i=1,mNAtoms
         idx(i) = i
      end do

      do i=1,mNAtoms
         rDist(i) = sqrt( (self(i)%cx-r(1))**2 + (self(i)%cy-r(2))**2 + &
                                    (self(i)%cz-r(3))**2 )
      end do

      call quickSortIndex(isGreaterEqual,isSmallerEqual,idx)

      if (rDist(idx(1)) < mNucThresh) then
         if (present(resultTypeCode)) then
            resultTypeCode = 10000 + self(idx(1))%elemIdx
         else if (present(resultIndexCode)) then
            resultIndexCode = 10000 + idx(1)
         else if (present(resultString)) then
            write(resultString,'(2a,i2)') "nucl  ",self(idx(1))%elem,idx(1)
         endif
      else if (self(idx(1))%za > 2 .and. rDist(idx(1)) < mCoreThresh/self(idx(1))%za ) then
         if (present(resultTypeCode)) then
            resultTypeCode = 20000 + self(idx(1))%elemIdx
         else if (present(resultIndexCode)) then
            resultIndexCode = 20000 + idx(1)
         else if (present(resultString)) then
            write(resultString,'(2a,i2)') "core  ",self(idx(1))%elem,idx(1)
         endif
      else if (rDist(idx(1)) + rDist(idx(2)) < mAtomDist(idx(1),idx(2)) + mBondDiff) then
         if (present(resultTypeCode)) then
            resultTypeCode = 30000 + self(idx(1))%elemIdx*100 + self(idx(2))%elemIdx
         else if (present(resultIndexCode)) then
            resultIndexCode = 30000 + idx(1)*100 + idx(2)
         else if (present(resultString)) then
            write(resultString,'(2a,i2,2x,a,i2)') "bond1 ",self(idx(1))%elem,idx(1),self(idx(2))%elem,idx(2)
         endif
      else if (rDist(idx(1))**2 + mAtomDist(idx(1),idx(2))**2 < rDist(idx(2)) ) then
         if (present(resultTypeCode)) then
            resultTypeCode = 40000 + self(idx(1))%elemIdx*100 + self(idx(2))%elemIdx
         else if (present(resultIndexCode)) then
            resultIndexCode = 40000 + idx(1)*100 + idx(2)
         else if (present(resultString)) then
            write(resultString,'(2a,i2,2x,a,i2)') "bond2 ",self(idx(1))%elem,idx(1),self(idx(2))%elem,idx(2)
         endif
      else
         if (present(resultTypeCode)) then
            resultTypeCode = 50000 + self(idx(1))%elemIdx
         else if (present(resultIndexCode)) then
            resultIndexCode = 50000 + idx(1)
         else if (present(resultString)) then
            write(resultString,'(2a,i2)') "lonep ",self(idx(1))%elem,idx(1)
         endif
      end if

      deallocate(rDist,idx)

   contains

      logical function isGreaterEqual(i,j)
         integer, intent(in) :: i,j
         isGreaterEqual = rDist(i) >= rDist(j)
      end function isGreaterEqual

      logical function isSmallerEqual(i,j)
         integer, intent(in) :: i,j
         isSmallerEqual = rDist(i) <= rDist(j)
      end function isSmallerEqual

   end subroutine atoms_whatPosition


   function atoms_whatType(self,x,y,z) result(resultv)
   !-------------------------------------------------!
      ! count the "type" 1..5 as defined in atoms_whatPosition for all electrons in coords "cd"
      ! return count of type as resultv(type)
      type(atom), intent(in)         :: self(:)
      real(r8), intent(in)             :: x(:),y(:),z(:)
      integer                        :: resultv(5)
      real(r8) r(3)
      integer i,rtc,t

      call assert(size(x)>0,"atoms_whatType: coord arg must be allocated")
      resultv = 0
      do i=1,size(x)
         r = (/ x(i), y(i), z(i) /)
         call atoms_whatPosition(self,r,resultTypeCode=rtc)
         t = rtc / 10000
         resultv(t) = resultv(t) + 1
      end do
   end function atoms_whatType


   subroutine atoms_spaceWarpTransform(self,other,x,y,z)
   !---------------------------------------------------!
      ! implements the space warp transformation from Cyrus Umrigar, e.g. PRB 61, R16291
      ! electron coordinates x,y,z are transformed from "other" nuclei to "self" nuclei
      ! in omega, any decaying function might be used. Currently x**(-4)
      type(atom), intent(in)  :: self(:)
      type(atom), intent(in)  :: other(:)
      real(r8), intent(inout)   :: x(:),y(:),z(:)
      real(r8)                  :: xnew(size(x)),ynew(size(x)),znew(size(x))
      real(r8)                  :: omega(size(self))
      real(r8)                  :: rai(size(self),size(x))
      integer i,a
      real(r8) sum
      real(r8), parameter       :: RMIN = 1.d-8
      call assert(size(self)==size(other),"atoms_spacewarptransform: sizes do not match")
      do i=1,size(x)
         do a=1,size(other)
            rai(a,i) = sqrt( (other(a)%cx-x(i))**2 + (other(a)%cy-y(i))**2 + (other(a)%cz-z(i))**2 )
            rai(a,i) = max(rai(a,i),RMIN)
         end do
      end do
      !!!write(999,*) ' space warp transformation:'
      do i=1,size(x)
         sum = 0
         do a=1,size(self)
            sum = sum + rai(a,i)**(-2)
         end do
         do a=1,size(self)
            omega(a) = rai(a,i)**(-2) / sum
         end do
         xnew(i) = x(i)
         ynew(i) = y(i)
         znew(i) = z(i)
         do a=1,size(self)
            xnew(i) = xnew(i) + (self(a)%cx - other(a)%cx) * omega(a)
            ynew(i) = ynew(i) + (self(a)%cy - other(a)%cy) * omega(a)
            znew(i) = znew(i) + (self(a)%cz - other(a)%cz) * omega(a)
         end do
         !!!write(999,*) i
         !!!write(999,'(3g20.11)') x(i),y(i),z(i)
         !!!write(999,'(3g20.11)') xnew(i),ynew(i),znew(i)
         !!!write(999,'(10g15.7)') maxval(omega(:))
      end do
      x = xnew; y = ynew; z = znew
   end subroutine atoms_spaceWarpTransform


! internal functions


   integer function getPSEIdx(c)
   !---------------------------!
   character(len=2), intent(in) :: c
   integer i

   do i=0,PSEMAX
      if (c==pse(i)) exit
   enddo

   if (i>PSEMAX) call error('getPSEIdx: Element not yet known')
   getPSEIdx = i
   end function getPSEIdx



   function getPSESymbol(i)
   !----------------------!
   character(len=2) :: getPSESymbol
   integer, intent(in) :: i

   getPSESymbol = pse(i)
   end function getPSESymbol



   subroutine calcSameAtomIndices(self)
   !----------------------------------!
   type(atom), intent(inout) :: self(:)
   integer i,f,idx
   integer sa(size(self))

   call assert(size(self)>0,'atom: calcJastrowIndices, empty atom array')
   sa = 0
   idx = 0
   do i=1,size(self)
      f = findSAValue(self(i)%elemIdx)
      if (f==0) then
         idx = idx + 1
         sa(idx) = self(i)%elemIdx
      end if
   end do

   ! now assign identical atoms the same index
   do i=1,size(self)
      f = findSAValue(self(i)%elemIdx)
      call assert(f>0,'calcSameAtomIndices: internal error')
      self(i)%sa = f
   end do

   contains
      integer function findSAValue(k)
         integer k,i
         do i=1,idx
            if (sa(i)==k) exit
         end do
         if (i>idx) then
            findSAValue = 0
         else
            findSAValue = i
         end if
      end function findSAValue
   end subroutine calcSameAtomIndices

end module atom_m
