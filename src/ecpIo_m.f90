! Copyright (C) 2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

MODULE ecpIo_m

   use kinds_m, only: r8
   use error_m
   use global_m
   use wfData_m, only: atoms, basis
   use ecp_m
   use parsing_m
   use utils_m, only: readFileParallel
   implicit none

   private
   public :: initecpparams, ecpinput, ecpinputex, ecpoutput, ecpwriteWF

   ! i/o routines for ECPs

CONTAINS


   subroutine initecpparams(lines,nl,ecp)
   !------------------------------------!
      integer, intent(in)           :: nl
      character(len=*), intent(in)  :: lines(nl)
      type(EcpType), intent(inout)  :: ecp

      integer  iflag,i,n,idx,io,anmbr,nGridPoints
      integer  oldlogmode
      real(r8)   ecpcutvalue
      logical  found

      oldlogmode = logmode
      call getinta(lines,nl,'verbose=',logmode,iflag)

      found = finda(lines,nl,'det_only')
      if (found) call ecp%setDetOnlyLocalisation()
      found = finda(lines,nl,'full_localisation')
      if (found) call ecp%setDetOnlyLocalisation(.false.)

      found = finda(lines,nl,'no_random_rotation')
      if (found) then 
         call ecp%setRandomRotation(.false.)
      else 
         found = finda(lines,nl,'random_rotation')
         if (found) call ecp%setRandomRotation(.true.)
      endif

      found = finda(lines,nl,'full_cutoff')
      if (found) call ecp%setCutoffMode("full_cutoff")
      found = finda(lines,nl,'nonlocal_cutoff')
      if (found) call ecp%setCutoffMode("nonlocal_cutoff")
      found = finda(lines,nl,'no_cutoff')
      if (found) call ecp%setCutoffMode("no_cutoff")

      ! careful: setCutoffThreshold _calculates_ all cutoff values
      ! depending on the cutoff mode! default is nonlocal_cutoff
      call getdbla(lines,nl,'cutoff=',ecpcutvalue,iflag)
      if (iflag == 0) then
         call ecp%setCutoffThreshold(ecpcutvalue)
      endif

      call getinta(lines,nl,'grid_points=',n,iflag)
      if (iflag == 0) call ecp%setGrid(n)

      call getinta(lines,nl,'grid_point_list=',n,iflag)
      idx = ifinda(lines,nl,'grid_point_list=')
      if (iflag == 0) then
         do i=1,n
            read(lines(idx+i),*,iostat=io) anmbr,nGridPoints
            if (io /= 0) call abortp('initecpparams: grid point list (atom, grid_points) missing. atom: idx in geom')
            call ecp%setGrid(nGridPoints,atom=anmbr)
         enddo
      endif

      if (MASTER .and. logmode >= 2) call outputParams(ecp)

      logmode = oldlogmode
   end subroutine initecpparams


   subroutine outputParams(ecp)
   !--------------------------!
      !!!type(EcpType), intent(in)  :: ecp     !!! getCutOffMode and getPseudoAtom -> function, then :in
      type(EcpType), intent(inout)  :: ecp

      integer  p
      real(r8) tmp
      logical ltmp
      character(len=30) s
      type(PseudoAtomData), allocatable :: pa(:)

      write(iul,'(a)',advance='no') ' ecp parameters: '

      ltmp = ecp%isDetOnlyLocalisation()
      if (ltmp) then
         write(iul,'(a)',advance='no') 'det_only '
      else
         write(iul,'(a)',advance='no') 'full_localisation '
      endif

      ltmp = ecp%isRandomRotation()
      if (ltmp) then
         write(iul,'(a)',advance='no') 'random_rotation '
      else
         write(iul,'(a)',advance='no') 'no_random_rotation '
      endif

      call ecp%getCutoffMode(s)
      write(iul,'(a)') trim(s)

      tmp = ecp%getCutoffThreshold()
      if (tmp > 0.d0) then
         write(iul,'(/a,g20.5)') ' cutoff threshold =',tmp
      endif
      ! careful: setCutoffThreshold _calculates_ all cutoff values
      ! depending on the cutoff mode! default is nonlocal_cutoff
      write(iul,'(/a)') ' pseudo atoms:'
      write(iul,'(a)') ' idx  atom  grid points  cutoff distances (A) '
      write(iul,'(a)') '-------------------------------------------- '
      call ecp%getPseudoAtoms(pa)
      do p=1,size(pa)
         if (pa(p)%cutoff > 0.d0) then
            write(iul,'(i4,2x,a4,6x,i4,6x,f10.5)') p,atoms(pa(p)%a)%elem,pa(p)%nGridPoints,pa(p)%cutoff*bohr2angs
         else
            write(iul,'(i4,2x,a4,6x,i4,10x,a)') p,atoms(pa(p)%a)%elem,pa(p)%nGridPoints,'n/a'
         endif
      enddo
      write(iul,*)
   end subroutine outputParams



   subroutine ecpinput(lines,nl,ecp)
   !--------------------------------
      ! ecpinput reads ECPs $ecp block (containing sequentially all pseudo atoms)
      character(len=*), intent(in)  :: lines(:)  ! lines of $ecp block
      integer, intent(in)           :: nl        ! actual # of lines
      type(EcpType), intent(inout)  :: ecp

      integer idx, io, lmax, p, a, nPseudoAtoms
      type(PseudoAtomData), allocatable :: pa(:)

      idx = 2
      read(lines(idx),*) nPseudoAtoms
      allocate(pa(nPseudoAtoms))
      idx = idx + 1

      do p = 1, nPseudoAtoms
         read(lines(idx),*,iostat=io) pa(p)%a, pa(p)%nCoreElecs, pa(p)%lmax
         if (io /= 0) call abortp('ecpinput: incorrect ECP format (new format expected)')
         if (pa(p)%a > size(atoms)) call abortp('ecpinput: illegal atom given')
         if (pa(p)%nCoreElecs <= 0) call abortp('ecpinput: # core electrons must be positive')
         if (pa(p)%lmax <= 0) call abortp('ecpinput: lmax must be positive')

         ne = ne - pa(p)%nCoreElecs
         a = pa(p)%a
         atoms(a)%za = atoms(a)%za - pa(p)%nCoreElecs
         atoms(a)%ecp = .true.

         lmax = pa(p)%lmax
         allocate(pa(p)%nterms(0:lmax))

         call readECPEntry(lines,idx,lmax,p,pa)
      enddo

      call ecp%setPseudoAtoms(pa)
   end subroutine ecpinput


   subroutine ecpinputex(ecp)
   !-------------------------
      ! ecpinputex reads ECPs from ECP library file
      type(EcpType), intent(inout)  :: ecp
      integer, parameter           :: MAXLINES=5000
      integer, parameter           :: MAXLEN=120
      character(len=MAXLEN)        :: lines(MAXLINES) ! lines array
      integer                      :: nl      ! actual # of lines
      integer a,nPseudoAtoms,p
      logical fileExists
      character(len=180) basispath, ecpfname
      character(len=3) ecpname
      type(PseudoAtomData), allocatable :: pa(:)

      call getAmolqcPath(basispath)
      call assert(len(trim(basispath)) < 176,"ecpinputex: amolqc path length exceeds definition")
      basispath = trim(basispath)//"/bib"

      allocate(pa(size(atoms)))
      nPseudoAtoms = 0

      if (basis == 'diff') then

         ! for each atom read individual ECP or no ECP at all
         p = 0   ! pseudo atom count
         do a = 1, size(atoms)
            ! check if first 3 characters are ECP name, cycle if not
            ecpname = atoms(a)%ba(1:3)
            if (ecpname /= 'LES' .and. ecpname /= 'SBK' .and. ecpname /= 'STU'   &
               .and. ecpname /= 'NEE' .and. ecpname /='BFD' .and. ecpname /='CRE' .and. ecpname /='MIT' &
               .and. ecpname /='MHE') cycle
            ecpfname = trim(basispath)//'/ecp/'//ecpname//'.ecp'
            if (MASTER) then
               inquire(file=ecpfname,exist=fileExists)
               if (.not.fileExists) call abortp('ecpinputex: ecp file not found')
            end if
            call readFileParallel(mytid,ecpfname,lines,nl)

            call internal_readecpfile()

         enddo
         nPseudoAtoms = p

      else if (basis /= 'default') then
         ! use same ECP and basis for all atoms
         ecpname = basis(1:3)

         if (ecpname == 'LES' .or. ecpname == 'SBK' .or. ecpname == 'STU'    &
            .or. ecpname == 'NEE' .or. ecpname =='BFD' .or. ecpname =='CRE' .or. ecpname =='MIT' &
            .or. ecpname == 'MHE') then
            ecpfname = trim(basispath)//'/ecp/'//ecpname//'.ecp'
            if (MASTER) then
               inquire(file=ecpfname, exist=fileExists)
               if (.not.fileExists) call abortp('(ecpinputex): ecp file not found')
            end if
            call readFileParallel(mytid,ecpfname,lines,nl)

            p = 0   ! pseudo atom count
            do a = 1, size(atoms)
               call internal_readecpfile()
            enddo
            nPseudoAtoms = p
         endif

      endif

      if (nPseudoAtoms > 0) then
         call ecp%setPseudoAtoms(pa(1:nPseudoAtoms))
      endif

   contains

      subroutine internal_readecpfile()
         integer idx, cha, nCoreElecs, lmax
         idx = 1
         do
            if (lines(idx)(1:4) == "####") then
               idx = idx + 1
               read(lines(idx),*) cha
               if (cha == 0) call abortp('(ecpinputex): ECP not found')
               read(lines(idx),*) cha,nCoreElecs,lmax
               if (cha == atoms(a)%elemIdx) then
                  ! correct entry for element found
                  p = p + 1
                  pa(p)%a = a
                  pa(p)%nCoreElecs = nCoreElecs
                  pa(p)%lmax = lmax
                  allocate(pa(p)%nterms(0:lmax))

                  ne = ne - nCoreElecs
                  atoms(a)%za = atoms(a)%za - nCoreElecs
                  atoms(a)%ecp = .true.

                  call readECPEntry(lines,idx,lmax,p,pa)

                  if (lines(idx)(2:5) == '****') exit
               endif
            else
               idx = idx + 1
            endif
         enddo
      end subroutine internal_readecpfile

   end subroutine ecpinputex


   subroutine readECPEntry(lines,idx,lmax,p,pa)
   !-------------------------------------------
      character(len=*), intent(in)        :: lines(:)  ! lines of ecp file
      integer, intent(inout)              :: idx       ! start line containing: a, nCoreElecs, lmax
                                                       ! updated on exit for next entry
      integer, intent(in)                 :: lmax      ! current lmax
      integer, intent(in)                 :: p         ! current pseudo atom
      type(PseudoAtomData), intent(inout) :: pa(:)     ! pseudo atom list         
      integer maxnt, nterms, idx1, nlk, l, n
      real(r8) alk, blk


      ! find max number of terms
      maxnt = 0
      read(lines(idx+1),*) nterms
      maxnt = max(maxnt,nterms)
      idx1 = idx + nterms + 2
      do l = 0, lmax - 1
         read(lines(idx1),*) nterms
         maxnt = max(maxnt,nterms)
         idx1 = idx1 + nterms + 1
      enddo

      allocate(pa(p)%alk(0:lmax,maxnt), pa(p)%blk(0:lmax,maxnt), pa(p)%nlk(0:lmax,maxnt))

      ! read lmax (local) component
      read(lines(idx+1),*) pa(p)%nterms(lmax)
      do n = 1, pa(p)%nterms(lmax)
         read(lines(idx+n+1),*) alk, nlk, blk
         pa(p)%alk(lmax,n) = alk
         pa(p)%nlk(lmax,n) = nlk
         pa(p)%blk(lmax,n) = blk
      enddo
      idx = idx + pa(p)%nterms(lmax) + 2

      ! read (non-local) l - lmax components l = 0 .. lmax-1
      do l = 0, lmax - 1
         read(lines(idx),*) pa(p)%nterms(l)
         do n = 1, pa(p)%nterms(l)
            read(lines(idx+n),*) alk, nlk, blk
            pa(p)%alk(l,n) = alk
            pa(p)%nlk(l,n) = nlk
            pa(p)%blk(l,n) = blk
         enddo
         idx = idx + pa(p)%nterms(l) + 1
      enddo
   end subroutine readECPEntry


   subroutine ecpoutput(iu,ecp)
   !---------------------------
      ! ecpoutput writes ECPs to open file unit 'iu'
      integer, intent(in)           :: iu
      type(EcpType), intent(inout)  :: ecp
      integer n, l, p, a
      type(PseudoAtomData),allocatable :: pa(:)

      write(iu,'(/10X,A)') ' effective core potentials:'

      ! copy pseudo atoms out off ecp
      call ecp%getPseudoatoms(pa)

      do p = 1, size(pa)
         a = pa(p)%a
         write(iu,'(a,i3,a,a2)')  '   atom ', a, ': ',atoms(a)%elem
         write(iu,'(a,i3,a,i2,a,i4)')  ' # core electrons =',pa(p)%nCoreElecs,  &
            '  l_core+1 =', pa(p)%lmax
         write(iu,*) '  V_l = sum_k   a_lk * r^(n_lk - 2) * exp(-b_lk*r^2)'
         write(iu,*) '  l      a_lk      n_lk       b_lk'
         l = pa(p)%lmax
         do n = 1, pa(p)%nterms(l)
            write(iu,'(I5,G16.7,i1,G16.7)') l, pa(p)%alk(l,n), pa(p)%nlk(l,n), pa(p)%blk(l,n)
         enddo
         do l = 0, pa(p)%lmax - 1
            do n = 1, pa(p)%nterms(l)
               write(iu,'(I5,G16.7,i1,G16.7)') l, pa(p)%alk(l,n), pa(p)%nlk(l,n), pa(p)%blk(l,n)
            enddo
         enddo
      enddo

   end subroutine ecpoutput


   subroutine ecpwriteWF(iu,ecp)
   !----------------------------
      ! ecpoutput writes ECPs to open file unit 'iu'
      integer, intent(in)           :: iu
      type(EcpType), intent(inout)  :: ecp
      integer n,l,p
      type(PseudoAtomData),allocatable :: pa(:)

      ! copy pseudo atoms out off ecp
      call ecp%getPseudoatoms(pa)

      write(iu,'(i4)') size(pa)
      do p = 1, size(pa)
         write(iu,'(i3,i3,i2)')  pa(p)%a, pa(p)%nCoreElecs, pa(p)%lmax
         l = pa(p)%lmax
         write(iu, '(i3)') pa(p)%nterms(l)
         do n = 1, pa(p)%nterms(l)
            write(iu,'(g17.10,i3,g17.10)') pa(p)%alk(l,n), pa(p)%nlk(l,n), pa(p)%blk(l,n)
         enddo
         do l = 0, pa(p)%lmax - 1
            do n = 1, pa(p)%nterms(l)
               write(iu,'(g17.10,i3,g17.10)') pa(p)%alk(l,n), pa(p)%nlk(l,n), pa(p)%blk(l,n)
            enddo
         enddo
      enddo

   end subroutine ecpwriteWF


END MODULE ecpIo_m
