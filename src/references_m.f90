! Copyright (C) 2012-2014, 2017-2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module references_m

   ! This module is the "public" interface to collecting reference electron configurations,
   ! in particular the local maxima of Psi**2.
   ! Since there are many useful ways to collect and sort the numerous maxima
   ! the implementation uses a base class and derived classes for different "modes"
   ! The "factory" method "create_refContainer" initializes the correct class
   ! and returns a polymorphic base class pointer allowing to access
   ! derived class methods via the base class pointer
   ! this module stores reference electron configurations
   ! reference positions (e.g. local maxima of psi**2) can be added (append or insert)
   ! and are internally sorted according to different criteria
   ! as output, .ref-files are generated containing the sorted reference configurations
   ! this module creates an abstract base class (referenceContainer) and a "create" method
   ! that creates the actual class. Only the base class methods need to be used in applications
   ! (polymorphic objects)
   ! The basis data structure is a list of list of reference electron configurations, requiring
   ! two indices for access

   ! this module also contains the reference_check and reference_optimize subroutines
   ! that are called directly from amolqc.f90
   !
   ! parallelism: this container should be kept on the MASTER (only)
   !   collect data (GATHER) on the MASTER and use insert/append for all gathered data
   !   no combining/merging is done in this module! local statistics is used only!

   use kinds_m, only: r8
   use global_m, only: getNElec,iul,MASTER,ne,logmode
   !!!use nextpsimaxModule, only: nextpsimax use minimizer_w_sing_module instead
   use utils_m, only: tokenize, readFileLocal, readFileParallel, findInList
   use parsing_m, only: getstra, getinta, getdbla
   use refBase_m
   use refSimple_m
   use refStr_m
   use refStrPerm_m
   use refVal_m
   use refValStr_m
   use refPos_m
   use refCtr_m
   use refDomain_m
   use refUtils_m, only: calcElocRef, calcEigenRef

   implicit none

   private
   public :: referenceContainer, create_refContainer, references_check, references_check1, references_optimize


contains

   subroutine create_refContainer(rc_p,typ,listLength,doSortFreq,elemLength,valueThreshold,simThreshold,sameThreshold, &
                                  refFile,exclFile,ire,simmaxThreshold)
      class(referenceContainer), pointer       :: rc_p
      character(len=3), intent(in)             :: typ
      integer, intent(in)                      :: listLength
      logical, intent(in)                      :: doSortFreq
      integer, optional, intent(in)            :: elemLength
      real(r8), optional, intent(in)             :: valueThreshold
      real(r8), optional, intent(in)             :: simThreshold
      real(r8), optional, intent(in)             :: sameThreshold
      real(r8), optional, intent(in)             :: simmaxThreshold
      integer, optional, pointer, intent(in)   :: ire(:,:)
      character(len=40), optional, intent(in)  :: refFile
      character(len=40), optional, intent(in)  :: exclFile
      class(refContainerSimple), pointer :: rcsi_p
      class(refContainerValue), pointer  :: rcv_p
      class(refContainerValStr), pointer :: rcvs_p
      class(refContainerStruct), pointer :: rcs_p
      class(refContainerCtr), pointer :: rcc_p
      class(refContainerStructPerm), pointer :: rcsp_p
      class(refContainerPos), pointer    :: rcp_p
      class(refContainerDom), pointer    :: rcd_p
      real(r8) dist2thr

      if (typ=="sim") then
         allocate(rcsi_p)
         call refcsi_create(rcsi_p,listLength,doSortFreq)
         rc_p => rcsi_p
      else if (typ=="val") then
         allocate(rcv_p)
         if (present(valueThreshold)) then
            call refcv_create(rcv_p,listLength,doSortFreq,valueThreshold)
         else
            call refcv_create(rcv_p,listLength,doSortFreq)
         end if
         rc_p => rcv_p
      else if (typ=="vst") then
         if (.not.present(elemLength)) call abortp("create_refContainer: elemLength missing")
         allocate(rcvs_p)
         if (present(valueThreshold)) then
            if (present(sameThreshold)) then
               call refcvs_create(rcvs_p,listLength,doSortFreq,elemLength,valueThreshold,sameThreshold=sameThreshold)
            else
               call refcvs_create(rcvs_p,listLength,doSortFreq,elemLength,valueThreshold)
            end if
         else
            if (present(sameThreshold)) then
               call refcvs_create(rcvs_p,listLength,doSortFreq,elemLength,sameThreshold=sameThreshold)
            else
               call refcvs_create(rcvs_p,listLength,doSortFreq,elemLength)
            end if
         end if
         rc_p => rcvs_p
      else if (typ=="str") then
         if (.not.present(elemLength)) call abortp("create_refContainer: elemLength missing")
         allocate(rcs_p)
         if (.not.present(sameThreshold)) call abortp("create_refContainer: mode str requires sameThreshold arg")
         if (present(exclFile)) then
            if (present(simThreshold)) then 
               call refcs_create(rcs_p,listLength,doSortFreq,elemLength,exclFile,simThreshold=simThreshold, &
                                 sameThreshold=sameThreshold)
            else
               call abortp("create_refContainer: illegal parameter combination given for mode=str")
            end if
         else
            if (present(simThreshold)) then 
               call refcs_create(rcs_p,listLength,doSortFreq,elemLength,simThreshold=simThreshold, &
                                 sameThreshold=sameThreshold)
            else
               call abortp("create_refContainer: illegal parameter combination given for mode=str")
            end if
         end if 
         rc_p => rcs_p
      else if (typ=="stp") then
         if (.not.present(elemLength)) call abortp("create_refContainer: elemLength missing")
         allocate(rcsp_p)
         if (present(simThreshold)) then
            if (present(sameThreshold)) then
               call refcsp_create(rcsp_p,listLength,doSortFreq,elemLength,simThreshold,sameThreshold)
            else
               call refcsp_create(rcsp_p,listLength,doSortFreq,elemLength,simThreshold)
            end if 
         else
            if (present(sameThreshold)) then
               call refcsp_create(rcsp_p,listLength,doSortFreq,elemLength,sameThreshold=sameThreshold)
            else
               call refcsp_create(rcsp_p,listLength,doSortFreq,elemLength)
            end if
         end if
         rc_p => rcsp_p
      else if (typ=="ctr") then
         if (.not.present(elemLength)) call abortp("create_refContainer: elemLength missing")
         allocate(rcc_p)
         if (.not.present(sameThreshold)) call abortp("create_refContainer: mode ctr requires sameThreshold arg")
         if (present(exclFile)) then
            if (present(simThreshold)) then 
               call refcc_create(rcc_p,listLength,doSortFreq,elemLength,refFile,simThreshold,sameThreshold,exclFile)
            else
               call abortp("create_refContainer: illegal parameter combination given for mode=ctr")
            end if
         else
            if (present(simThreshold)) then 
               call refcc_create(rcc_p,listLength,doSortFreq,elemLength,refFile,simThreshold,sameThreshold)
            else
               call abortp("create_refContainer: illegal parameter combination given for mode=str")
            end if
         end if 
         rc_p => rcc_p
      else if (typ=="dom") then
         if (.not.present(elemLength)) call abortp("create_refContainer: elemLength missing")
         allocate(rcd_p)
         if (.not.present(sameThreshold)) call abortp("create_refContainer: mode dom requires sameThreshold arg")
         if (present(exclFile)) then
            if (present(simThreshold)) then 
               call refcd_create(rcd_p,listLength,doSortFreq,elemLength,refFile,simThreshold,sameThreshold,exclFile)
            else
               call abortp("create_refContainer: illegal parameter combination given for mode=ctr")
            end if
         else
            if (present(simThreshold)) then 
               call refcd_create(rcd_p,listLength,doSortFreq,elemLength,refFile,simThreshold,sameThreshold)
            else
               call abortp("create_refContainer: illegal parameter combination given for mode=str")
            end if
         end if 
         rc_p => rcd_p
      else if (typ=="pos") then

         if (.not.present(refFile)) call abortp("create_refContainer: refFile for type pos required")
         if (present(ire).and.present(exclFile)) call abortp("create_refContainer: ire and exclFile exclude each other")

         allocate(rcp_p)

         if (present(ire)) then

            if (present(simThreshold)) then
               if (present(simmaxThreshold)) then
                  call refcp_create(rcp_p,listLength,doSortFreq,refFile,simThreshold,simmaxThreshold,ire=ire)
               else
                  call refcp_create(rcp_p,listLength,doSortFreq,refFile,simThreshold,simThreshold,ire=ire)
               end if
            else
               call refcp_create(rcp_p,listLength,doSortFreq,refFile,ire=ire)
            end if 
 
         else if (present(exclFile)) then

            if (present(simThreshold)) then
               if (present(simmaxThreshold)) then
                  call refcp_create(rcp_p,listLength,doSortFreq,refFile,simThreshold,simmaxThreshold,exclFile=exclFile)
               else
                  call refcp_create(rcp_p,listLength,doSortFreq,refFile,simThreshold,simThreshold,exclFile=exclFile)
               end if
            else
               call refcp_create(rcp_p,listLength,doSortFreq,refFile,exclFile=exclFile)
            end if 

         else ! neither ire nore exclFile given

            if (present(simThreshold)) then
               if (present(simmaxThreshold)) then
                  call refcp_create(rcp_p,listLength,doSortFreq,refFile,simThreshold,simmaxThreshold)
               else
                  call refcp_create(rcp_p,listLength,doSortFreq,refFile,simThreshold,simThreshold)
               end if
            else
               call refcp_create(rcp_p,listLength,doSortFreq,refFile)
            end if 

         end if

         rc_p => rcp_p

      else
         call abortp("create_refContainer: illegal type")
      end if
   end subroutine create_refContainer

   subroutine references_check(lines,nl)
      character(len=120), intent(in) :: lines(:)
      integer, intent(in)            :: nl
      integer, parameter             :: MAXBLOCK=500000
      integer, parameter             :: MAXLEN=60
      character(len=MAXLEN),allocatable :: refLines(:)
      character(len=80)              :: fname,evname
      real(r8)                         :: f,EL,vpot,meanDist,maxDist
      real(r8)                         :: fgrad(getNElec(),3)
      real(r8)                         :: nucThresh,coreThresh,bondThresh,hh,r(3)
      real(r8), pointer                :: eval(:)=>null(),evec(:,:)=>null()
      integer                        :: a, i, j, jj, iflag, iter, n, nrl, verb, iflag1, iflag2, iflag3
      integer                        :: naNuc, nbNuc    ! alpha, beta elecs at nucleus (not varied)
      integer                        :: aNucElec(getNNuc()), bNucElec(getNNuc()) ! nuc a -> elec i (at nuc)
      integer                        :: slist(getNElec())   ! slist(i) == a : elec i is at nuc a (0: not at nuc)
      integer                        :: refpar(2)
      type(refContainerSimple) :: refs
      type(reference), pointer       :: rp,rp0
      type(refl_vlist)               :: rl
      type(refl_vlist), pointer      :: p,p0
      logical writeEV,calcDiff,withExcl
      character(len=20)              :: resultStr
      character(len=40)              :: exclFile
      character(len=120)             :: line
      character(len=20)              :: words(20)
      character(len=1)               :: emode
      integer, parameter :: iu=14
      integer io, nWords, exclMode, nnuc
      integer, allocatable           :: exclist(:)

      if (MASTER) then

         call getstra(lines,nl,'ref_file=',fname,iflag)
         if (iflag /= 0) call abortp('$analyze_refs: ref_file=fname required')
         nucThresh = 0.01d0
         call getdbla(lines,nl,'nuc_thresh=',nucThresh,iflag1)
         bondThresh = 0.4d0
         call getdbla(lines,nl,'bond_thresh=',bondThresh,iflag2)
         coreThresh = 2.d0
         call getdbla(lines,nl,'core_thresh=',coreThresh,iflag3)
         hh = 0.001d0
         call getdbla(lines,nl,'h=',hh,iflag)
         refpar(1) = 0
         call getinta(lines,nl,'ref_n=',refpar(1),iflag)
         refpar(2) = 100
         call getinta(lines,nl,'maxrefs=',refpar(2),iflag)
         writeEV = finda(lines,nl,'eigenvec')
         calcDiff = finda(lines,nl,'calc_diff')
         withExcl = .false.
         call getstra(lines,nl,'excl_file=',exclFile,iflag2)
         if (iflag2==0) withExcl = .true.
         if (withExcl) then
            ! read file containing coded list of excluded position (as determined by atoms_m.f90:atoms_whatPosition)
            open(iu,file=exclFile,status='old',iostat=io)
            call assert(io==0,' exclusion file '//trim(exclFile)//' does not exist')
            read(iu,'(a)') line
            call tokenize(line,words,nWords)
            if (nWords==0) call abortp("refcs_create: illegal format in exclusion file")
            read(words(1),*) n
            exclMode = 0
            if (nWords>=2) then
               read(words(2),'(a1)') emode
               if (emode=='i') exclMode = 1
            end if 
            if (nWords >= 5) then
               read(words(3),*) nucThresh
               read(words(4),*) coreThresh
               read(words(5),*) bondThresh
               call atoms_initPositionThresholds(nucThresh, coreThresh, bondThresh)
            end if 
            allocate(exclist(n))        
            do i = 1, n
               read(iu,*) exclist(i)
            end do
            close(iu)
            call sort(exclist)
         end if 

         verb = logmode
         call getinta(lines, nl, 'verbose=', verb, iflag)

         allocate(refLines(MAXBLOCK))
         call readFileLocal(fname, refLines, nrl)
         call refs%readFromFile(refpar, refLines, nrl, .false.)
         deallocate(refLines)

         if (iflag1==0.or.iflag2==0.or.iflag3==0) then
            call atoms_initPositionThresholds(nucThresh, coreThresh, bondThresh)
         end if 

         if (verb >= 1) then
            write(iul,'(/a/)') '  * * *  analyze references  * * *'
            write(iul,*) ' reference file = ',fname
            write(iul,'(a,i5,a,f13.5)') ' with size ',refs%getSize()
            if (refpar(1)==0) then
               write(iul,'(a)') ' analyzing all (first) references'
            else
               write(iul,'(a,i4,a)') ' analyzing all references of ',refpar(1),'-th entry'
            end if
            write(iul,'(a15,f13.5)') ' diff h=',hh
            write(iul,'(3(a15,f13.5))') ' nuc_thresh=',nucThresh,' core_thresh=',coreThresh, &
                                           ' bond_thresh=',bondThresh
            write(iul,'(/a)') ' position of nuclei (in bohr):'
            do a = 1, getNNuc()
               write(iul,'(i2,1x,a2,3f12.5)') a, atoms(a)%elem, atoms(a)%cx, atoms(a)%cy, atoms(a)%cz
            end do
         end if

         if (writeEV) then
            evname = trim(basename)//'.ev'
            open(iu,file=evname)
            write(iu,'(i5,a)') ncenter, " nuclei:"
            do i = 1, getNNuc()
               write(iu,'(i4,1x,a2,1x,3f14.7)') i, atoms(i)%elem, atoms(i)%cx, atoms(i)%cy, atoms(i)%cz
            end do
            write(iu,'(i5,1x,a)') refs%getSize(), ' references:'
         end if

         naNuc = 0; nbNuc = 0
         aNucElec = 0; bNucElec = 0
         do n = 1, refs%getSize()
            call refs%getElemPtr(n, p)
            rp => p%elem(1)
            call calcElocRef(rp, f, EL, vpot, fgrad)
            if (verb >= 1) write(iul,'(/a,i4,3(a,f12.5))') 'ref', n, '  f=', f, &
               '  E_loc=', EL, '  V_pot=', vpot
            !!call findNucElecs(thresh,rp%x,rp%y,rp%z,naNuc,aNucElec,nbNuc,bNucElec)

            if (verb >= 2) then
               write(iul,*) ' coords:'
               do i = 1, getNElec()
                  r = (/ rp%x(i), rp%y(i), rp%z(i) /)
                  call atoms_whatPosition(atoms, r, resultStr)
                  !!if (i <= getNAlpha()) then
                  !!   a = findInList(aNucElec,i)
                  !!else
                  !!   a = findInList(bNucElec,i)
                  !!end if
                  !!if (a > 0) then
                  !!   write(iul,'(i5,3F12.5,2x,a,i2,2a)') i,rp%x(i),rp%y(i),rp%z(i), &
                  !!                                       trim(atoms(a)%elem)//'(',a,')  :',resultStr
                  !!else
                  write(iul,'(i5,3F12.5,2a)') i, rp%x(i), rp%y(i), rp%z(i), ' : ', resultStr
                  !!end if
               end do

               if (calcDiff) then
                  write(iul,*)
                  do j = 1, n - 1
                     call refs%getElemPtr(j, p0)
                     rp0 => p0%elem(1)
                     if (withExcl) then
                        call calcRefDifference(rp0,rp,2,meanDist,maxDist,exclist=exclist,exclmode=exclMode)
                     else
                        call calcRefDifference(rp0,rp,2,meanDist,maxDist)
                     endif
                     write(iul,'(a,2i3,F12.5)') 'max dist:', n, j, maxDist
                  end do
                  write(iul,*)
               end if

               if (verb >= 3) then
                  write(iul,*) ' grad:'
                  do i = 1, getNElec()
                     write(iul,'(i5,3F12.5)') i, fgrad(i,1), fgrad(i,2), fgrad(i,3)
                  end do
               end if
            end if
            if (writeEV) then
               write(iu,'(a,i5,a,f15.5,a,i6,4i4)') "ref:",n,"  1 F(ref):", f,"  found:",rp%count, &
                  naNuc, getNAlpha(), nbNuc, getNelec()
               write(iu,'(i5)') size(rp%x)
               do i = 1, size(rp%x)
                  write(iu,'(3f13.6)') rp%x(i), rp%y(i), rp%z(i)
               end do
            end if
            call calcEigenRef(rp, verb, eval, evec, slist, nnuc, hh, nucThresh)
            call assert(associated(eval) .and. associated(evec),'(references_check): eigenvalue calculation failed')
            if (verb >= 1) then
               write(iul,*) ' eigenvalues:'
               do i = 1, size(eval), 5
                  write(iul,'(5g15.5)') eval(i:min(i+4,size(eval)))
               end do
            end if
            if (writeEV) then
               write(iu,'(i5,a)') size(eval)," eigenvectors:"
               do i = 1, size(eval)
                  write(iu,'(a,2i5,f15.5)') "ev"//":", n, i, eval(i)
                  write(iu,'(i5)') getNElec()
                  j = 1
                  do jj = 1, getNElec()
                     if (slist(jj) > 0) then
                        write(iu, '(3f13.6)') 0.d0, 0.d0, 0.d0
                     else
                        write(iu,'(3f13.6)') evec(j : j+2, i)
                        j = j + 3
                     end if
                  end do
               end do
            end if
         end do
         call refs%destroy()
         if (writeEV) close(iu)
      end if
   end subroutine references_check

   subroutine references_check1(lines,nl)
      character(len=120), intent(in) :: lines(:)
      integer, intent(in)            :: nl
      integer                        :: nmax,mmax,iflag,iflag1,iflag2,iflag3
      integer                        :: verb,ref,noel,nomax
      integer                        :: idx,nre,nrefs,i,j,n,a
      integer, pointer               :: ire(:,:) => null()
      real(r8)                         :: ftol,tolSim,tolSame,tolDist2,hdist,f
      real(r8)                         :: nucThresh,coreThresh,bondThresh
      character(len=3)               :: maxMode
      character(len=40)              :: exclFile,refFile
      character(len=120)             :: line
      character(len=20)              :: words(20)
      integer, parameter :: iu=14
      integer io,nWords,exclMode
      integer, allocatable           :: exclist(:)
      real(r8) x(ne),y(ne),z(ne)
      type(reference) :: r
      class(referenceContainer), pointer  :: rc_p => null()         ! list of stored references

      if (MASTER) then

         call getstra(lines,nl,'ref_file=',refFile,iflag)
         if (iflag /= 0) call abortp('$analyze_refs: ref_file=fname required')
         nmax = 30
         call getinta(lines,nl,'nmax=',nmax,iflag)          ! max # (main list) entries ref file
         mmax = 5
         call getinta(lines,nl,'mmax=',mmax,iflag)        ! max # of sublists in ref fie
         ftol = 1.d-3
         call getdbla(lines,nl,'tol_fctn=',ftol,iflag)      ! same function value if abs diff is smaller
         tolSim = 1.d-1
         call getdbla(lines,nl,'tol_sim=',tolSim,iflag)      ! similar structure if max distance is smaller
         tolSame = 1.d-2
         call getdbla(lines,nl,'tol_same=',tolSame,iflag)    ! same structure if max distance is smaller
         tolDist2 = tolSim
         call getdbla(lines,nl,'tol_simmax=',tolDist2,iflag) ! type=pos: same position if distance is smaller
         hdist = -1.d0
         call getdbla(lines,nl,'H_dist=',hdist,iflag)
         nucThresh = 0.01d0
         bondThresh = 0.4d0
         coreThresh = 2.d0

         call getstra(lines,nl,'max_mode=',maxMode,iflag)
         if (iflag/=0) call abortp("references_check1 requires max_mode")
         if  (maxMode=="val") then
            call create_refContainer(rc_p,maxMode,nmax,.false.,elemLength=1,valueThreshold=ftol)
         else if (maxMode=="vst") then
            call create_refContainer(rc_p,maxMode,nmax,.false.,elemLength=mmax,valueThreshold=ftol, &
                                     sameThreshold=tolSame)
         else if (maxMode=="str") then
            call create_str()     ! internal subroutine
         else if (maxMode=="stp") then
            call create_refContainer(rc_p,maxMode,nmax,.false.,elemLength=mmax,simThreshold=tolSim, &
                                     sameThreshold=tolSame)
         else if (maxMode=="pos") then
            call create_pos()     ! internal subroutine
         else
            call abortp("$analyze_refs: illegal value for max_mode=[val|vst|str|stp|pos]")
         end if

         ! default values above "maxMode". Values here have preference over values in excl_file
         call getdbla(lines,nl,'nuc_thresh=',nucThresh,iflag1)
         call getdbla(lines,nl,'bond_thresh=',bondThresh,iflag2)
         call getdbla(lines,nl,'core_thresh=',coreThresh,iflag3)
         if (iflag1==0.or.iflag2==0.or.iflag3==0) then
            call atoms_initPositionThresholds(nucThresh,coreThresh,bondThresh)
         end if 

         verb = logmode
         call getinta(lines,nl,'verbose=',verb,iflag)
         if (verb >= 1) then
            write(iul,'(/a/)') '  * * *  analyze references  * * *'
            write(iul,*) ' reference file = ',refFile
            write(iul,'(3(a15,i5,8x))') 'nmax=',nmax,'mmax=',mmax
            write(iul,'(3(a15,f13.5))') 'tol_fctn=',ftol,'tol_sim=',tolSim
            write(iul,'(3(a15,f13.5))') 'tol_same=',tolSame,'tol_simmax=',tolDist2
            write(iul,'(3(a15,f13.5))') ' nuc_thresh=',nucThresh,' core_thresh=',coreThresh, &
                                           ' bond_thresh=',bondThresh
            write(iul,'(/a)') ' position of nuclei (in bohr):'
            do a=1,getNNuc()
               write(iul,'(i2,1x,a2,3f12.5)') a,atoms(a)%elem,atoms(a)%cx,atoms(a)%cy,atoms(a)%cz
            end do
         end if

         call r%new(ne)
         !! read ref file and insert all entries
         open(iu,file=refFile,status='old',iostat=io)
         call assert(io==0,' ref file '//trim(refFile)//' does not exist')
         read(iu,*) n   ! ncenter
         do i=1,n
            read(iu,*) j  ! index,elem,x,y,z
         end do
         ref=1
         read(iu,*) nomax
         do
            read(iu,'(A)',iostat=io) line
            call tokenize(line,words,n)
            read(words(5),*) f
         if (io /= 0) exit
            read(iu,*) noel
            call assert(noel==ne,"refcp_create: wrong number of electrons in ref file")
            do i=1,noel
               read(iu,*) x(i),y(i),z(i)
            enddo
            call r%set(x,y,z,f)
            call rc_p%insert(r)
            ref = ref+1
         enddo
         close(iu)
         call r%destroy()

         if (verb >= 2) write(iul,'(/i5,2a/)') ref," references read from file ",refFile
         call rc_p%writeShortList(iul,'Max')
         call rc_p%destroy()

      end if

   contains

         subroutine create_str()
            call getstra(lines,nl,'excl_file=',exclfile,iflag2)
            if (iflag2 == 0) then
               call create_refContainer(rc_p,maxMode,nmax,.false.,exclFile=exclfile,elemLength=mmax, &
                                 simThreshold=tolSim,sameThreshold=tolSame)
            else
               call create_refContainer(rc_p,maxMode,nmax,.false.,elemLength=mmax, &
                                 simThreshold=tolSim,sameThreshold=tolSame)
            end if 
         end subroutine create_str

         subroutine create_pos()
            call getstra(lines,nl,'ref_file=',reffile,iflag)
            call assert(iflag == 0,'$find_maxima: max_mode=pos requires ref_file')
            nre = 0
            call getinta(lines,nl,'ignore_ref_elecs=',nre,iflag)
            if (iflag==0 .and. nre>0) then
               do idx=1,nl
                  if (index(lines(idx),'ignore_ref_elecs')>0) exit
               end do
               call assert(nl > idx,'$find_maxima: not enough lines for ignore_ref_elecs')
               idx = idx + 1
               read(lines(idx),*) nrefs
               call assert(nrefs>0,"$find_maxima: line after ignore_ref_elecs must contain # of references")
               call assert(nl >= idx+nrefs,'$find_maxima: not enough lines for ignore_ref_elecs')
               allocate(ire(nre,nrefs))
               do i=1,nrefs
                  read(lines(idx+i),*) (ire(j,i),j=1,nre)
                  call sort(ire(:,i))
               end do
            else
               allocate(ire(1,1))
               ire(1,1) = 0
            end if
            call getstra(lines,nl,'excl_file=',exclfile,iflag2)

            if (iflag2 == 0) then ! excl_file given

               call create_refContainer(rc_p,maxMode,nmax,.false.,refFile=reffile,exclFile=exclfile, &
                        simThreshold=tolSim,simmaxThreshold=tolDist2)

            else if (ire(1,1) /= 0) then ! ignore reference electron list given

               call create_refContainer(rc_p,maxMode,nmax,.false.,refFile=reffile,ire=ire, &
                                   simThreshold=tolSim,simmaxThreshold=tolDist2)

            else  ! neither ire nore exclFile given

               call create_refContainer(rc_p,maxMode,nmax,.false.,refFile=reffile, &
                                   simThreshold=tolSim,simmaxThreshold=tolDist2)

            end if
         end subroutine create_pos  

   end subroutine references_check1



   subroutine references_optimize(lines,nl)
      character(len=120), intent(in) :: lines(:)
      integer, intent(in)            :: nl
      character(len=80)              :: fname,newName
      integer, parameter             :: MAXBLOCK=50000
      integer, parameter             :: MAXLEN=60
      character(len=MAXLEN)          :: refLines(MAXBLOCK)
      real(r8)                         :: x(getNElec()), y(getNElec()), z(getNElec())
      real(r8)                         :: hThresh,EL,vpot,f,thresh
      real(r8)                         :: eps(3)
      real(r8)                         :: fgrad(getNElec(),3)
      integer                        :: a,i,iflag,n,verb,maxiter,nCore,ierr,nrl
      integer                        :: naNuc, nbNuc    ! alpha, beta elecs at nucleus (not varied)
      integer                        :: aNucElec(getNNuc()), bNucElec(getNNuc()) ! nuc a -> elec i (at nuc)
      integer                        :: refpar(2)
      logical                        :: doSpaceWarpTransform
      type(refContainerSimple)       :: refs
      type(reference), pointer       :: rp
      type(refl_vlist)               :: rl
      type(refl_vlist), pointer      :: p

      call abortp("reference_optimize: requires update to minimizer class")

      if (MASTER) then

         call getstra(lines,nl,'ref_file=',fname,iflag)
         if (iflag /= 0) call abortp('$optimize_refs: ref_file=fname required')
         call getstra(lines,nl,'new_file=',newName,iflag)
         if (iflag /= 0) call abortp('$optimize_refs: new_file=fname required')
         refpar(1) = 0
         call getinta(lines,nl,'ref_n=',refpar(1),iflag)
         refpar(2) = 100
         call getinta(lines,nl,'maxrefs=',refpar(2),iflag)
         verb = logmode
         call getinta(lines,nl,'verbose=',verb,iflag)
         doSpaceWarpTransform = .true.
         if (finda(lines,nl,'no_transform')) doSpaceWarpTransform = .false.

         call readFileLocal(fname,refLines,nrl)
         call refs%readFromFile(refpar,refLines,nrl,doSpaceWarpTransform)

         !!!call refs%writeSimpleList("test","test")

         call getinta(lines,nl,'max_detail=',verb,iflag)
         maxiter = 100
         call getinta(lines,nl,'max_iter=',maxiter,iflag)
         hThresh = -1.d0
         call getdbla(lines,nl,'H_dist=',hThresh,iflag)
         eps(1) = 1.d-5
         call getdbla(lines,nl,'max_grad=',eps(1),iflag)
         eps(2) = 1.d-4
         call getdbla(lines,nl,'max_dist=',eps(2),iflag)
         eps(3) = 1.d-5
         call getdbla(lines,nl,'max_f=',eps(3),iflag)
         thresh = 0.002d0
         call getdbla(lines,nl,'nuc_thresh=',thresh,iflag)

         if (verb >= 1) write(iul,'(/a/)') '  * * *  optimize references  * * *'

         if (verb >= 2) then
            write(iul,*) ' reference file = ',fname
            write(iul,'(a,i5,a,f13.5)') ' with size ',refs%getSize()
            if (refpar(1)==0) then
               write(iul,'(a)') ' analyzing all (first) references'
            else
               write(iul,'(a,i4,a)') ' analyzing all references of ',refpar(1),'-th entry'
            end if
            write(iul,'(3X,A)') 'using LBFGS-B'
            write(iul,'(3X,A)') 'maximization with 1s core electrons fixed at nucleus'
            write(iul,'(1X,A21,I6,9X,A21,F10.4)') 'max_iter =',maxiter
            write(iul,'(2(1X,A21,G12.2,3X))') " conv. grad tol =",eps(1)," conv. dist tol =",eps(2)
            write(iul,'(2(1X,A21,G12.2,3X))') " conv. func tol =",eps(3)," H nuc when dist <",hThresh
            write(iul,'(2(1X,A21,G12.2,3X))') " nuc thresh =",thresh
            write(iul,*)
         end if

         if (verb >=1 .and. verb < 3) then
            write(iul,'(3x,a5,a13,a5,a6)') 'ref #','final f','ierr','nCore'
         end if

         do n=1,refs%getSize()
            call refs%getElemPtr(n,p)
            rp => p%elem(1)
            !!!!!! use minimizer class ! call nextpsimax(rp%x,rp%y,rp%z,maxiter,hThresh,eps,verb,f,nCore,ierr)
            rp%count = 0
            rp%f = f
            if (verb >=1 .and. verb < 3) then
               write(iul,'(5x,i3,f13.6,3x,i2,3x,i3)') n,f,ierr,nCore
            end if
            if (verb >= 3 .and. ierr < 4) then
               write(iul,'(a,i4,a,f12.5,2(a,i4))') 'ref ',n,': final f=',f,' ierr=',ierr,' nCore=',nCore
               call calcElocRef(rp,f,EL,vpot,fgrad)
               rp%f = f
               write(iul,*) ' coords:'
               call findNucElecs(thresh, rp%x, rp%y, rp%z, naNuc, nbNuc, aNucElec, bNucElec)
               do i = 1, getNElec()
                  if (i <= getNAlpha()) then
                     a = findInList(aNucElec, i)
                  else
                     a = findInList(bNucElec, i)
                  end if
                  if (a > 0) then
                     write(iul,'(i5,3F12.5,2x,a,i2,a)') i,x(i),y(i),z(i),trim(atoms(a)%elem)//'(',a,')'
                  else
                     write(iul,'(i5,3F12.5)') i,x(i),y(i),z(i)
                  end if
               end do
               write(iul,*) ' grad:'
               do i=1,getNElec()
                  write(iul,'(i5,3F12.5)') i,fgrad(i,1),fgrad(i,2),fgrad(i,3)
               end do
            end if
         end do

         fname = newName(1:len(trim(newName))-4)   ! remove '.ref'
         call refs%writeSimpleList(fName,"Opt")
         call refs%destroy()

      end if
   end subroutine references_optimize


end module references_m
