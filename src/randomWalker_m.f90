! Copyright (C) 1996, 1998-1999, 2006, 2014-2015, 2017-2018 Arne Luechow
! Copyright (C) 2008 Annett Schwarz
! Copyright (C) 2013, 2015-2016 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

MODULE randomWalker_m

  ! defines type Random Walker
  ! note: ensure consistency of data type (x,y,z consistent with psi,eLocal)
  ! always directly update data members when position changes

! $Id: randomWalker_m.f90, father of Type rwp will now be alloc. wit 1 father[1]
!      per default and set to zero
! $Id: randomWalker_m.f90,v 1.2 2008/02/11 09:37:01 annett Exp $

! $Log: randomWalker_m.f90,v $
! Revision 1.2  2008/02/11 09:37:01  annett
! Energieoptimierung hinzugefuegt
!
! Revision 1.1.1.1  2007/04/25 13:42:20  luechow
! QMC program amolqc. rewritten in 2006. AL
!

  use kinds_m, only: r8
  use elocData_m
  use eloc_m, only: eloc
  use eConfigs_m
  use mpiInterface_m, only: myMPISendDouble, myMPIReceiveDouble

  implicit none

  ! "class variables"
  integer, private       :: dNElecs = 0       ! same size for all RWalker
  integer, private       :: dFathers = 0      ! same size for all (descendent weight gen)
  integer, private       :: dRWCount = 0      ! count created walkers
  integer, private       :: dEE=0              ! number of electrons**2
  integer, private       :: dNE=0              ! ncenter*nelectrons
  integer, private       :: dNcenters=0        ! ncenter
  logical, private       :: dEpart = .false.   ! Enegery Partitioning

  type RandomWalker
     private
     real(r8), pointer     :: x(:) => null()      ! electron position (x(ne))
     real(r8), pointer     :: y(:) => null()      ! electron position (y(ne))
     real(r8), pointer     :: z(:) => null()      ! electron position (z(ne))
     real(r8)              :: w      = 1
     real(r8)              :: phi    = 0          ! psi as phi*exp(u)
     real(r8)              :: u      = 0          !
     real(r8)              :: eLocal = 0          !  H Psi / Psi
     real(r8)              :: eLocalSave=0.0      ! local energy before tmove
     logical             :: tmovedrw=.false.    ! walker was acepted tmove walker
     !real(r8)              :: ecp_nonlocalneg = 0 !  negative part of non-local ecp
     real(r8)              :: ecp_nonlocal = 0 !   part of non-local ecp
     real(r8)              :: vPot   = 0          ! V  (with nucl-nucl repulsion)
     real(r8)              :: vee   = 0           ! V  (with nucl-nucl repulsion)
     real(r8), pointer     :: xdrift(:) => null() ! d_x Psi / Psi
     real(r8), pointer     :: ydrift(:) => null() ! d_y Psi / Psi
     real(r8), pointer     :: zdrift(:) => null() ! d_z Psi / Psi
     integer             :: persist = 0         ! persistency count
     integer, pointer    :: fathers(:) => null()  ! keep ancestors for descendent weighting
     real(r8), pointer     :: Ekin_epart(:) => null()   ! Kinetic energy for each electron
     real(r8), pointer     :: Vne_epart(:,:) => null()    ! Matrix for Nucleus - Electron interaction
     real(r8), pointer     :: Vee_epart(:,:) => null()    !   ..   ..   Electron - Electron .
  end type RandomWalker

  ! overload the assignment operator
  interface assignment(=)
     module procedure rw_assign
  end interface

  interface pos
     module procedure posn
     module procedure pos3n
  end interface pos

  interface resetTo
     module procedure resetTon
     module procedure resetTo3n
  end interface resetTo

CONTAINS

    !-------- class variable methods -------------------------------

    !--------------------------------!
    subroutine setNumberOfElectrons(n)
        !--------------------------------!
        integer, intent(in) :: n
        call assert(dRWCount==0,"setNumberOfElectrons: only before RW creation")
        dNElecs = n
        dEE = n*n
        dNE = n*dNCenters
    end subroutine setNumberOfElectrons

    !--------------------------------!
    subroutine setNumberOfCenters(n)
        !--------------------------------!
        integer, intent(in) :: n
        call assert(dRWCount==0,"setNumberofcenters: only before RW creation")
        dNcenters = n
        dNE = n*dNElecs
    end subroutine setNumberOfCenters

    !--------------------
    subroutine setEpart(b)
        !--------------------
         logical,intent(in) :: b
        dEpart = b
    end subroutine setEpart

  !---------------------------------------------!
  subroutine setDescendentWeightingGenerations(n)
  !---------------------------------------------!
    integer, intent(in) :: n
    call assert(dRWCount==0,"setDescendentWeightingGenerations: only before RW creation")
    dFathers = n
  end subroutine setDescendentWeightingGenerations

  pure integer function nElecs()
    nElecs = dNElecs
  end function nElecs

  pure integer function nDWGen()
    nDWGen = dFathers
  end function nDWGen

  pure integer function sizeOfRandomWalker()
    if(dEpart)then
      sizeOfRandomWalker = 6*dNElecs + 6 + dFathers+dNElecs+dNElecs**2+dNElecs*dNCenters
    else
      sizeOfRandomWalker = 6*dNElecs + 6 + dFathers
    endif

  end function sizeOfRandomWalker


!-------- set methods -------------------------------


  !-------------------!
  subroutine rw_new(rw)
  !-------------------!
    ! allocate arrays only, x,y,z, eloc and phi are set to zero!

    type(RandomWalker), intent(inout) :: rw
    integer alstat,n,m

    call assert(dNElecs>0,"RandomWalker: rw_new: number of electrons not set")
    ! call assert(dFathers>0,"RandomWalker: rw_new: number of fathers not set")

    n = dNElecs
    m = dNcenters

    rw%w = 1
    rw%phi = 0
    rw%u = 0
    rw%eLocal = 0
    rw%ecp_nonlocal = 0
    !rw%ecp_nonlocalneg= 0
    rw%vpot = 0
    rw%persist = 0

    allocate(rw%x(n),rw%y(n),rw%z(n),stat=alstat)
    call assert(alstat==0,"RandomWalker: rw_new: allocate 1 failed")
    rw%x = 0; rw%y = 0; rw%z = 0
    allocate(rw%xdrift(n),rw%ydrift(n),rw%zdrift(n),stat=alstat)
    call assert(alstat==0,"RandomWalker: rw_new: allocate 2 failed")
    rw%xdrift = 0; rw%ydrift = 0; rw%zdrift = 0
    if(dFathers>0) then
       allocate(rw%fathers(dFathers),stat=alstat)
       call assert(alstat==0,"RandomWalker: rw_new: allocate 3 failed")
       rw%fathers = 0
    endif
    if(dEpart) then
      allocate(rw%Ekin_epart(n),rw%Vne_epart(m,n),rw%Vee_epart(n,n),stat=alstat)
      call assert(alstat==0,"RandomWalker: rw_new: allocate 4 failed")
      rw%Ekin_epart = 0; rw%Vne_epart=0; rw%Vee_epart=0
    endif
    dRWCount = dRWCount + 1
  end subroutine rw_new


  !-----------------------!
  subroutine rw_destroy(rw)
  !-----------------------!
    ! Deallocate all allocated arrays within rw
    type(RandomWalker), intent(inout) :: rw
    integer alstat

    deallocate(rw%x,rw%y,rw%z, stat=alstat)
    call assert(alstat==0,"RandomWalker: rw_destroy: allocate 1 failed")
    deallocate(rw%xdrift,rw%ydrift,rw%zdrift, stat=alstat)
    call assert(alstat==0,"RandomWalker: rw_destroy: allocate 2 failed")
    if (associated(rw%fathers)) deallocate(rw%fathers, stat=alstat)
    call assert(alstat==0,"RandomWalker: rw_destroy: allocate 3 failed")
    if (associated(rw%Ekin_epart)) deallocate(rw%Ekin_epart,rw%Vne_epart,rw%Vee_epart, stat=alstat)
    call assert(alstat==0,"RandomWalker: rw_destroy: allocate 4 failed")

    dRWCount = dRWCount - 1
  end subroutine rw_destroy

  !---------------------------------!
  pure logical function rw_isAllocated(rw)
  !---------------------------------!
     type(RandomWalker), intent(in) :: rw
     rw_isAllocated = associated(rw%x)
  end function rw_isAllocated


  !---------------------------!
  elemental subroutine rw_assign(rw2,rw1)
  !---------------------------!

    type(RandomWalker),intent(inout) :: rw2
    type(RandomWalker),intent(in)    :: rw1
    integer alstat,n,f,m

    n = dNElecs
    m = dNCenters
    f = dFathers

    if (.not. rw_isAllocated(rw1)) then
      ! uninitialized walker, we can't copy anything
      return
    endif

    if (.not. associated(rw2%x)) then
       allocate(rw2%x(n),rw2%y(n),rw2%z(n),rw2%xdrift(n),rw2%ydrift(n),rw2%zdrift(n),stat=alstat)
       if(dFathers>0) allocate(rw2%fathers(f),stat=alstat)
       if(dEpart) then
           allocate(rw2%Ekin_epart(n),rw2%Vne_epart(m,n),rw2%Vee_epart(n,n),stat=alstat)
           rw2%Ekin_epart = 0; rw2%Vne_epart=0; rw2%Vee_epart=0
       endif
    endif

    rw2%x = rw1%x
    rw2%y = rw1%y
    rw2%z = rw1%z
    rw2%w = rw1%w
    rw2%phi     = rw1%phi
    rw2%u       = rw1%u
    rw2%eLocal  = rw1%eLocal
    !rw2%ecp_nonlocalneg  = rw1%ecp_nonlocalneg
    rw2%ecp_nonlocal  = rw1%ecp_nonlocal
    rw2%vPot    = rw1%vPot
    rw2%xdrift  = rw1%xdrift
    rw2%ydrift  = rw1%ydrift
    rw2%zdrift  = rw1%zdrift
    rw2%persist = rw1%persist
    if(dFathers>0) rw2%fathers = rw1%fathers
    if(dEpart)then
      rw2%Ekin_epart = rw1%Ekin_epart
      rw2%Vne_epart  = rw1%Vne_epart
      rw2%Vee_epart  = rw1%Vee_epart
    endif

  end subroutine rw_assign


  !-----------------------------------------!
  subroutine writeRandomWalkerToStream(rw,iu)
  !-----------------------------------------!
    ! write random walker to open output stream
    type(RandomWalker),intent(in)    :: rw
    integer, intent(in)              :: iu

    write(iu) rw%x,rw%y,rw%z
    write(iu) rw%w,rw%phi,rw%u,rw%eLocal,rw%vPot
    write(iu) rw%xdrift,rw%ydrift,rw%zdrift
    write(iu) rw%persist
    if(dFathers>0) write(iu) rw%fathers
    if(dEpart)then
      write(iu) rw%Ekin_epart
      write(iu) rw%Vne_epart
      write(iu) rw%Vee_epart
    endif

  end subroutine writeRandomWalkerToStream


  !------------------------------------------!
  subroutine readRandomWalkerFromStream(rw,iu)
  !------------------------------------------!
    ! reads random walker to open output stream, overwriting current values
    ! allocating arrays if necessary
    type(RandomWalker),intent(inout) :: rw
    integer, intent(in)              :: iu
    integer n,alstat,m

    if (.not.associated(rw%x)) then
       n = dNElecs
       m = dNCenters
       allocate(rw%x(n),rw%y(n),rw%z(n),stat=alstat)
        call assert(alstat==0,"RandomWalker: readRandomWalkerFromStream allocate 1 failed")
        allocate(rw%xdrift(n),rw%ydrift(n),rw%zdrift(n),stat=alstat)
        call assert(alstat==0,"RandomWalker: readRandomWalkerFromStream allocate 2 failed")
        if(dFathers>0) then
         allocate(rw%fathers(dFathers),stat=alstat)
        call assert(alstat==0,"RandomWalker: readRandomWalkerFromStream allocate 3 failed")
        endif
       if(dEpart) then
         allocate(rw%Ekin_epart(n),rw%Vne_epart(m,n),rw%Vee_epart(n,n),stat=alstat)
         call assert(alstat==0,"RandomWalker: readRandomWalkerFromStream allocate 4 failed")
         rw%Ekin_epart = 0; rw%Vne_epart=0; rw%Vee_epart=0
       endif
    endif

    read(iu) rw%x,rw%y,rw%z
    read(iu) rw%w,rw%phi,rw%u,rw%eLocal,rw%vPot
    read(iu) rw%xdrift,rw%ydrift,rw%zdrift
    read(iu) rw%persist
    if(dFathers>0) read(iu) rw%fathers
    if(dEpart)then
      read(iu) rw%Ekin_epart
      read(iu) rw%Vne_epart
      read(iu) rw%Vee_epart
    endif

  end subroutine readRandomWalkerFromStream


  !----------------------------------------!
  subroutine writeRandomWalkerToVector(rw,v)
  !----------------------------------------!

    type(RandomWalker),intent(in)      :: rw
    real(r8), intent(inout)              :: v(:)
    integer n,i

    call assert(size(v)==sizeOfRandomWalker(),"RandomWalker: vector has wrong size")
    n=dNElecs
    i=1;   v(i:i+n-1)=rw%x
    i=i+n; v(i:i+n-1)=rw%y
    i=i+n; v(i:i+n-1)=rw%z
    i=i+n; v(i)=rw%w; v(i+1)=rw%phi; v(i+2)=rw%u; v(i+3)=rw%eLocal; v(i+4)=rw%vPot
    i=i+5; v(i:i+n-1)=rw%xdrift
    i=i+n; v(i:i+n-1)=rw%ydrift
    i=i+n; v(i:i+n-1)=rw%zdrift
    i=i+n; v(i)=rw%persist
    if (dFathers > 0) then
       i = i + 1; v(i:i+dFathers-1) = rw%fathers
    end if
    if (dEpart) then
      i=i+dfathers; v(i:i+n-1) = rw%Ekin_epart
      i=i+n; v(i:i+dEE-1) = reshape(rw%Vee_epart,(/dEE/))
      i=i+dEE;v(i:i+dNE-1) = reshape(rw%Vne_epart,(/dNE/))
    endif

  end subroutine writeRandomWalkerToVector


  !-----------------------------------------!
  subroutine readRandomWalkerFromVector(rw,v)
  !-----------------------------------------!

    type(RandomWalker),intent(inout)   :: rw
    real(r8), intent(inout)              :: v(:)
    integer n,i,alstat,m

    call assert(size(v)==sizeOfRandomWalker(),"RandomWalker: vector has wrong size")

    n=dNElecs
    m=dNCenters
    if (.not.associated(rw%x)) then
       allocate(rw%x(n),rw%y(n),rw%z(n),stat=alstat)
       if (alstat /= 0) call abortp("RandomWalker:fromVector: allocate 1 failed")
       allocate(rw%xdrift(n),rw%ydrift(n),rw%zdrift(n),stat=alstat)
       if (alstat /= 0) call abortp("RandomWalker:fromVector: allocate 2 failed")
       allocate(rw%fathers(dFathers),stat=alstat)
       if (alstat /= 0) call abortp("RandomWalker:fromVector: allocate 3 failed")
       if(dFathers>0) then
            allocate(rw%fathers(dFathers),stat=alstat)
            call assert(alstat==0,"RandomWalker: readRandomWalkerFromVectora allocate 4 failed")
        endif
       if(dEpart) then
         allocate(rw%Ekin_epart(n),rw%Vne_epart(m,n),rw%Vee_epart(n,n),stat=alstat)
         call assert(alstat==0,"RandomWalker: readRandomWalkerFromVector allocate 5 failed")
         rw%Ekin_epart = 0; rw%Vne_epart=0; rw%Vee_epart=0
       endif

    endif

    i=1;   rw%x=v(i:i+n-1)
    i=i+n; rw%y=v(i:i+n-1)
    i=i+n; rw%z=v(i:i+n-1)
    i=i+n; rw%w=v(i); rw%phi=v(i+1); rw%u=v(i+2); rw%eLocal=v(i+3); rw%vPot=v(i+4)
    i=i+5; rw%xdrift=v(i:i+n-1)
    i=i+n; rw%ydrift=v(i:i+n-1)
    i=i+n; rw%zdrift=v(i:i+n-1)
    i=i+n; rw%persist=nint(v(i))
    i=i+1
    if(dFathers>0) rw%fathers=INT(v(i:i+dFathers-1))
    if(dEpart)then
       i=i+dfathers; rw%Ekin_epart = v(i:i+n-1)
       i=i+n; rw%Vee_epart = reshape(v(i:i+dEE-1),(/dNElecs,dNElecs/))
       i=i+dEE; rw%Vne_epart = reshape(v(i:i+dNE-1),(/dNCenters,dNElecs/))
    endif
  end subroutine readRandomWalkerFromVector


  !----------------------------------------!
  subroutine sendRandomWalkerTo(rw,node,tag)
  !----------------------------------------!

    type(RandomWalker), intent(in)      :: rw
    integer, intent(in)                 :: node
    integer, intent(in)                 :: tag
    real(r8)                              :: v(sizeOfRandomWalker())
    integer n

    n = sizeOfRandomWalker()
    call writeRandomWalkerToVector(rw,v)
    call myMPISendDouble(v,n,node,tag)

  end subroutine sendRandomWalkerTo


  !---------------------------------------------!
  subroutine receiveRandomWalkerFrom(rw,node,tag)
  !---------------------------------------------!

    type(RandomWalker), intent(inout)   :: rw
    integer, intent(in)                 :: node
    integer, intent(in)                 :: tag
    real(r8)                              :: v(sizeOfRandomWalker())
    integer n

    n = sizeOfRandomWalker()
    call myMPIReceiveDouble(v,n,node,tag)
    call readRandomWalkerFromVector(rw,v)

  end subroutine receiveRandomWalkerFrom


  !----------------------------!
  subroutine recalculateEloc(rw)
  !----------------------------!

    ! call eloc and updates components. This is necessary when the wave function
    ! has been modified (as in jastrow optimization)

    type(RandomWalker), intent(inout) :: rw
    type(eConfigArray)  :: ec

    call eConfigArray_new(ec,dNElecs,1)
    call eConfigArray_set(ec,1,rw%x,rw%y,rw%z)
    call eloc(0,ec,'none')

    rw%phi     = elPhi(1)
    rw%u       = elU(1)
    rw%eLocal  = elEloc(1)
    rw%ecp_nonlocal = elECPPotNl
    rw%vPot    = elVee(1) + elVen(1) + vpot0         ! with nucl-nucl repulsion
    rw%vee     = elVee(1)
    rw%xdrift  = elxDrift(1:dNElecs,1)
    rw%ydrift  = elyDrift(1:dNElecs,1)
    rw%zdrift  = elzDrift(1:dNElecs,1)
    if(dEpart)then
      rw%Ekin_epart =  elEkin_epart(1:dNElecs,1)
      rw%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,1)
      rw%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,1)
    endif
    call eConfigArray_destroy(ec)
  end subroutine recalculateEloc


  !-------------------------------------!
  subroutine rw_recalculateElocBlock(rwb,twoLevelStep,doCalc)
  !-------------------------------------!

    ! This version recalculates a whole block of walker
    ! call eloc and updates components. This is necessary when the wave function
    ! has been modified (as in jastrow optimization)

    type(RandomWalker), intent(inout) :: rwb(:)
    integer, intent(in), optional :: twoLevelStep
    logical, intent(in), optional :: doCalc(:) ! for twoLevel
    integer w
    type(eConfigArray)  :: ec

    call eConfigArray_new(ec,dNElecs,size(rwb))

    do w=1,size(rwb)
       call eConfigArray_set(ec,w,rwb(w)%x,rwb(w)%y,rwb(w)%z)
    enddo
    call eloc(0,ec,'none',twoLevelStep=twoLevelStep,doCalc=doCalc)

    do w=1,size(rwb)
       if(present(doCalc)) then
        if(.not. doCalc(w)) cycle
       endif
       rwb(w)%phi     = elPhi(w)
       rwb(w)%u       = elU(w)
       rwb(w)%eLocal  = elEloc(w)
       rwb(w)%ecp_nonlocal = elECPPotNl
       rwb(w)%vPot    = elVen(w) + elVee(w) + vpot0         ! with nucl-nucl repulsion
       rwb(w)%vee     = elVee(w)
       rwb(w)%xdrift  = elxDrift(1:dNElecs,w)
       rwb(w)%ydrift  = elyDrift(1:dNElecs,w)
       rwb(w)%zdrift  = elzDrift(1:dNElecs,w)
       if(dEpart)then
          rwb(w)%Ekin_epart =  elEkin_epart(1:dNElecs,w)
          rwb(w)%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,w)
          rwb(w)%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,w)
       endif
    enddo

    call eConfigArray_destroy(ec)
  end subroutine rw_recalculateElocBlock


  !----------------------------------!
  subroutine resetTon(rw,x,y,z,optType)
  !----------------------------------!

    ! sets position vector and calculates eloc data member, resets the others

    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: x(:),y(:),z(:)  ! position
    character(len=9), optional        :: optType
    type(eConfigArray)  :: ec
    call assert(size(x)>=dNElecs,"RandomWalker:resetTon: illegal dimension")

    call eConfigArray_new(ec,dNElecs,1)
    call eConfigArray_set(ec,1,x,y,z)
    ! We only need to calculate vNlk if we have ECP! and we do optimization.
    !!if (present(optType)) then
    !!  if(use_ecp) then
    !!    call eloc(0,ec,optType,CalcvNlk=.true.)
    !!  else
    !!   call eloc(0,ec,optType)
    !! endif
    !!else
       call eloc(0,ec,'none')
    !!endif

    rw%x       = x(1:dNElecs)
    rw%y       = y(1:dNElecs)
    rw%z       = z(1:dNElecs)
    rw%w       = 1
    rw%phi     = elPhi(1)
    rw%u       = elU(1)
    rw%eLocal  = elEloc(1)
    rw%ecp_nonlocal = elECPPotNl
    rw%vPot    = elVee(1) + elVen(1) + vpot0         ! with nucl-nucl repulsion
    rw%vee     = elVee(1)
    rw%xdrift  = elxDrift(1:dNElecs,1)
    rw%ydrift  = elyDrift(1:dNElecs,1)
    rw%zdrift  = elzDrift(1:dNElecs,1)
    rw%persist = 0
    if (dFathers > 0) rw%fathers = 0
    if(dEpart)then
      rw%Ekin_epart =  elEkin_epart(1:dNElecs,1)
      rw%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,1)
      rw%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,1)
    endif
    call eConfigArray_destroy(ec)
  end subroutine resetTon

   !----------------------------------!
   subroutine resetTo3n(rw,r,optType)
   !----------------------------------!
      type(RandomWalker), intent(inout) :: rw
      real(r8), intent(in) :: r(:)
      character(len=9), optional :: optType
      real(r8), allocatable :: x(:), y(:), z(:)
      integer :: i, len

      call assert(MOD(SIZE(r),3)==0,"RandomWalker:resetTo3n: illegal dimension")

      len = SIZE(r)/3
      allocate(x(len), y(len), z(len))

      do i = 1, len
         x(i) = r(3*i-2)
         y(i) = r(3*i-1)
         z(i) = r(3*i)
      end do

      call resetTon(rw,x,y,z,optType)
   end subroutine resetTo3n

  !----------------------------------!
  subroutine resetTo_without_Calc(rw,x,y,z)
  !----------------------------------!

    ! sets position vector and eloc data member, resets the others without call to eloc
    ! can be used if we alrady called eloc and just need an update for rw

    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: x(:),y(:),z(:)  ! position

    call assert(size(x)>=dNElecs,"RandomWalker:resetTo: illegal dimension")


    rw%x       = x(1:dNElecs)
    rw%y       = y(1:dNElecs)
    rw%z       = z(1:dNElecs)
    rw%w       = 1
    rw%phi     = elPhi(1)
    rw%u       = elU(1)
    rw%eLocal  = elEloc(1)
    rw%ecp_nonlocal = elECPPotNl
    rw%vPot    = elVee(1) + elVen(1) + vpot0         ! with nucl-nucl repulsion
    rw%vee     = elVee(1)
    rw%xdrift  = elxDrift(1:dNElecs,1)
    rw%ydrift  = elyDrift(1:dNElecs,1)
    rw%zdrift  = elzDrift(1:dNElecs,1)
    rw%persist = 0
    if (dFathers > 0) rw%fathers = 0
    if(dEpart)then
      rw%Ekin_epart =  elEkin_epart(1:dNElecs,1)
      rw%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,1)
      rw%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,1)
    endif

  end subroutine resetTo_without_Calc

  !----------------------------------------!
  subroutine rw_resetToBlock(rwb,ec)
  !----------------------------------------!

    ! sets position vector and calculates eloc data member, resets the others

    type(RandomWalker), intent(inout) :: rwb(:)
    type(eConfigArray), intent(inout) :: ec
    integer w
    call assert(size(rwb)==eConfigArray_size(ec),'rw_resetToBlock: sizes do not match')

    call eloc(0,ec,'none')

    do w=1,size(rwb)
       call eConfigArray_get(ec,w,rwb(w)%x,rwb(w)%y,rwb(w)%z)
       rwb(w)%w       = 1
       rwb(w)%phi     = elPhi(w)
       rwb(w)%u       = elU(w)
       rwb(w)%eLocal  = elEloc(w)
       rwb(w)%ecp_nonlocal = elECPPotNl
       rwb(w)%vPot    = elVee(w) + elVen(w) + vpot0         ! with nucl-nucl repulsion
       rwb(w)%vee     = elVee(w)
       rwb(w)%xdrift  = elxDrift(:,w)
       rwb(w)%ydrift  = elyDrift(:,w)
       rwb(w)%zdrift  = elzDrift(:,w)
       rwb(w)%persist = 0
       if (dFathers > 0) rwb(w)%fathers = 0
       if(dEpart)then
         rwb(w)%Ekin_epart =  elEkin_epart(1:dNElecs,w)
         rwb(w)%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,w)
         rwb(w)%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,w)
       endif
    enddo
  end subroutine rw_resetToBlock


  !------------------------!
  subroutine setTo(rw,x,y,z)
  !------------------------!

    ! sets position vector and calculates eloc data member

    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: x(:),y(:),z(:)  ! position
    type(eConfigArray)  :: ec
    call assert(size(x)>=dNElecs,"RandomWalker:setTo: illegal dimension")

    call eConfigArray_new(ec,dNElecs,1)

    call eConfigArray_set(ec,1,x,y,z)

    call eloc(0,ec,'none')
    rw%x       = x(1:dNElecs)
    rw%y       = y(1:dNElecs)
    rw%z       = z(1:dNElecs)
    rw%phi     = elPhi(1)
    rw%u       = elU(1)
    rw%eLocal  = elEloc(1)
    rw%ecp_nonlocal = elECPPotNl
    rw%vPot    = elVee(1) + elVen(1) + vpot0         ! with nucl-nucl repulsion
    rw%vee     = elVee(1)
    rw%xdrift  = elxDrift(1:dNElecs,1)
    rw%ydrift  = elyDrift(1:dNElecs,1)
    rw%zdrift  = elzDrift(1:dNElecs,1)
    if(dEpart)then
      rw%Ekin_epart =  elEkin_epart(1:dNElecs,1)
      rw%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,1)
      rw%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,1)
    endif

    call eConfigArray_destroy(ec)
  end subroutine setTo



  subroutine rw_setToBlock(rwb, ec, twoLevelStep, doCalc, tmove, tau, isMoved)
  !--------------------------------------------------------------------------!

    ! sets position vector and calculates eloc data member at that position

    type(RandomWalker), intent(inout) :: rwb(:)
    type(eConfigArray), intent(inout) :: ec      ! when T move accepted: new positions
    integer, intent(in), optional :: twoLevelStep
    logical, intent(in), optional :: doCalc(eConfigArray_size(ec))
    integer, intent(in), optional :: tmove       ! triggers T move calculation
    real(r8), intent(in), optional  :: tau         ! time step for T moves
    logical, intent(inout), optional :: isMoved(:)  ! T moves accepted
    integer w
    call assert(size(rwb)==eConfigArray_size(ec),'rw_resetToBlock: sizes do not match')

    if (present(twoLevelStep)) then
       call eloc(0,ec,'none',twoLevelStep=twoLevelStep,doCalc=doCalc)

       do w = 1, size(rwb)
          if (present(doCalc)) then
             if (.not. doCalc(w)) cycle
          endif
          call eConfigArray_get(ec,w,rwb(w)%x,rwb(w)%y,rwb(w)%z)
          rwb(w)%phi     = elPhi(w)
          rwb(w)%u       = elU(w)
          rwb(w)%eLocal  = elEloc(w)
          rwb(w)%ecp_nonlocal=elECPPotNl
          rwb(w)%vPot    = elVee(w) + elVen(w) + vpot0         ! with nucl-nucl repulsion
          rwb(w)%vee     = elVee(w)
          rwb(w)%xdrift  = elxDrift(:,w)
          rwb(w)%ydrift  = elyDrift(:,w)
          rwb(w)%zdrift  = elzDrift(:,w)
          if(dEpart)then
             rwb(w)%Ekin_epart =  elEkin_epart(1:dNElecs,w)
             rwb(w)%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,w)
             rwb(w)%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,w)
          endif
       enddo

    else
       ! first set positions (ec might change in T moves)
       do w = 1, size(rwb)
          call eConfigArray_get(ec,w,rwb(w)%x,rwb(w)%y,rwb(w)%z)
       end do
       if (present(tmove)) then
          call eloc(0, ec, 'none', tmove=tmove, tau=tau, isMoved=isMoved)
       else
          call eloc(0, ec, 'none')
       end if
       ! then set results
       do w = 1, size(rwb)
          rwb(w)%phi     = elPhi(w)
          rwb(w)%u       = elU(w)
          rwb(w)%eLocal  = elEloc(w)
          rwb(w)%ecp_nonlocal = elECPPotNl
          rwb(w)%vPot    = elVee(w) + elVen(w) + vpot0         ! with nucl-nucl repulsion
          rwb(w)%vee     = elVee(w)
          rwb(w)%xdrift  = elxDrift(:,w)
          rwb(w)%ydrift  = elyDrift(:,w)
          rwb(w)%zdrift  = elzDrift(:,w)
          if (dEpart) then
             rwb(w)%Ekin_epart =  elEkin_epart(1:dNElecs,w)
             rwb(w)%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,w)
             rwb(w)%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,w)
          end if
       end do
    end if
  end subroutine rw_setToBlock


  !-----------------------------!
  subroutine setOneTo(rw,n,x,y,z)
  !-----------------------------!

    ! sets position of n-th electron and calculates eloc data members
    ! NOTE: CALL eloc FOR ALL ELECTRONS

    type(RandomWalker), intent(inout) :: rw
    integer, intent(in)               :: n
    real(r8), intent(in)                :: x,y,z  ! position n-th electron
    type(eConfigArray)                :: ec

    call assert(n>0 .and. n<=dNElecs,"RandomWalker:setOneTo: wrong electron")

    rw%x(n)    = x
    rw%y(n)    = y
    rw%z(n)    = z
    call eConfigArray_new(ec,dNElecs,1)
    call eConfigArray_set(ec,1,rw%x,rw%y,rw%z)
    call eloc(0,ec,'none')

    rw%phi     = elPhi(1)
    rw%u       = elU(1)
    rw%eLocal  = elEloc(1)
    rw%ecp_nonlocal = elECPPotNl
    rw%vPot    = elVee(1) + elVen(1) + vpot0         ! with nucl-nucl repulsion
    rw%vee     = elVee(1)
    rw%xdrift  = elxDrift(1:dNElecs,1)
    rw%ydrift  = elyDrift(1:dNElecs,1)
    rw%zdrift  = elzDrift(1:dNElecs,1)
    if(dEpart)then
      rw%Ekin_epart =  elEkin_epart(1:dNElecs,1)
      rw%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,1)
      rw%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,1)
    endif
    call eConfigArray_destroy(ec)
  end subroutine setOneTo


  !--------------------------!
  subroutine moveAll(rw,x,y,z)
  !--------------------------!

    ! moves position vector by x,y,z and calculates eloc data member

    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: x(:),y(:),z(:)  ! position
    type(eConfigArray)                :: ec

    call assert(size(x)>=dNElecs,"RandomWalker:recalculateElec: wrong dimension")
    rw%x       = rw%x + x(1:dNElecs)
    rw%y       = rw%y + y(1:dNElecs)
    rw%z       = rw%z + z(1:dNElecs)

    call eConfigArray_new(ec,dNElecs,1)
    call eConfigArray_set(ec,1,rw%x,rw%y,rw%z)
    call eloc(0,ec,'none')

    rw%phi     = elPhi(1)
    rw%u       = elU(1)
    rw%eLocal  = elEloc(1)
    rw%ecp_nonlocal = elECPPotNl
    rw%vPot    = elVee(1) + elVen(1) + vpot0         ! with nucl-nucl repulsion
    rw%vee     = elVee(1)
    rw%xdrift  = elxDrift(1:dNElecs,1)
    rw%ydrift  = elyDrift(1:dNElecs,1)
    rw%zdrift  = elzDrift(1:dNElecs,1)
    if(dEpart)then
      rw%Ekin_epart =  elEkin_epart(1:dNElecs,1)
      rw%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,1)
      rw%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,1)
    endif

    call eConfigArray_destroy(ec)
  end subroutine moveAll


  !----------------------------!
  subroutine moveOne(rw,n,x,y,z)
  !----------------------------!

    ! sets position of n-th electron and calculates eloc data member

    type(RandomWalker), intent(inout) :: rw
    integer, intent(in)               :: n
    real(r8), intent(in)                :: x,y,z  ! position n-th electron
    type(eConfigArray)                :: ec

    call assert(n>0 .and. n<=dNElecs,"RandomWalker:setOneTo: wrong electron")

    rw%x(n)    = rw%x(n) + x
    rw%y(n)    = rw%y(n) + y
    rw%z(n)    = rw%z(n) + z

    call eConfigArray_new(ec,dNElecs,1)
    call eConfigArray_set(ec,1,rw%x,rw%y,rw%z)
    call eloc(n,ec,'none')

    rw%phi     = elPhi(1)
    rw%u       = elU(1)
    rw%eLocal  = elEloc(1)
    rw%ecp_nonlocal= elECPPotNl
    rw%vPot    = elVee(1) + elVen(1) + vpot0         ! with nucl-nucl repulsion
    rw%vee     = elVee(1)
    rw%xdrift  = elxDrift(1:dNElecs,1)
    rw%ydrift  = elyDrift(1:dNElecs,1)
    rw%zdrift  = elzDrift(1:dNElecs,1)
    if(dEpart)then
      rw%Ekin_epart =  elEkin_epart(1:dNElecs,1)
      rw%Vne_epart  =  elVne_epart(1:dNcenters,1:dNElecs,1)
      rw%Vee_epart  =  elVee_epart(1:dNElecs,1:dNElecs,1)
    endif

    call eConfigArray_destroy(ec)
  end subroutine moveOne

  !-----------------------------!
  subroutine setfather(rw,i,n)
  !-----------------------------!
    type(RandomWalker), intent(inout) :: rw
    integer, intent(in)               :: n,i

    call assert(i>0.and.i<=dFathers,"RandomWalker:setFather: wrong number")
    rw%fathers(i) = n

  end subroutine setfather

  !------------------------!
  subroutine setWeight(rw,w)
  !------------------------!
    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: w

    rw%w = w
  end subroutine setWeight

  !------------------------!
  subroutine eLocalSave(rw,el) ! only for Tmove
  !------------------------!
    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: el
    rw%eLocalSave = el
  end subroutine eLocalSave

  !------------------------!
  subroutine eLocalSet(rw,el,tmoved) ! only for Tmove
  !------------------------!
    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: el
    logical,intent(in)                :: tmoved
    rw%eLocal = el
    rw%tmovedrw = tmoved
  end subroutine eLocalSet

  !------------------------!
  subroutine resetTrw(rw) ! only for Tmove
  !------------------------!
    type(RandomWalker), intent(inout) :: rw

    rw%tmovedrw = .false.
  end subroutine resetTrw

  !-----------------------------!
  subroutine multiplyWeight(rw,w)
  !-----------------------------!
    type(RandomWalker), intent(inout) :: rw
    real(r8), intent(in)                :: w

    rw%w = rw%w * w
  end subroutine multiplyWeight


  !-----------------------------!
  subroutine resetPersistency(rw)
  !-----------------------------!
    type(RandomWalker), intent(inout) :: rw

    rw%persist = 0
  end subroutine resetPersistency

  !-----------------------------!
  subroutine setPersist(rw,p)
  !-----------------------------!
    type(RandomWalker), intent(inout) :: rw
    integer                           :: p
    rw%persist = p
  end subroutine setPersist

  !----------------------------!
  subroutine incrPersistency(rw)
  !----------------------------!
    type(RandomWalker), intent(inout) :: rw

    rw%persist = rw%persist + 1
  end subroutine incrPersistency



!------------ access functions -----------------------

  pure real(r8) function wgt(rw)
    type(RandomWalker), intent(in) :: rw
    wgt = rw%w
  end function wgt

  pure real(r8) function psi(rw)
    type(RandomWalker), intent(in) :: rw
    psi = rw%phi * exp(rw%u)
  end function psi

  pure real(r8) function lnpsi2(rw)
    type(RandomWalker), intent(in) :: rw
    lnpsi2 = 2*( log(abs(rw%phi)) + rw%u)
  end function lnpsi2

  pure real(r8) function phi(rw)
    type(RandomWalker), intent(in) :: rw
    phi = rw%phi
  end function phi

  pure real(r8) function ju(rw)
    type(RandomWalker), intent(in) :: rw
    ju = rw%u
  end function ju

  pure real(r8) function E_local(rw)            ! local energy H Psi / Psi
    type(RandomWalker), intent(in) :: rw
    E_local = rw%eLocal
  end function E_local

  pure real(r8) function E_localSaveGet(rw)            ! For Tmove
    type(RandomWalker), intent(in) :: rw
    E_localSaveGet = rw%eLocalSave
  end function E_localSaveGet

  pure logical function isTmoved(rw)                 !  For Tmove
    type(RandomWalker), intent(in) :: rw
    isTmoved = rw%tmovedrw
  end function isTmoved


    pure real(r8) function ECP_nl(rw)            ! local energy H Psi / Psi
    type(RandomWalker), intent(in) :: rw
    ECP_nl = rw%ecp_nonlocal
  end function ECP_nl

  pure real(r8) function E_pot(rw)              ! local potential energy
    type(RandomWalker), intent(in) :: rw
    E_pot = rw%vPot
  end function E_pot

  pure real(r8) function E_ee(rw)              ! local potential energy
    type(RandomWalker), intent(in) :: rw
    E_ee = rw%vee
  end function E_ee

  pure real(r8) function E_kin(rw)              ! local kinetic energy
    type(RandomWalker), intent(in) :: rw
    E_kin = rw%eLocal - rw%vPot
  end function E_kin

  subroutine getNablaPsi(rw,NablaPsi)
    type(RandomWalker), intent(in) :: rw
    real(r8), intent(out)            :: NablaPsi(:)
    integer :: i
    call assert(size(NablaPsi)==3*dNElecs,"RandomWalker:getNablaPsi wrong dimension")
    do i=1,dNElecs
       NablaPsi(3*i-2) = rw%xdrift(i) * psi(rw)
       NablaPsi(3*i-1) = rw%ydrift(i) * psi(rw)
       NablaPsi(3*i)   = rw%zdrift(i) * psi(rw)
    end do
  end subroutine getNablaPsi

  pure real(r8) function LaplPsi(rw)            ! Laplacian(Psi)
    type(RandomWalker), intent(in) :: rw
    LaplPsi = - 2.d0 * E_kin(rw) * psi(rw)
  end function LaplPsi

  pure real(r8) function HPsi(rw)               ! H(Psi)
    type(RandomWalker), intent(in) :: rw
    HPsi = rw%eLocal * psi(rw)
  end function HPsi

  pure integer function persist(rw)
    type(RandomWalker), intent(in) :: rw
    persist = rw%persist
  end function persist

  pure integer function father(rw,n)
    type(RandomWalker), intent(in) :: rw
    integer, intent(in)            :: n
    father = rw%fathers(n)
  end function father


  pure subroutine posn(rw,x,y,z)                     ! position vector
    type(RandomWalker), intent(in) :: rw
    real(r8), intent(out)            ::  x(:),y(:),z(:)
    x = rw%x; y = rw%y; z = rw%z
  end subroutine posn

  pure subroutine pos3n(rw,r)                     ! 3n position vector
     type(RandomWalker), intent(in) :: rw
     real(r8), intent(out) :: r(:)
     integer :: i

     !call assert(SIZE(r) == 3 * dNElecs, "pos3n: illegal size of r")

     do i = 1, dNElecs
        r(3*i-2) = rw%x(i)
        r(3*i-1) = rw%y(i)
        r(3*i)   = rw%z(i)
     end do
  end subroutine pos3n

  ! note: this should be "const pointer" which is impossible in fortran
  ! possibly breaks consistency
  subroutine posPtr(rw,x,y,z)                       ! position vector
    type(RandomWalker), intent(inout) :: rw
    real(r8), pointer ::  x(:),y(:),z(:)
    x => rw%x; y => rw%y; z => rw%z
  end subroutine posPtr


  pure subroutine drift(rw,xdrift,ydrift,zdrift)    ! drift vector
    type(RandomWalker), intent(in) :: rw
    real(r8), intent(out)            :: xdrift(:),ydrift(:),zdrift(:)
    xdrift = rw%xdrift; ydrift = rw%ydrift; zdrift = rw%zdrift
  end subroutine drift

  ! note: this should be "const pointer" which is impossible in fortran
  ! possibly breaks consistency
  subroutine driftPtr(rw,xdrift,ydrift,zdrift)      ! drift vector
    type(RandomWalker), intent(inout) :: rw
    real(r8), pointer ::  xdrift(:),ydrift(:),zdrift(:)
    xdrift => rw%xdrift; ydrift => rw%ydrift; zdrift => rw%zdrift
  end subroutine driftPtr

  ! ENERGY PARTIONING
  pure subroutine Ekini(rw,kini)
    type(RandomWalker),intent(in) :: rw
    real(r8),intent(out)            :: kini(:)
    kini = rw%Ekin_epart
  end subroutine Ekini

  pure subroutine EVne(rw,vnei)
    type(RandomWalker),intent(in) :: rw
    real(r8),intent(out)            :: vnei(:,:)
    vnei = rw%Vne_epart
  end subroutine EVne

  pure subroutine EVee(rw,veei)
    type(RandomWalker),intent(in) :: rw
    real(r8),intent(out)            :: veei(:,:)
    veei = rw%Vee_epart
  end subroutine EVee

  subroutine display(rw,mode,iu)
    type(RandomWalker), intent(in) :: rw
    integer, intent(in)            :: mode
    integer, intent(in), optional  :: iu

    if (present(iu)) then
       write(iu,*) "displaying random walker:"
       write(iu,'(3(A,G15.8))') "psi=",psi(rw)," E_loc=",rw%eLocal," w=",rw%w
       write(iu,'(2(A,G15.8))') "phi=",rw%phi,"u=",rw%u
       write(iu,'(2(A,G15.8))') "vpot=",rw%vPot,"ecp=",elECPPot
       select case (mode)
       case (2)
          write(iu,'(A)') ' x      = '
          write(iu,'(10F15.8)') rw%x
          write(iu,'(A)') ' y      = '
          write(iu,'(10F15.8)') rw%y
          write(iu,'(A)') ' z      = '
          write(iu,'(10F15.8)') rw%z
       end select
    else
       write(*,*) "displaying random walker"
       write(*,'(3(A,G15.8))') "psi=",psi(rw)," E_loc=",rw%eLocal," w=",rw%w
       write(*,'(2(A,G15.8))') "phi=",rw%phi,"u=",rw%u
       write(*,'(2(A,G15.8))') "vpot=",rw%vPot,"ecp=",elECPPot
       select case (mode)
       case (2)
          write(*,'(A)') ' x      = '
          write(*,'(10F15.8)') rw%x
          write(*,'(A)') ' y      = '
          write(*,'(10F15.8)') rw%y
          write(*,'(A)') ' z      = '
          write(*,'(10F15.8)') rw%z
       end select
    endif
  end subroutine display


  subroutine getNearestNucleusDist(rw,Rnnu,Nnu)
  !-------------------------------------------!

    type(RandomWalker), intent(in) :: rw
    real(r8), intent(out)            :: Rnnu(:)
    integer, intent(out)           :: Nnu(:)
    integer                        :: i,a
    real(r8)                         :: rai
    real(r8)                         :: x(ne),y(ne),z(ne)

    call pos(rw,x,y,z)

    do i=1,ne
      Rnnu(i) = 10000.0
      do a=1,ncenter
        rai = sqrt( (x(i)-atoms(a)%cx)**2 + (y(i)-atoms(a)%cy)**2 + (z(i)-atoms(a)%cz)**2 )
        if (rai < Rnnu(i)) then
           Rnnu(i) = rai
           Nnu(i) = a
        endif
      enddo
    enddo

  end subroutine getNearestNucleusDist



END MODULE randomWalker_m
