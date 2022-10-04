! Copyright (C) 1996, 2012-2013, 2015-2016, 2019 Arne Luechow
! Copyright (C) 2013 Alexander Sturm
!
! SPDX-License-Identifier: GPL-3.0-or-later

module initialPositions_m

  use kinds_m, only: r8
  use error_m
  use random_m, only: myran, mygran
  use atom_m
  use wfData_m
  use aos_m, only: ao1calc,ao1splcalc
  use mos_m, only: mo1calc,mat
  use multiDet_m, only: mclist
  use eloc_m, only: wf_getNCoreElecs
  implicit none

  real(r8), allocatable :: mM(:,:,:)           ! matrix for each MO
  real(r8), allocatable :: mMu(:,:)            ! mean for each MO
  real(r8)              :: mMOScale=1.d0       ! scale factor for MO covs
  integer, allocatable:: mMOList(:)          ! MO number is which index?

  private
  public :: init_LMOSampling, destroy_LMOSampling, createRandomElectronPositions, lmoPosition, getLMOSize

contains


  !--------------------------------------------------!
  subroutine createRandomElectronPositions(mode,x,y,z)
  !--------------------------------------------------!

! create random electron positions
! very simple here: just add gaussian random number to the
! coordinate cx,cy,cz, of the atom(=center) for as many electron as
! the nuclear charge za. Improvement: The variance of the gaussian
! should reflect the atomic radius.

    integer, intent(in) :: mode
    real(r8), intent(out) :: x(:)
    real(r8), intent(out) :: y(:)
    real(r8), intent(out) :: z(:)

    call assert(size(x)==ne,'createRandomElectronPositions: array size and ne differ')

    select case(mode)
    case(1)
       call simpleGaussianPositions(x,y,z)
    case(2)
       call densityPositions(x,y,z)
    case(3)
       call lmoPositions(x,y,z)
    case default
       call abortp('createRandomElectronPositions: illegal mode')
    end select
  end subroutine createRandomElectronPositions


  !---------------------------------------!
  subroutine simpleGaussianPositions(x,y,z)
  !---------------------------------------!

! create random walker
! very simple here: just add gaussian random number to the
! coordinate cx,cy,cz, of the atom(=center) for as many electron as
! the nuclear charge za. Improvement: The variance of the gaussian
! should reflect the atomic radius.

    real(r8), intent(out) :: x(:)
    real(r8), intent(out) :: y(:)
    real(r8), intent(out) :: z(:)

    integer :: i,ii,ii1,ii2,a

    ii  = 0                      ! electron counter
    ii1 = 0                      ! electron counter ALPHA
    ii2 = nalpha                 ! electron counter BETA
    do a=1,ncenter
       do i=1,atoms(a)%za - atoms(a)%pa  ! partial charges
          ii = ii + 1            ! next electron
          if (mod(ii,2)==1) then ! alpha electron
             ii1    = ii1 + 1
             x(ii1) = atoms(a)%cx + mygran()
             y(ii1) = atoms(a)%cy + mygran()
             z(ii1) = atoms(a)%cz + mygran()
          else if (ii2>=ne) then  ! alpha electron
             ii1    = ii1 + 1
             x(ii1) = atoms(a)%cx + mygran()
             y(ii1) = atoms(a)%cy + mygran()
             z(ii1) = atoms(a)%cz + mygran()
          else                     ! beta electron
             ii2    = ii2 + 1
             x(ii2) = atoms(a)%cx + mygran()
             y(ii2) = atoms(a)%cy + mygran()
             z(ii2) = atoms(a)%cz + mygran()
          endif
       enddo
    enddo
  end subroutine simpleGaussianPositions



  !--------------------------------!
  subroutine densityPositions(x,y,z)
  !--------------------------------!
    real(r8), intent(out) :: x(:)
    real(r8), intent(out) :: y(:)
    real(r8), intent(out) :: z(:)

    real(r8), parameter :: sq2 = sqrt(2.d0)
    real(r8), parameter :: sq3 = sqrt(3.d0)
    ! note matrix elements are constructed column major, i.e. transposed matrix is given
    ! alpha (a) and beta (b) positions are given separately.
    ! K shell: line with distances from origin 1
    real(r8), parameter :: vKshella(3,1) = reshape( (/ 0d0, 0d0, 1d0 /), (/3,1/) )
    real(r8), parameter :: vKshellb(3,1) = reshape( (/ 0d0, 0d0,-1d0 /), (/3,1/) )
    ! L and M shell: cube  with distances from origin sqrt(2)
    ! order: two intertwined tetrahedra, first for alpha, second for beta
    real(r8), parameter :: vLMshella(3,4) = reshape( (/ 1d0, 0d0, 1d0, &
                                                    -1d0, 0d0, 1d0, &
                                                     0d0, 1d0,-1d0, &
                                                     0d0,-1d0,-1d0 /), (/3,4/) )
    real(r8), parameter :: vLMshellb(3,4) = reshape( (/-1d0, 0d0,-1d0, &
                                                     1d0, 0d0,-1d0, &
                                                     0d0, 1d0, 1d0, &
                                                     0d0,-1d0, 1d0  /), (/3,4/) )
    ! N shell (18 e): octahedron and points over vertices  with distances from origin sqrt(2)
    real(r8), parameter :: vNshella(3,9) = reshape( (/ 0d0, 0d0, sq2, &
                                                     0d0, 0d0,-sq2, &
                                                     sq2, 0d0, 0d0, &
                                                    -1d0, 1d0, 0d0, &
                                                    -1d0,-1d0, 0d0, &
                                                    -1d0, 0d0, 1d0, &
                                                    -1d0, 0d0,-1d0, &
                                                     0d0, 1d0,-1d0, &
                                                     0d0,-1d0, 1d0  /), (/3,9/) )
    real(r8), parameter :: vNshellb(3,9) = reshape( (/ 0d0, sq2, 0d0, &
                                                     0d0,-sq2, 0d0, &
                                                    -sq2, 0d0, 0d0, &
                                                     1d0, 0d0, 1d0, &
                                                     1d0, 0d0,-1d0, &
                                                     1d0, 1d0, 0d0, &
                                                     1d0,-1d0, 0d0, &
                                                     0d0, 1d0, 1d0, &
                                                     0d0,-1d0,-1d0  /), (/3,9/) )

    ! atom shell mean radii estimated (less than 1/2 or twice ) from ELF atom shell radii (Savin/Kohout,IJQC 60,875)
    real(r8), parameter :: radiusK(36) = (/ 1.0,                                       0.9,  &
                                          0.7,  0.5,  0.35, 0.25, 0.22, 0.18, 0.16,  0.14,  &
                                          0.12, 0.11, 0.10, 0.09, 0.08, 0.075, 0.07, 0.07,  &
                                          0.07, 0.07, &
                                          0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04, 0.04, &
                                                      0.04, 0.04, 0.04, 0.03, 0.03, 0.03  /)
    real(r8), parameter :: radiusL(36) = (/ 0.0,                                       0.0, &
                                          3.0,  2.0,  1.5,  1.2,  1.0,  0.8,  0.7,   0.6, &
                                          1.2,  0.96, 0.80, 0.69, 0.60, 0.54, 0.49,  0.44, &
                                          0.40, 0.37, &
                                          0.35, 0.32, 0.30, 0.29, 0.27, 0.26, 0.24, 0.23, 0.22, 0.17, &
                                                      0.20, 0.20, 0.19, 0.18, 0.17, 0.17  /)
    real(r8), parameter :: radiusM(36) = (/ 0.0,                                       0.0, &
                                          0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0, &
                                          4.2,  3.3,  2.8,  2.3,  2.0,  1.8,  1.6,   1.4, &
                                          2.0,  1.6,  &
                                          1.5,  1.4,  1.4,  1.3,  1.3,  1.2,  1.2,  1.1, 1.1, 0.9, &
                                                      0.96, 0.87, 0.80, 0.74, 0.69, 0.64 /)
    real(r8), parameter :: radiusN(36) = (/ 0.0,                                       0.0, &
                                          0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0, &
                                          0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0, &
                                          5.0,  4.9,  &
                                          4.8,  4.6,  4.4,  4.3,  4.1,  4.0,  3.9,  3.8,  3.8, 3.6, &
                                                      3.2,  2.9,  2.6,  2.4,  2.2,  2.0  /)
    integer :: idxa,idxb,a,n,MS2
    integer :: zz, nCore
    real(r8)  :: rot(3,3)
    real(r8)  :: v(3)

    idxa = 1; idxb = nalpha+1
    MS2 = 0
    do a=1,ncenter
       nCore = wf_getNCoreElecs(a)
       if (nCore > 0) then     ! we have a PP
         zz = atoms(a)%za + nCore
         n = zz - atoms(a)%pa
         select case (nCore)
         case (2)            ! He core
            select case (n)
            case (3:9)
               call fillOpenLShell()
            case (10:18)
               call fillClosedLShell()
               call fillOpenMShell1()
            case (19:36)
               call fillClosedLShell()
               call fillClosedMShell()
               call fillOpenNShell()
            case default
               call abortp('densityPositions: He core, available for Li - Kr')
            end select
         case (10)            ! Ne core
            select case (n)
            case (11:18)
               call fillOpenMShell1()
            case (19:36)
               call fillClosedMShell()
               call fillOpenNShell()
            case default
               call abortp('densityPositions: Ne core, available for Na - Kr')
            end select
         case (18)            ! Ar core
            select case (n)
            case (19:36)
               call fillOpenNShell()
            case default
               call abortp('densityPositions: Ar core, available for K - Kr')
            end select
         case default
            call abortp('densityPositions: only He, Ne, and Ar core implemented')
         end select
       else                                   ! no PP
         zz = atoms(a)%za
         n = zz - atoms(a)%pa
         select case (n)
         case (0)
            continue
         case (1)
            call fillOpenKShell()  ! internal function
         case (2:9)
            call fillClosedKShell()
            call fillOpenLShell()
         case (10:18)
            call fillClosedKShell()
            call fillClosedLShell()
            call fillOpenMShell1()
         case (19:36)
            call fillClosedKShell()
            call fillClosedLShell()
            call fillClosedMShell()
            call fillOpenNShell()
         case default
            call abortp('densityPositions: no PP implemented only up to Kr')
         end select
       end if
    enddo
    call assert(MS2+1==mult,'densityPositions: start=density is not implemented for open K and L shells.'//&
                            ' Use start=gaussian.')
    call assert(idxa==nalpha+1 .and. idxb==ne+1,'densityPositions: internal error2')

  contains

     subroutine getRandomRotationMatrix(rotmat)
     real(r8), intent(inout) :: rotmat(:,:)
     real(r8) theta,phi,psi,costh,sinth,cosph,sinph,cosps,sinps

     ! random Euler angles
     theta = Pi*myran(); costh = cos(theta); sinth = sin(theta)
     phi = 2*Pi*myran(); cosph = cos(phi); sinph = sin(phi)
     psi = 2*Pi*myran(); cosps = cos(psi); sinps = sin(psi)

     ! (random) rotation matrix according to Goldstein, recall transposed input with reshape
     rotmat = reshape( (/ cosps*cosph-costh*sinph*sinps,  cosps*sinph+costh*cosph*sinps, sinth*sinps, &
                         -sinps*cosph-costh*sinph*cosps,  -sinps*sinph+costh*cosph*cosps, sinth*cosps, &
                          sinth*sinph,                    -sinth*cosph,                   costh /), (/3,3/) )
     end subroutine

     subroutine fillClosedShell(nn, shella, shellb, radius)
        integer, intent(in) :: nn
        real(r8), intent(in) :: shella(:,:), shellb(:,:)
        real(r8), intent(in) :: radius

        integer :: i
        real(r8) :: r
        if (nn > size(shella,2) .or. nn > size(shellb,2)) then
           call abortp("initial positions (start=density): illegal number of electrons in shell")
        endif 
        call getRandomRotationMatrix(rot)
        do i=1, nn
           r = radius*(1d0 + max(mFac*mygran(),-0.9d0))
           v = r/sq2 * matmul(rot,shella(:,i))
           x(idxa) = v(1) + atoms(a)%cx
           y(idxa) = v(2) + atoms(a)%cy
           z(idxa) = v(3) + atoms(a)%cz
           r = radius*(1d0 + max(mFac*mygran(),-0.9d0))     ! test if new r necessary
           v = r/sq2 * matmul(rot,shellb(:,i))
           x(idxb) = v(1) + atoms(a)%cx
           y(idxb) = v(2) + atoms(a)%cy
           z(idxb) = v(3) + atoms(a)%cz
           idxa = idxa + 1; idxb = idxb + 1
        enddo
     end subroutine fillClosedShell

     subroutine fillOpenShell(nn, shella, shellb, radius)
        integer, intent(in) :: nn
        real(r8), intent(in) :: shella(:,:), shellb(:,:)
        real(r8), intent(in) :: radius

        integer :: alpha, beta, i, idx
        real(r8) :: r, coord(3)

        call getRandomRotationMatrix(rot)
        alpha = 0
        beta = 0
        do i = 1, nn
           if (MS2 < mult) then
              idx = idxa + alpha
              alpha = alpha + 1
              if (alpha > size(shella,2)) then
                 call abortp("initial position (start=density): illegal open shell")
              endif
              MS2 = MS2 + 1
              coord(:) = shella(:, alpha)
           else
              idx = idxb + beta
              beta = beta + 1
              MS2 = MS2 - 1
              coord(:) = shellb(:, beta)
           endif

           r = radius * (1d0 + max(mFac*mygran(), -0.9d0))
           v = r/sq2 * matmul(rot, coord(:))

           x(idx) = v(1) + atoms(a)%cx
           y(idx) = v(2) + atoms(a)%cy
           z(idx) = v(3) + atoms(a)%cz
        enddo

        idxa = idxa + alpha
        idxb = idxb + beta
     end subroutine fillOpenShell

     subroutine fillOpenKShell()
        ! vKshella twice for consistency with previous version.
        ! since open K shell means one electron only, vKshella,vKshellb
        ! would give (due to random rotation) statistically equivalent
        ! results.
        call fillOpenShell(n, vKshella, vKshella, sq2*radiusK(zz))
     end subroutine fillOpenKShell

     subroutine fillClosedKShell()
        call fillClosedShell(1, vKshella, vKshellb, sq2*radiusK(zz))
     end subroutine fillClosedKShell

     subroutine fillOpenLShell()
        call fillOpenShell(n - 2, vLMshella, vLMshellb, radiusL(zz))
     end subroutine fillOpenLShell

     subroutine fillClosedLShell()
        call fillClosedShell(4, vLMshella, vLMshellb, radiusL(zz))
     end subroutine fillClosedLShell

     subroutine fillOpenMShell1()
        ! this uses cube arrangement for up to 8 electrons in M shell (Ar)
        call fillOpenShell(n - 10, vLMshella, vLMshellb, radiusM(zz))
     end subroutine fillOpenMShell1

     subroutine fillClosedMShell()
        call fillClosedShell(4, vLMshella, vLMshellb, radiusM(zz))
     end subroutine fillClosedMShell

     subroutine fillOpenNShell()
        call fillOpenShell(n - 18, vNshella, vNshellb, radiusN(zz))
     end subroutine fillOpenNShell
  end subroutine densityPositions



  subroutine init_LMOSampling(n,moScale)

#ifdef MPI
     use MPI_F08
#endif
     integer, intent(in) :: n       ! n+1 grid points in x,y, and z direction
     real(r8), intent(in)  :: moScale ! scaling factor for multivariate gaussians
     integer alstat
     real(r8) dist
     real(r8), allocatable :: xsum(:),ysum(:),zsum(:)
     real(r8), allocatable :: x2sum(:),y2sum(:),z2sum(:)
     real(r8), allocatable :: xysum(:),xzsum(:),yzsum(:),wsum(:)
     real(r8)              :: x(ne),y(ne),z(ne),rai(ncenter,ne)
     real(r8)              :: cov(3,3), ev(3,3), lambda(3), sigma(3,3), work(12)
     real(r8) ax,bx,ay,by,az,bz,hx,hy,hz
     real(r8) value, r(3)
     integer, allocatable:: MOIdxList(:)       ! MO index list
     integer a, i, ii, j, k, maxMOIdx, MOIdxLen, mo, lwork, ierr
     character(len=40)    :: str

     call assert(n>0 .and. moScale>0,'initLMOSampling: illegal argument')

     mMOScale = moScale

     if (MASTER .and. logmode >= 2) then
        write(iul,'(/a,i5)') '  - -  LMO sampling init  - - '
        write(iul,'(a,i4,a,g10.3)') ' grid with ',n,' points and MO scale factor = ',mMOScale
     end if

     dist = 2.d0
     call atoms_getBox(atoms,dist,ax,bx,ay,by,az,bz)
     hx = (bx-ax) / n
     hy = (by-ay) / n
     hz = (bz-az) / n

     maxMOIdx = maxval(mclist(:,1))

     if (logmode>=3) then
        write(iul,*) ' box:'
        write(iul,'(2(g12.3))') ax*bohr2angs,bx*bohr2angs,ay*bohr2angs,by*bohr2angs, &
                                az*bohr2angs,bz*bohr2angs
        write(iul,'(3(g12.3))') hx*bohr2angs,hy*bohr2angs,hz*bohr2angs
        write(iul,'(a,i6)') 'maxMOIdx = ',maxMOIdx
     end if

     allocate(mMOList(ne),mM(3,3,maxMOIdx),mMu(3,maxMOIdx), &
              MOIdxList(maxMOIdx),xsum(maxMOIdx),ysum(maxMOIdx), &
              zsum(maxMOIdx),x2sum(maxMOIdx),y2sum(maxMOIdx),z2sum(maxMOIdx), &
              xysum(maxMOIdx),xzsum(maxMOIdx),yzsum(maxMOIdx),wsum(maxMOIdx),stat=alstat)
     call assert(alstat==0,'init_LMOSampling: allocation failed')

     xsum = 0; ysum = 0; zsum = 0
     x2sum = 0; y2sum = 0; z2sum = 0
     xysum = 0; xzsum = 0; yzsum = 0
     wsum = 0
     x = 0; y = 0; z = 0
     rai = 0

     ! create list of all occurring MOs in first determinant (=mclist(:,1))
     MOIdxList = 0
     ii = 0
     do i=1,ne
        if ( .not.(any(MOIdxList(1:ii)==mclist(i,1))) ) then
           ii = ii + 1
           MOIdxList(ii) = mclist(i,1)
           mMOList(i) = ii
        end if
     end do
     MOIdxLen = ii

     if (MASTER .and. logmode >= 3) then
        write(iul,*) ' MOIdxList = '
        write(iul,'(10i4)') (MOIdxList(i),i=1,MOIdxLen)
     end if


     if (nproc < 3) then

        call calcLocal()

     else

        call calcParallel()

     end if

     deallocate(xsum,ysum,zsum,x2sum,y2sum,z2sum,xysum,xzsum,yzsum,wsum)

  CONTAINS

     subroutine calcLocal()

     integer i,j,k,a,ii

     do i=0,n
        x(1) = ax + i*hx
        do j=0,n
           y(1) = ay + j*hy
           do k=0,n
              z(1) = az + k*hz
              do a=1,ncenter
                 rai(a,1) = sqrt((x(1)-atoms(a)%cx)**2 + (y(1)-atoms(a)%cy)**2  &
                          + (z(1)-atoms(a)%cz)**2)
              end do
              if (spline) then
                 call ao1splcalc(1,x,y,z,rai)
              else
                 call ao1calc(1,x,y,z,rai)
              end if
              call mo1calc(1)
              do ii=1,MOIdxLen
                 mo = MOIdxList(ii)
                 value = mat(mo,1,1)**2
                 xsum(ii) = xsum(ii) + value*x(1)
                 x2sum(ii) = x2sum(ii) +value*x(1)*x(1)
                 ysum(ii) = ysum(ii) + value*y(1)
                 y2sum(ii) = y2sum(ii) + value*y(1)*y(1)
                 zsum(ii) = zsum(ii) + value*z(1)
                 z2sum(ii) = z2sum(ii) + value*z(1)*z(1)
                 xysum(ii) = xysum(ii) + value*x(1)*y(1)
                 xzsum(ii) = xzsum(ii) + value*x(1)*z(1)
                 yzsum(ii) = yzsum(ii) + value*y(1)*z(1)
                 wsum(ii) = wsum(ii) + value
              end do
           end do
        end do
     end do

     do ii=1,MOIdxLen
        mMu(1,ii) = xsum(ii) / wsum(ii)
        mMu(2,ii) = ysum(ii) / wsum(ii)
        mMu(3,ii) = zsum(ii) / wsum(ii)
     end do
     if (logmode >= 3) then
        write(iul,*) ' mu for each MO in list:'
        do ii=1,MOIdxLen
           mo = MOIdxList(ii)
           write(iul,'(2i5,3f10.3)') ii,mo,mMu(:,ii)*bohr2angs
        end do
        write(iul,*) ' MO types:'
        do ii=1,MOIdxLen
           mo = MOIdxList(ii)
           r(:) = mMu(:,ii)
           call atoms_whatPosition(atoms,r,str)
           write(iul,*) mo, str
        end do
        write(iul,*) ' cov matrices, eigenvalues and vectors:'
     end if

     do ii=1,MOIdxLen
        cov(1,1) = x2sum(ii)/wsum(ii) - mMu(1,ii)**2
        cov(2,2) = y2sum(ii)/wsum(ii) - mMu(2,ii)**2
        cov(3,3) = z2sum(ii)/wsum(ii) - mMu(3,ii)**2
        cov(1,2) = xysum(ii)/wsum(ii) - mMu(1,ii)*mMu(2,ii)
        cov(1,3) = xzsum(ii)/wsum(ii) - mMu(1,ii)*mMu(3,ii)
        cov(1,2) = yzsum(ii)/wsum(ii) - mMu(2,ii)*mMu(3,ii)
        cov(2,1) = cov(1,2)
        cov(3,1) = cov(1,3)
        cov(3,2) = cov(2,3)
        if (logmode >= 3) write(iul,'(3f10.3)') ((cov(i,j)*bohr2angs,j=1,3),i=1,3)
        lwork = size(work)
        ev = cov   
        call dsyev('V', 'U', 3, ev, 3, lambda, work, lwork, ierr)
        if (ierr /= 0) call abortp('init_lmoSampling: diagonalization failed')
        if (logmode >= 3) then
           write(iul,'(a,3g11.3)') ' lambda: ',lambda(:)
           write(iul,*) 'ev:'
           write(iul,'(3g11.3)') ((ev(i,j),j=1,3),i=1,3)
        end if
        call assert(minval(lambda) > 0,'init_lmoSampling: non positive eigenvalues')
        sigma = 0
        sigma(1,1) = sqrt(lambda(1)); sigma(2,2) = sqrt(lambda(2)); sigma(3,3) = sqrt(lambda(3))
        mM(:,:,ii) = matmul(ev,sigma)
     end do

     end subroutine calcLocal



     subroutine calcParallel()

     ! dynamic master/worker distribution of 3d gridpoint calculation

     integer i,j,k,a,ii,count,counter,allTasks,idx(2),srcID,vmax,workercounter
     integer, allocatable :: procIdx(:,:)
     real(r8), allocatable :: values(:),allValues(:)
#ifdef MPI
     integer                :: ierr
     type(MPI_STATUS)       :: status
#endif

     vmax = 10*MOIdxLen
     allocate(values(vmax),allValues(vmax))
     values = 0; allValues = 0

     if (MASTER) then    ! ------------ MASTER PART ------------------------------------

     allocate(procIdx(2,0:nproc-1))
     allTasks = (n+1)**2
     if (logmode>=3) write(iul,'(a,i5,a,i5,a)') ' allTasks=',allTasks,' on ',nproc,' nodes'

     counter = 0
     OUTER: do i=0,n
        MIDDLE: do j=0,n
           counter = counter + 1
           idx(1)=i; idx(2)=j
           procIdx(:,counter) = idx(:)   ! save indices
#ifdef MPI
           call MPI_SEND(idx,2,MPI_INTEGER,counter,counter,MPI_COMM_WORLD,ierr)
#endif
           if (counter == nproc-1) exit OUTER
         end do MIDDLE
      end do OUTER

      if (logmode >= 4) write(iul,*) ' master has distributed work up to: ',idx(1),idx(2)

      ! receive other results
      do ii=1,allTasks
         if (logmode >= 4) write(iul,*) ' master recv loop: ',ii
#ifdef MPI
         call MPI_RECV(values,vmax,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
         call MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,count,ierr)
         srcID = status%MPI_SOURCE
         if (logmode >= 4) write(iul,*) " master recv:",srcID,count,status%MPI_TAG
         call assert(count==vmax,' master has not received vmax data')
#else
         srcID = 0
#endif
         allValues = allValues + values
         if (counter < allTasks) then
            ! send slave new work
            counter = counter + 1
            idx(2) = idx(2) + 1
            if (idx(2)>n) then
               idx(2) = 0
               idx(1) = idx(1)+1
            end if
            procIdx(:,srcID) = idx
#ifdef MPI
            if (logmode >= 4) write(iul,'(a,2i8,i5)') "master send:",ii,counter,srcID
            call MPI_SEND(idx,2,MPI_INTEGER,srcID,counter,MPI_COMM_WORLD,ierr)
#endif
         else
            ! tell slave that he is done
#ifdef MPI
            if (logmode >= 4) write(iul,'(a,2i8,6i5)') "master send end code:",ii,srcID
            idx(1) = -1; idx(2) = -1
            call MPI_SEND(idx,2,MPI_INTEGER,srcID,counter,MPI_COMM_WORLD,ierr)
#endif
         end if
      end do

      else   !  -------  SLAVE PART  ------------------

      workercounter = 0
      do
         if (logmode >= 4) write(iull,*) "slave recv:",mytid
#ifdef MPI
         call MPI_RECV(idx,2,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
         call MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,count,ierr)
#endif
      if (idx(1)<0) exit

         ! do assigned work
         workercounter = workercounter + 1
#ifdef MPI
         if (logmode >= 4) write(iull,'(a,i5,i8,3i6)') "slave recvd:",mytid,workercounter,idx,status%MPI_TAG
#endif
         x(1) = ax + idx(1)*hx
         y(1) = ay + idx(2)*hy
         values = 0
         xsum=0; x2sum=0; ysum=0; y2sum=0; zsum=0; z2sum=0; xysum=0; xzsum=0; yzsum=0; wsum=0
         do k=0,n
            z(1) = az + k*hz
            do a=1,ncenter
              rai(a,1) = sqrt((x(1)-atoms(a)%cx)**2 + (y(1)-atoms(a)%cy)**2  &
                       + (z(1)-atoms(a)%cz)**2)
            end do
            if (spline) then
              call ao1splcalc(1,x,y,z,rai)
            else
              call ao1calc(1,x,y,z,rai)
            end if
            call mo1calc(1)
            do ii=1,MOIdxLen
              mo = MOIdxList(ii)
              value = mat(mo,1,1)**2
              xsum(ii) = xsum(ii) + value*x(1)
              x2sum(ii) = x2sum(ii) +value*x(1)*x(1)
              ysum(ii) = ysum(ii) + value*y(1)
              y2sum(ii) = y2sum(ii) + value*y(1)*y(1)
              zsum(ii) = zsum(ii) + value*z(1)
              z2sum(ii) = z2sum(ii) + value*z(1)*z(1)
              xysum(ii) = xysum(ii) + value*x(1)*y(1)
              xzsum(ii) = xzsum(ii) + value*x(1)*z(1)
              yzsum(ii) = yzsum(ii) + value*y(1)*z(1)
              wsum(ii) = wsum(ii) + value
            end do
         end do

         values(1:MOIdxLen)              = xsum
         values(MOIdxLen+1:2*MOIdxLen)   = x2sum
         values(2*MOIdxLen+1:3*MOIdxLen) = ysum
         values(3*MOIdxLen+1:4*MOIdxLen) = y2sum
         values(4*MOIdxLen+1:5*MOIdxLen) = zsum
         values(5*MOIdxLen+1:6*MOIdxLen) = z2sum
         values(6*MOIdxLen+1:7*MOIdxLen) = xysum
         values(7*MOIdxLen+1:8*MOIdxLen) = xzsum
         values(8*MOIdxLen+1:9*MOIdxLen) = yzsum
         values(9*MOIdxLen+1:10*MOIdxLen)= wsum

#ifdef MPI
         if (logmode >= 4) write(iull,'(a,2i5,g12.3)') "slave send:",mytid,idx(1),values(1)
         call MPI_SEND(values,vmax,MPI_DOUBLE_PRECISION,0,workercounter,MPI_COMM_WORLD,ierr)
#endif

      end do

      endif   ! MASTER/SLAVE

      ! now master calculates final data and broadcasts values

      if (MASTER) then

      xsum  = allValues(1:MOIdxLen)
      x2sum = allValues(MOIdxLen+1:2*MOIdxLen)
      ysum  = allValues(2*MOIdxLen+1:3*MOIdxLen)
      y2sum = allValues(3*MOIdxLen+1:4*MOIdxLen)
      zsum  = allValues(4*MOIdxLen+1:5*MOIdxLen)
      z2sum = allValues(5*MOIdxLen+1:6*MOIdxLen)
      xysum = allValues(6*MOIdxLen+1:7*MOIdxLen)
      xzsum = allValues(7*MOIdxLen+1:8*MOIdxLen)
      yzsum = allValues(8*MOIdxLen+1:9*MOIdxLen)
      wsum  = allValues(9*MOIdxLen+1:10*MOIdxLen)

      do ii=1,MOIdxLen
         mMu(1,ii) = xsum(ii) / wsum(ii)
         mMu(2,ii) = ysum(ii) / wsum(ii)
         mMu(3,ii) = zsum(ii) / wsum(ii)
      end do
      if (logmode >= 3) then
         write(iul,*) ' mu for each MO in list:'
         do ii=1,MOIdxLen
            mo = MOIdxList(ii)
            write(iul,'(2i5,3f10.3)') ii,mo,mMu(:,ii)*bohr2angs
         end do
         write(iul,*) ' MO types:'
         do ii=1,MOIdxLen
            mo = MOIdxList(ii)
            r(:) = mMu(:,ii)
            call atoms_whatPosition(atoms,r,str)
            write(iul,*) mo, str
         end do
      end if

      do ii=1,MOIdxLen
         cov(1,1) = x2sum(ii)/wsum(ii) - mMu(1,ii)**2
         cov(2,2) = y2sum(ii)/wsum(ii) - mMu(2,ii)**2
         cov(3,3) = z2sum(ii)/wsum(ii) - mMu(3,ii)**2
         cov(1,2) = xysum(ii)/wsum(ii) - mMu(1,ii)*mMu(2,ii)
         cov(1,3) = xzsum(ii)/wsum(ii) - mMu(1,ii)*mMu(3,ii)
         cov(1,2) = yzsum(ii)/wsum(ii) - mMu(2,ii)*mMu(3,ii)
         cov(2,1) = cov(1,2)
         cov(3,1) = cov(1,3)
         cov(3,2) = cov(2,3)
         if (logmode >= 3) then
            write(iul,*) 'result ii,mo=',ii,MOIdxList(ii)
            write(iul,*) ' cov matrices:'
            write(iul,'(3f10.3)') ((cov(i,j)*bohr2angs,j=1,3),i=1,3)
         end if
         lwork = size(work)
         ev = cov   
         call dsyev('V', 'U', 3, ev, 3, lambda, work, lwork, ierr)
         if (ierr /= 0) call abortp('init_lmoSampling: diagonalization failed')
         if (logmode >= 3) then
            write(iul,'(a,3g11.3)') ' lambda: ',lambda(:)
            write(iul,*) 'ev:'
            write(iul,'(3g11.3)') ((ev(i,j),j=1,3),i=1,3)
         end if
         call assert(minval(lambda)>0,'init_lmoSampling: non positive eigenvalues')
         sigma = 0
         sigma(1,1) = sqrt(lambda(1)); sigma(2,2) = sqrt(lambda(2)); sigma(3,3) = sqrt(lambda(3))
         mM(:,:,ii) = matmul(ev,sigma)
         if (logmode >= 3) then
            write(iul,*) 'mM:'
            write(iul,'(3g11.3)') ((mM(i,j,ii),j=1,3),i=1,3)
         end if
      end do

      end if ! MASTER

#ifdef MPI
      call MPI_BCAST(mMu(1,1),3*MOIdxLen,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mM(1,1,1),9*MOIdxLen,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

      if (MASTER) deallocate(procIdx)
      deallocate(values,allValues)

      if (logmode >= 4) then
         write(iul,*) ' after broadcast: mu and m for each MO in list:'
         do ii=1,MOIdxLen
            mo = MOIdxList(ii)
            write(iul,'(2i5,3f10.3)') ii,mo,mMu(:,ii)*bohr2angs
         end do
         do ii=1,MOIdxLen
            mo = MOIdxList(ii)
            write(iul,'(2i5,3f10.3)') ii,mo,mM(1,:,ii)
            write(iul,'(2i5,3f10.3)') ii,mo,mM(2,:,ii)
            write(iul,'(2i5,3f10.3)') ii,mo,mM(3,:,ii)
         end do
      end if


      end subroutine calcParallel

  end subroutine init_LMOSampling



  subroutine destroy_LMOSampling()
     integer alstat
     deallocate(mMOList,mM,mMu,stat=alstat)
     call assert(alstat==0,'destroy_LMOSampling: deallocation failed')
  end subroutine destroy_LMOSampling


  subroutine lmoPositions(x,y,z)
    real(r8), intent(out) :: x(:)
    real(r8), intent(out) :: y(:)
    real(r8), intent(out) :: z(:)
    real(r8) :: g(3,1),g1(3,1)
    integer i,ii

    call assert(allocated(mM).and.allocated(mMu),'lmoPositions: not initialized')

    do i=1,ne
       ii = mMOList(mclist(i,1))
       g(1,1) = mMOScale*mygran()
       g(2,1) = mMOScale*mygran()
       g(3,1) = mMOScale*mygran()
       g1 = matmul(mM(:,:,ii),g) + mMu(:,ii:ii)
       x(i) = g1(1,1)
       y(i) = g1(2,1)
       z(i) = g1(3,1)
    end do
  end subroutine lmoPositions


  subroutine lmoPosition(mo, r)
     integer, intent(in) :: mo
     real(r8), intent(inout) :: r(3) 
     real(r8) :: g(3,1),g1(3,1)

     call assert(allocated(mM).and.allocated(mMu),'lmoPosition: not initialized')
     
     g(1,1) = mMOScale*mygran()
     g(2,1) = mMOScale*mygran()
     g(3,1) = mMOScale*mygran()
     g1 = matmul(mM(:,:,mo),g) + mMu(:,mo:mo)
     r = g1(1:3, 1)
  end subroutine lmoPosition


  function getLMOSize() result(n)
     integer :: n
     n = size(mMu, 2)
  end function getLMOSize


end module initialPositions_m
