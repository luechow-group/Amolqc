! Copyright (C) 1996, 1999 Arne Luechow
! Copyright (C) 1999 Sebastian Manten
!
! SPDX-License-Identifier: GPL-3.0-or-later

! properties_m.f90 : module for calculations of properties/expectation values

MODULE properties_m

   use kinds_m, only: r8
   use error_m, only: error
   use rWSample_m
   use wfData_m
   use mpiInterface_m, only: myMPIAllReduceSumDouble

   implicit none

   private
   public :: propInit,setCores,future_init,futurewlk,propOutput, &
          propPrint,propCalculate,deallocFWArrays

   ! Future Walking Control Input
   integer,allocatable :: mtopsave(:) 
   integer :: dprops  = 0                     
   integer :: stprops = 0
   integer :: tau2    = 0
   integer :: notau2  = 0
   integer :: counter = 0
   integer,allocatable  :: bufcounter(:)
   integer,allocatable  :: bufoffset(:)

   integer                 :: mPropcount = 9   ! number of simultaneous properties

   real(r8)  :: ttwgt                    ! accumulative total weight
   real(r8)  :: ttwgt2(10)               ! accumulative total w.-sons

   ! distributed multipoles/polarizabilities
   !      integer,parameter    :: ndctrmax=12 ! max # of centers
   integer,private        :: mNcenter
   real(r8),allocatable     :: mDcenter(:,:)   ! x,y,z coords of global center (0), and
                                            ! distributed centers (1..ndctrmax)
   ! data structure for determination of the region (getRegion)
   real(r8),public,allocatable :: mRegpoint(:,:)
   real(r8),public,allocatable :: mRegradius(:)
   !Arrays for Props. fathers and Sons-k
   real(r8),allocatable  :: dmasum(:,:)
   real(r8),allocatable  :: dmasum2(:,:,:)
   real(r8),allocatable  :: dmaest(:,:)
   real(r8),allocatable  :: dmaest2(:,:,:)
   real(r8),allocatable  :: dmasumd(:,:,:,:)

   integer                 :: mSize = 0         ! Parameter for Array allocation
   type(RWSample), pointer :: mSample => null()


CONTAINS


   subroutine propInit(lines,nl)
      character(len=120), intent(in) :: lines(:)
      integer, intent(in)           :: nl
      integer                       :: iflag

      !gegen select case austauschen
      if(mProptype == 'dma'.or.mProptype == 'totals') then
         mPropcount = 9
         call setCenter()
         call setCores()
         if(mProptype == "totals") mNcenter = 0
      else
         call error('propInit: currently just dma+totals')
      endif
   end subroutine propInit

   subroutine propOutput(iul)
      integer,intent(in) ::  iul
      integer            ::  i,j

      selectcase(mProptype)
      case('dma') 
        write(iul,'(/A)') ' Total and Distr. Multipole Analysis:' 
      case('totals') 
        write(iul,'(/A)') ' Total Multipoles:'
      case default 
        call error("propOutput: illegal proptype")    
      endselect

      write(iul,'(A,I5,A)') ' collecting data every ',stprops,' steps'
      write(iul,'(A,I5)') ' no of datasets collected per block ',dprops
      write(iul,'(A,I5)') ' no of output-steps one block  ',notau2
      write(iul,'(A,I5)') ' interval between output-steps ',tau2

      if(mProptype == 'dma') then
         write(iul,'(I5,A)') mNcenter,' distributed centers'
         write(iul,'(A/5(3F10.6/))') ' distributed centers in bohr:',       & 
                        ((mDcenter(i,j),i=1,3),j=0,mNcenter)
         write(iul,'(A/5(4F10.6/))') ' Region definition points and radii', &
                     ((mRegpoint(i,j),i=1,3),mRegradius(j),j=1,mNcenter)
      endif
   end subroutine propOutput


   subroutine future_init(sample,mStprops,mTau2,mNotau2)
      type(RWSample), intent(in),target :: sample
      integer, intent(in)               :: mStprops, mTau2, mNotau2  
      integer :: j, alstat, mDprops    

      mSample => sample

      if(.not.associated(mSample,sample)) call error('futurewlk: mSample not associated')

      mSize    = getMaxSampleSize(mSample) ! To ensure that the array is big enough -> Branching
      mDprops  = nDWGen()

      dprops   = mDprops    
      stprops  = mStprops 
      tau2     = mTau2  
      notau2   = mNotau2

      call allocArrays()
      call resetArrays()
          
      do j = 1,dprops
         bufoffset(j) = (j-1)*stprops         
      enddo   
      counter  = 0

      mSample => null()

   end subroutine future_init       

   subroutine futurewlk(sample)

      type(RWSample), intent(inout), target :: sample
      integer :: j,t,lastRW


      if(associated(mSample,sample)) call error('futurewlk: mSample => inconsistent')  
      mSample => sample 
      if(.not.associated(mSample,sample)) call error('futurewlk: mSample not associated')

      counter = counter + 1
      do j = 1,dprops 
         if (counter.gt.bufoffset(j)) then
            bufcounter(j) = bufcounter(j) + 1
            if (bufcounter(j) == 1) then
               mtopsave(j) = getSampleSize(mSample)
               call SaveFathers(mSample,j)
               call propCollect(j)
            endif   
            if (bufcounter(j) == (stprops*dprops)) then ! one step before new counter-start  
               bufcounter(j) = 0
               lastRW = mtopsave(j) 
               call propCollectw(lastRW,j,10)
            elseif (mod(bufcounter(j),tau2) == 0) then    
               if ((bufcounter(j)/tau2) <= notau2) then
                  t = bufcounter(j)/tau2
                  lastRW = mtopsave(j)
                  call propCollectw(lastRW,j,t)
               endif
            endif   
         endif
      enddo
           
      mSample => null()      

   end subroutine futurewlk

   subroutine propCollect(ix)
   !     input parameter:
   integer  ix
   type(randomWalker), pointer       :: rwp  => null()
   integer  ps,a,i,j,wi,idx,ierr,m
   real(r8)   xx,x0,yy,y0,zz,z0,rr2,r02
   real(r8)   x(nmax),y(nmax),z(nmax)  ! x,y,z coord of walker
   real(r8)   w,twgt ! walker weight(s)
   real(r8)   dsuma(0:mPropcount,0:mNcenter),dsum(0:mPropcount,0:mNcenter)

   ! distributed dipole and quadrupole moments
   ! order: 0=q, 1,2,3 = x,y,z, 4-9=xx,yy,zz,xy,xz,yz

   dsuma   = 0d0

   if(mProptype == 'dma'.or.mProptype == 'totals') then
   do a=1,ncenter

   x0 = atoms(a)%cx - mDcenter(1,0)
   y0 = atoms(a)%cy - mDcenter(2,0)
   z0 = atoms(a)%cz - mDcenter(3,0)
   r02= x0*x0 + y0*y0 + z0*z0
   rr2= xx*xx + yy*yy + zz*zz

   ! total moments
   dsuma(0,0)   = dsuma(0,0)   + atoms(a)%za    ! charge
   dsuma(1,0)   = dsuma(1,0)   + atoms(a)%za*x0 ! dipole
   dsuma(2,0)   = dsuma(2,0)   + atoms(a)%za*y0
   dsuma(3,0)   = dsuma(3,0)   + atoms(a)%za*z0
   dsuma(4,0)   = dsuma(4,0)   + 0.5d0*atoms(a)%za*(3d0*x0*x0 - r02) ! quadpole
   dsuma(5,0)   = dsuma(5,0)   + 0.5d0*atoms(a)%za*(3d0*y0*y0 - r02)
   dsuma(6,0)   = dsuma(6,0)   + 0.5d0*atoms(a)%za*(3d0*z0*z0 - r02)
   dsuma(7,0)   = dsuma(7,0)   + 1.5d0*atoms(a)%za*x0*y0
   dsuma(8,0)   = dsuma(8,0)   + 1.5d0*atoms(a)%za*x0*z0
   dsuma(9,0)   = dsuma(9,0)   + 1.5d0*atoms(a)%za*y0*z0

   if(mProptype == 'dma') then

    idx = getRegion(atoms(a)%cx,atoms(a)%cy,atoms(a)%cz)
    xx = atoms(a)%cx - mDcenter(1,idx)
    yy = atoms(a)%cy - mDcenter(2,idx)
    zz = atoms(a)%cz - mDcenter(3,idx)

   ! distributed moments (contribution to region idx)
   dsuma(0,idx) = dsuma(0,idx) + atoms(a)%za
   dsuma(1,idx) = dsuma(1,idx) + atoms(a)%za*xx
   dsuma(2,idx) = dsuma(2,idx) + atoms(a)%za*yy
   dsuma(3,idx) = dsuma(3,idx) + atoms(a)%za*zz
   dsuma(4,idx) = dsuma(4,idx) + 0.5d0*atoms(a)%za*(3d0*xx*xx - rr2)
   dsuma(5,idx) = dsuma(5,idx) + 0.5d0*atoms(a)%za*(3d0*yy*yy - rr2)
   dsuma(6,idx) = dsuma(6,idx) + 0.5d0*atoms(a)%za*(3d0*zz*zz - rr2)
   dsuma(7,idx) = dsuma(7,idx) + 1.5d0*atoms(a)%za*xx*yy
   dsuma(8,idx) = dsuma(8,idx) + 1.5d0*atoms(a)%za*xx*zz
   dsuma(9,idx) = dsuma(9,idx) + 1.5d0*atoms(a)%za*yy*zz
   endif
   enddo

   rwp => getFirst(mSample)
       
   m = getSampleSize(mSample)
   	     
   do ps = 1,m 

   call pos(rwp,x,y,z)
   w = wgt(rwp)

   dsum = 0d0
   do i=1,ne
      
      x0 = x(i) - mDcenter(1,0)
      y0 = y(i) - mDcenter(2,0)
      z0 = z(i) - mDcenter(3,0)            
      r02= x0*x0 + y0*y0 + z0*z0
      rr2= xx*xx + yy*yy + zz*zz

      ! total moments
      dsum(0,0)   = dsum(0,0) - 1.d0 ! charge
      dsum(1,0)   = dsum(1,0) - x0   ! dipole
      dsum(2,0)   = dsum(2,0) - y0
      dsum(3,0)   = dsum(3,0) - z0
      dsum(4,0)   = dsum(4,0) - 0.5d0*(3d0*x0*x0 - r02) ! quadpole
      dsum(5,0)   = dsum(5,0) - 0.5d0*(3d0*y0*y0 - r02)
      dsum(6,0)   = dsum(6,0) - 0.5d0*(3d0*z0*z0 - r02)
      dsum(7,0)   = dsum(7,0) - 1.5d0*x0*y0
      dsum(8,0)   = dsum(8,0) - 1.5d0*x0*z0
      dsum(9,0)   = dsum(9,0) - 1.5d0*y0*z0
      ! distributed moments (contribution to region idx)
    if(mProptype == 'dma') then

      idx = getRegion(x(i),y(i),z(i))
      xx = x(i) - mDcenter(1,idx)
      yy = y(i) - mDcenter(2,idx)
      zz = z(i) - mDcenter(3,idx)

      dsum(0,idx) = dsum(0,idx) - 1.d0
      dsum(1,idx) = dsum(1,idx) - xx
      dsum(2,idx) = dsum(2,idx) - yy
      dsum(3,idx) = dsum(3,idx) - zz
      dsum(4,idx) = dsum(4,idx) - 0.5d0*(3d0*xx*xx - rr2)
      dsum(5,idx) = dsum(5,idx) - 0.5d0*(3d0*yy*yy - rr2)
      dsum(6,idx) = dsum(6,idx) - 0.5d0*(3d0*zz*zz - rr2)
      dsum(7,idx) = dsum(7,idx) - 1.5d0*xx*yy
      dsum(8,idx) = dsum(8,idx) - 1.5d0*xx*zz
      dsum(9,idx) = dsum(9,idx) - 1.5d0*yy*zz
    endif
   enddo

   dmasumd(:,:,ps,ix) = dsuma(:,:)  + dsum(:,:)

      dmasum(:,:) = dmasum(:,:) + w*dmasumd(:,:,ps,ix)
      ttwgt       = ttwgt       + w

   if (.not.isNext(mSample)) exit
   rwp => getNext(mSample)

   enddo
   endif ! total or dma
   rwp => null()

   end subroutine propCollect


   subroutine propCollectw(lastW,ix,ntau2)
   !     input parameter:
   integer  ix,ntau2,alstat
   integer  wi,j
   integer  lastW                      ! last walker
   real(r8),allocatable  ::  ws(:)       ! weight-sum of sons

   allocate(ws(mSize),stat=alstat)
   if (alstat /= 0) call error("Properties::proCollectw: allocate ws failed")
   ws(:) = 0

   call GetSonWeight(ix,mSample,ws) 

   do wi = 1,lastW
      dmasum2(:,:,ntau2) = dmasum2(:,:,ntau2)           &
                          + ws(wi)*dmasumd(:,:,wi,ix) 
      ttwgt2(ntau2)      = ttwgt2(ntau2) + ws(wi)
   enddo

   deallocate(ws)
   end subroutine propCollectw


   subroutine propCalculate
   integer :: i,ierr,j
   integer :: nr1,nr2,nr3,nr4
   real(r8)  ::sumttwgt,sumttwgt2(10)

   ! the following two lines are using sequence association
   call myMPIAllReduceSumDouble(dmasum(1,1),dmaest(1,1),size(dmasum))
   call myMPIAllReduceSumDouble(dmasum2(1,1,1),dmaest2(1,1,1),size(dmasum2))
   call myMPIAllReduceSumDouble(ttwgt,sumttwgt,1)
   call myMPIAllReduceSumDouble(ttwgt2,sumttwgt2,size(ttwgt2))

   if (nproc == 1) then
   dmaest    = dmasum
   sumttwgt  = ttwgt
   dmaest2   = dmasum2
   sumttwgt2 = ttwgt2
   endif 

   if (mytid == 0) then
      dmaest(:,:)  = dmaest(:,:)/sumttwgt
      do i = 1,notau2
         dmaest2(:,:,i) = dmaest2(:,:,i)/sumttwgt2(i)
      enddo   
      dmaest2(:,:,10) = dmaest2(:,:,10)/sumttwgt2(10)         
   endif

   end subroutine propCalculate

   subroutine propPrint
   integer :: i,j,k,nd
   !    real(r8)  :: tmp1,tmp2,tmp3

   if (mytid == 0) then
   write(iul,'(/A)') ' multipole moments as expectation values'
   write(iul,*)      ''
   write(iul,'(/A)') ' without descended weighting'

   write(iul,'(/A,F20.10)') 'total charge = ',dmaest(0,0)
   write(iul,'(/A,3(/A,F8.5,A,F8.5,A))') ' total dipole moments:',  &
   'mu_x = ',dmaest(1,0),' a.u. = ',dmaest(1,0)*debye,' D',         &
   'mu_y = ',dmaest(2,0),' a.u. = ',dmaest(2,0)*debye,' D',         &
   'mu_z = ',dmaest(3,0),' a.u. = ',dmaest(3,0)*debye,' D'
   write(iul,'(/A,6(/A,F8.5))') ' total quadrupole moments [a.u.] =', &                           
   'xx = ',dmaest(4,0),'yy = ',dmaest(5,0),                         &   
   'zz = ',dmaest(6,0),'xy = ',dmaest(7,0),                         &
   'xz = ',dmaest(8,0),'yz = ',dmaest(9,0) 
   if(mProptype == 'dma') then
   do nd=1,mNcenter
   write(iul,'(/A,F20.10)') 'partial charge = ',dmaest(0,nd)
   write(iul,'(/A,I5,3(/2(A,F8.5),A))')                             & 
    ' distrib. dipole moments =',nd,                                & 
    'mu_x = ',dmaest(1,nd),' a.u. = ',dmaest(1,nd)*debye,' D',      &
    'mu_y = ',dmaest(2,nd),' a.u. = ',dmaest(2,nd)*debye,' D',      &  
    'mu_z = ',dmaest(3,nd),' a.u. = ',dmaest(3,nd)*debye,' D'      
   write(iul,'(/A,I5,6(/A,F8.5))')                                  &
    ' distrib. quadrupole moments [a.u.] =',nd,                     &
    'xx = ',dmaest(4,nd),'yy = ',dmaest(5,nd),                      &
    'zz = ',dmaest(6,nd),'xy = ',dmaest(7,nd),                      &
    'xz = ',dmaest(8,nd),'yz = ',dmaest(9,nd)                      
   enddo
   endif 

   write(iul,*)      ''
   write(iul,*)      ''
   write(iul,'(/A)') ' with descended weighting'

   do k = 1,10
   if ((k.eq.10).or.(k.le.notau2)) then
   write(iul,*) ''
   write(iul,*) ''
   if(k.eq.10) then 
   write(iul,'(A20,I8)') 'number of steps:',stprops*dprops  
   else   
   write(iul,'(A20,I8)') 'number of steps:',k*tau2  
   endif   


   write(iul,'(/A,F20.10)') 'total charge = ',dmaest2(0,0,k)
   write(iul,'(/A,3(/A,F8.5,A,F8.5,A))') ' total dipole moments:',    &
   'mu_x =',dmaest2(1,0,k),' a.u. =',dmaest2(1,0,k)*debye,' D',       &
   'mu_y =',dmaest2(2,0,k),' a.u. =',dmaest2(2,0,k)*debye,' D',       &
   'mu_z =',dmaest2(3,0,k),' a.u. =',dmaest2(3,0,k)*debye,' D'
   write(iul,'(/A,6(/A,F8.5))') 'total quadrupole moments [a.u.] =',  &
   'xx = ',dmaest2(4,0,k),'yy = ',dmaest2(5,0,k),                     &
   'zz = ',dmaest2(6,0,k),'xy = ',dmaest2(7,0,k),                     &
   'xz = ',dmaest2(8,0,k),'yz = ',dmaest2(9,0,k)                      
   if(mProptype == 'dma') then
   do nd=1,mNcenter
   write(iul,'(/A,F20.10)') 'partial charge = ',dmaest2(0,nd,k)
   write(iul,'(/A,I5,3(/2(A,F8.5),A))') 'distrib. dipole moments =',nd,  & 
   'mu_x =',dmaest2(1,nd,k),' a.u. =',dmaest2(1,nd,k)*debye,' D',     &
   'mu_y =',dmaest2(2,nd,k),' a.u. =',dmaest2(2,nd,k)*debye,' D',     &
   'mu_z =',dmaest2(3,nd,k),' a.u. =',dmaest2(3,nd,k)*debye,' D'
   write(iul,'(/A,I5,6(/A,F8.5))')                                    & 
   'distrib. quadrupole moments [a.u.] =',nd,                       &
   'xx = ',dmaest2(4,nd,k),'yy = ',dmaest2(5,nd,k),                  &
   'zz = ',dmaest2(6,nd,k),'xy = ',dmaest2(7,nd,k),                  &
   'xz = ',dmaest2(8,nd,k),'yz = ',dmaest2(9,nd,k)

   enddo
   endif ! dma if

   endif  ! if k 1to10
   enddo ! k do

   endif !mytid
   end subroutine propPrint


   integer function getRegion(x,y,z)
   ! This is no good solution for the region problem -> s. Baders AIM
   real(r8), intent(in):: x,y,z

   integer :: i,iMin
   real(r8)  :: distance(ncenter),minDist

   do i=1,ncenter
   distance(i) = sqrt( (x-mRegpoint(1,i))**2                        &  
    + (y-mRegpoint(2,i))**2 + (z-mRegpoint(3,i))**2 )                &
    / mRegradius(i)
   enddo
   minDist = distance(1)
   iMin    = 1
   do i=2,ncenter
   if (distance(i) < minDist) then
      minDist = distance(i)
      iMin = i
   endif
   enddo
   getRegion = iMin
   end function getRegion



!-------------Set/Get/Alloc Routines ---------------------------------------------------------------!

   subroutine setCores()
      integer          :: j,alstat
       allocate(mDcenter(3,0:ncenter),stat=alstat)
       if (alstat /= 0) call error("Properties:propInput: allocate dCenter failed")

        mDcenter(:,:) = 0.0d0

       do j=0,ncenter
         if(j == 0) then
          mDcenter(1,j) = 0.0d0
          mDcenter(2,j) = 0.0d0
          mDcenter(3,j) = 0.0d0
         else
          mDcenter(1,j) = atoms(j)%cx
          mDcenter(2,j) = atoms(j)%cy
          mDcenter(3,j) = atoms(j)%cz
         endif
      enddo 
   end subroutine setCores

   subroutine setCenter()
      mNcenter = ncenter
   end subroutine setCenter


   subroutine allocArrays()
      integer :: iflag,alstat
      allocate(mtopsave(dprops),stat=alstat)
      if (alstat /= 0) call error("Properties:future_init: allocate topsave failed") 
      allocate(bufcounter(dprops),stat=alstat)
      if (alstat /= 0) call error("Properties:future_init: allocate bufcount failed")
      allocate(bufoffset(dprops),stat=alstat)
      if (alstat /= 0) call error("Properties:future_init: allocate bufcounter failed")
      allocate(dmasum(0:mPropcount,0:mNcenter),stat=alstat)
      if (alstat /= 0) call error("Properties:future_init: allocate dmasum failed")
      allocate(dmasum2(0:mPropcount,0:mNcenter,10),stat=alstat)
      if (alstat /= 0) call error("Properties:future_init: allocate dmasum2 failed")
      allocate(dmaest(0:mPropcount,0:mNcenter),stat=alstat)
      if (alstat /= 0) call error("Properties:future_init: allocate dmaest failed")
      allocate(dmaest2(0:mPropcount,0:mNcenter,10),stat=alstat) 
      if (alstat /= 0) call error("Properties:future_init: allocate dmaest2 failed") 
      allocate(dmasumd(0:mPropcount,0:mNcenter,mSize,dprops),stat=alstat)
      if (alstat /= 0) call error("Properties:future_init: allocate dmasumd failed")
   end subroutine allocArrays 


   subroutine deallocFWArrays()
      integer :: iflag,alstat
      deallocate(mtopsave)      
      deallocate(bufcounter)
      deallocate(bufoffset)
      deallocate(dmasum)
      deallocate(dmasum2)
      deallocate(dmaest)
      deallocate(dmaest2) 
      deallocate(dmasumd)
   end subroutine deallocFWArrays 


   subroutine resetArrays() 
      mtopsave(:)    = 0 
      bufcounter(:)  = 0
      bufoffset(:)   = 0              
      dmasum(:,:)      = 0.0d0
      dmasum2(:,:,:)   = 0.0d0
      dmaest(:,:)      = 0.0d0
      dmaest2(:,:,:)   = 0.0d0
      dmasumd(:,:,:,:) = 0.0d0
   end subroutine resetArrays

   
END MODULE properties_m



subroutine propInput(proptype,ncenters,iu)
   use error_m
   use properties_m
   implicit none
   integer, intent(in)          :: iu       ! open file unit for reading input
   integer, intent(in)          :: ncenters
   character(len=6), intent(in) :: proptype
   integer                      :: i,j,alstat

   if(proptype == 'dma') then
      allocate(mRegpoint(3,ncenters),mRegradius(ncenters),stat=alstat)
      if (alstat /= 0) call error("Properties:propInput: allocate reg-Arrays failed failed")
      mRegradius(:) = 0.0d0
      mRegpoint(:,:) = 0.0d0 
      do j=1,ncenters  
       read(iu,*) mRegpoint(1,j), mRegpoint(2,j), mRegpoint(3,j),mRegradius(j)
      enddo
   endif
end subroutine propInput


