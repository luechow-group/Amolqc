! Copyright (C) 2006-2007, 2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

MODULE rWStatistics_m

! statistics (mean,stddev,cov) for random walkers positions (i.e. 3*ne dim)
! in cartesian or spherical coordinates

  use kinds_m, only: r8, i8
  use randomWalker_m
  use vectorStatistics_m, vs_add => add

  implicit none

  private
  public  :: rwStatInit, rwStatReset, rwAddToStat, rwStatMean, rwStatStdDev, rwStatPrint

  type(vectorStat),save  :: mVStat                ! vector statistics object
  character*1            :: mType='c'             ! type for rw statistics


CONTAINS

  !-------------------------!
  subroutine rwStatInit(type)
  !-------------------------!

    character*1, intent(in) :: type

    if (type == 's') then
       mType = 's'           ! statistics for spherical coordinates
    else
       mType = 'c'           ! statistics for cartesian coordinates
    endif
    
    call new(mVStat,3*ne)
    call reset(mVStat)
    
  end subroutine rwStatInit


  !----------------------!
  subroutine rwStatReset()
  !----------------------!
    
    call reset(mVStat)
    
  end subroutine rwStatReset


  !-------------------------!
  subroutine rwAddToStat(rwp)
  !-------------------------!

    type(RandomWalker), pointer :: rwp
    real(r8)    :: v(3*ne),vmin,vmax
    real(r8)    :: vx(ne), vy(ne), vz(ne)
    integer   :: i

    call pos(rwp,vx,vy,vz)
    if (mType == 'c') then
       do i=1,ne
          v(3*i-2) = vx(i)
          v(3*i-1) = vy(i)
          v(3*i)   = vz(i)
       enddo
    else
       do i=1,ne
          call cartesian2spherical(vx(i),vy(i),vz(i),v(3*i-2),v(3*i-1),v(3*i))
       enddo
    endif

!      vmin = minval(v)
!      vmax = maxval(v)
!      if (vmin < -10.0 .or. vmax > 10.) then
!         write(*,*) 'rwAddToStat:', currw
!         write(*,'(6(G10.3))') (rwArray(currw)%x(i),rwArray(currw)%y(i),
!     .               rwArray(currw)%z(i),i=1,ne)
!         write(*,'(6(G10.3))') (v(i),i=1,3*ne)
!      endif

    call vs_add(mVStat,v)

  end subroutine rwAddToStat


  !-------------------!
  function rwStatMean()
  !-------------------!

    real(r8)     :: rwStatMean(3*ne)
      
    rwStatMean = mean(mVStat)
    
  end function rwStatMean


  !---------------------!
  function rwStatStdDev()
  !---------------------!

    real(r8)     :: rwStatStdDev(3*ne)
      
    rwStatStdDev = sqrt(variance(mVStat))

  end function rwStatStdDev

  !------------------------------!
  subroutine rwStatPrint(iu,step)
  !------------------------------!
      
    integer, intent(in) :: iu
    integer(i8), intent(in) :: step
    integer             :: i,j
    real(r8)              :: vmean(3*ne), vstddev(3*ne), vcov(3*ne,3*ne)

    vmean = mean(mVStat)
    vstddev = sqrt(variance(mVStat))
    write(iu,'(A14,I15)') 'step: ',step
    if (mType == 'c') then
       write(iu,*) ' Random Walker position statistics in cartesian coords'
    else
       write(iu,*) ' Random Walker position statistics in spherical coords'
    endif
    write(iu,*) ' data count = ',dataCount(mVStat)
    write(iu,*) ' mean +/- std dev (individual) in a.u.:'
    write(iu,'(3(G12.3,A,G12.3))') (vmean(i),'+-',vstddev(i),i=1,3*ne)
    vcov = covariance(mVStat)
    write(iu,*) ' cov(i,j) in a.u.:'
    do i=1,3*ne
       write(iu,fmt='(I3,1X)',advance='no') i
       write(iu,'(10(G12.3))') (vcov(i,j),j=1,i)
    enddo
    
  end subroutine rwStatPrint


      
END MODULE rWStatistics_m


!====================
! utility subroutines
!====================




subroutine cartesian2spherical(x,y,z,r,theta,phi)
  use kinds_m, only: r8
  implicit none
  real(r8) x,y,z,r,theta,phi
  r = sqrt(x**2+y**2+z**2)
  theta = acos(z/r)
  phi = atan2(y,x)
end subroutine cartesian2spherical


subroutine spherical2cartesian(r,theta,phi,x,y,z)
  use kinds_m, only: r8
  implicit none
  real(r8) x,y,z,r,theta,phi
  x = r*sin(theta)*cos(phi)
  y = r*sin(theta)*sin(phi)
  z = r*cos(theta)
end subroutine spherical2cartesian
      


