! Copyright (C) 2006-2007 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later


! simple statistics (mean, covariance) for vectors

!============================
module vectorStatistics_m
!============================

! $Id: vectorStatistics_m.f90,v 1.1.1.1 2007/04/25 13:42:19 luechow Exp $

! $Log: vectorStatistics_m.f90,v $
! Revision 1.1.1.1  2007/04/25 13:42:19  luechow
! QMC program amolqc. rewritten in 2006. AL
!

!
!  TODO: how to check size of a matrix?  assert in covariance
!

  use kinds_m, only: r8
  use error_m
  implicit none

  type vectorStat
     private
     real(r8), pointer   :: dAverageSum(:)
     real(r8), pointer   :: dCovarianceSum(:,:)
     integer           :: dVectorSize
     integer           :: dDataCounter = 0
  end type vectorStat

  interface assignment(=)
     module procedure vectorStat_assign
  end interface

  interface operator(+)
     module procedure vectorStat_add
  end interface

contains

  !--------------------!
  subroutine new(self,n)
  !--------------------!

    type(vectorStat), intent(inout) :: self
    integer, intent(in)             :: n       ! dimension of vectors

    call assert(n>0,"vectorStat_new: n>0 required")
    self%dVectorSize = n
    allocate(self%dAverageSum(n),self%dCovarianceSum(n,n))
  end subroutine new

  !---------------------!
  subroutine delete(self)
  !---------------------!

    type(vectorStat), intent(inout) :: self
    if (associated(self%dAverageSum)) deallocate(self%dAverageSum)
    if (associated(self%dCovarianceSum)) deallocate(self%dCovarianceSum)
  end subroutine delete

  !-----------------------!
  subroutine add(self,data)
  !-----------------------!

    type(vectorStat), intent(inout) :: self
    real(r8), intent(in)              :: data(:)
    integer                         :: i,j
    call assert(size(data)==self%dVectorSize,"vectorStat_add: wrong size")
    self%dAverageSum = self%dAverageSum + data
    do i=1,self%dVectorSize
       do j=1,i
          self%dCovarianceSum(j,i) = self%dCovarianceSum(j,i) + data(j)*data(i)
       enddo
    enddo
    self%dDataCounter = self%dDataCounter + 1

  end subroutine add


  !--------------------!
  subroutine reset(self)
  !--------------------!

    type(vectorStat), intent(inout) :: self
    self%dAverageSum = 0
    self%dCovarianceSum = 0
    self%dDataCounter = 0
  end subroutine reset


  !-------------------------------------!
  subroutine vectorStat_assign(self,stat)
  !-------------------------------------!

    type(vectorStat), intent(out) :: self
    type(vectorStat), intent(in)  :: stat
    self%dAverageSum = stat%dAverageSum
    self%dCovarianceSum = stat%dCovarianceSum
    self%dDataCounter = stat%dDataCounter
  end subroutine vectorStat_assign


  !------------------------------------------------!
  type(vectorStat) function vectorStat_add(self,add)
  !------------------------------------------------!

    type(vectorStat), intent(in)  :: self
    type(vectorStat), intent(in)  :: add

    call new(vectorStat_add, size(self%dAverageSum))
    vectorStat_add%dAverageSum = self%dAverageSum + add%dAverageSum
    vectorStat_add%dCovarianceSum = self%dCovarianceSum + add%dCovarianceSum
    vectorStat_add%dDataCounter = self%dDataCounter + add%dDataCounter
  end function vectorStat_add


  ! access functions


  !-----------------!
  function mean(self)
  !-----------------!

    type(vectorStat), intent(in)  :: self
    real(r8)                        :: mean(self%dVectorSize)
    call assert(self%dDataCounter>0,"vectorStat_mean: no data available")
    mean = self%dAverageSum / self%dDataCounter
  end function mean


  !---------------------!
  function variance(self)
  !---------------------!

    type(vectorStat), intent(in)  :: self
    real(r8)                        :: variance(self%dVectorSize)
    integer                       :: i
    call assert(self%dDataCounter>0,"vectorStat_variance: no data available")
    do i=1,self%dVectorSize
       variance(i) = self%dCovarianceSum(i,i)/self%dDataCounter &
                     - (self%dAverageSum(i)/self%dDataCounter)**2
    enddo
  end function variance

  !-----------------------!
  function covariance(self)
  !-----------------------!

    type(vectorStat), intent(in)  :: self
    real(r8)                        :: covariance(self%dVectorSize,self%dVectorSize)
    real(r8)                        :: mn(self%dVectorSize)
    integer                       :: i,j
    call assert(self%dDataCounter>0,"vectorStat_covariance: no data available")
    mn = mean(self)
    do i=1,self%dVectorSize
       do j=1,i
          covariance(j,i) = self%dCovarianceSum(j,i)/self%dDataCounter - mn(j)*mn(i)
          covariance(i,j) = covariance(j,i)
       enddo
    enddo
  end function covariance

  !-----------------------!
  function stdDevMean(self)
  !-----------------------!

    type(vectorStat), intent(in)  :: self
    real(r8)                        :: stdDevMean(self%dVectorSize)
    integer                       :: i
    call assert(self%dDataCounter>1,"vectorStat_stdDevMean: not enough data")
    stdDevMean = variance(self)
    do i=1,self%dVectorSize
       stdDevMean(i) = sqrt( stdDevMean(i) / (self%dDataCounter - 1) )
    enddo
  end function stdDevMean


  !------------------------------!
  integer function dataCount(self)
  !------------------------------!

    type(vectorStat), intent(in)  :: self
    dataCount = self%dDataCounter
  end function dataCount


  !----------------------------
  subroutine serialize(self,iu)
  !----------------------------

    ! write object unformatted to open file unit iu
    type(vectorStat), intent(in)  :: self
    integer, intent(in)           :: iu    ! open file unit
    write(iu) self%dVectorSize,self%dDataCounter
    write(iu) self%dAverageSum,self%dCovarianceSum
  end subroutine serialize


  !------------------------------
  subroutine deSerialize(self,iu)
  !------------------------------

    ! read object unformatted from open file unit iu
    type(vectorStat), intent(inout)  :: self
    integer, intent(in)              :: iu    ! open file unit
    read(iu) self%dVectorSize,self%dDataCounter
    if (associated(self%dAverageSum)) deallocate(self%dAverageSum)
    if (associated(self%dCovarianceSum)) deallocate(self%dCovarianceSum)
    allocate(self%dAverageSum(self%dVectorSize),self%dCovarianceSum(self%dVectorSize,self%dVectorSize))
  end subroutine deSerialize


  !--------------------------
  subroutine debugPrint(self)
  !--------------------------
    
    type(vectorStat), intent(in)  :: self
    print*,"debugPrint:dataCounter=",self%dDataCounter
    print*,"debugPrint:vectorSize =",self%dVectorSize
    write(*,'(A,/,100(10G10.3))') "debugPrint:averageSum:",self%dAverageSum
    write(*,'(A,/,100(10G10.3))') "debugPrint:covarianceSum=",self%dCovarianceSum
  end subroutine debugPrint

end module vectorStatistics_m

