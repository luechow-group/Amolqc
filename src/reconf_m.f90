! Copyright (C) 2009, 2013, 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

Module reconfg_m

use kinds_m, only: r8
use error_m
use rWSample_m
implicit none


  private
  public      :: qmc_reconfInit,qmc_reconfDestroy,qmc_reconfNew,qmc_reconf1,qmc_reconf2
   
  integer     :: mN = 0
  integer     :: mMode = 0     ! 1: reset weights to mean  2: reset weights to one
  integer, allocatable :: mIdxPlus(:),mIdxMinus(:)
  real(r8), allocatable  :: mWPlus(:),mWMinus(:)
  real(r8), allocatable  :: mWSumPlus(:),mWSumMinus(:)

  type Weighta
    integer :: ip = 0 
    real(r8)  :: wp = 0.d0
    integer :: in = 0
    real(r8)  :: wn = 0.d0
  end type Weighta



contains

  subroutine qmc_reconfInit(sample,mode)
    type(RWSample), intent(inout) :: sample
    integer, intent(in) :: mode

    mN = getSampleSize(sample)
    mMode = 1
    if (mode==4 .or. mode==6) mMode = 2
    allocate(mIdxPlus(mN),mIdxMinus(mN),mWPlus(mN),mWMinus(mN),mWSumPlus(mN),mWSumMinus(mN))
    mIdxPlus = 0
    mIdxMinus = 0
    mWPlus = 0
    mWMinus = 0
    mWSumPlus = 0
    mWSumMinus = 0
  end subroutine qmc_reconfInit

  subroutine qmc_reconfDestroy()
    deallocate(mIdxPlus,mIdxMinus,mWPlus,mWMinus,mWSumPlus,mWSumMinus)
  end subroutine qmc_reconfDestroy

  subroutine qmc_reconfNew(sample)
  !------------------------------!
    ! implements Caffarel's reconfiguration as described by Stuart Rothstein JCP 2004
    ! notation as in paper
    ! this implementation deletes killed walkers (with w=0) and surely duplicates walkers
    ! with relative weight >= 2
    ! the search for walkers to duplicate/replace can be sped up with binary search
    type(RWSample), intent(inout) :: sample
    type(RandomWalker), pointer   :: rwp => null()
    integer i,k,nPlus,nMinus,Nreconf,kPlus,kMinus
    real(r8) Wmean,wPSum,wMSum,wi,xi

    call assert(getSampleSize(sample)==mN,"qmc_reconfNew: sample size has changed")

    Wmean = getSampleWeight(sample) / mN
    !!print*,"start with Wmean=",Wmean

    wPSum = 0; wMSum = 0
    nPlus = 0; nMinus = 0
    rwp => getFirst(sample)
    do i=1,mN
      wi = wgt(rwp)/Wmean
      if (wi > 1.d0) then
        nPlus = nPlus + 1
        mIdxPlus(nPlus) = i
        mWPlus(nPlus) = wi - 1.d0
        wPSum = wPSum + wi - 1.d0
        mWSumPlus(nPlus) = wPSum
      else if (wi < 1.d0) then
        nMinus = nMinus + 1
        mIdxMinus(nMinus) = i
        mWMinus(nMinus) = 1.d0 - wi
        wMSum = wMSum + 1.d0 - wi
        mWSumMinus(nMinus) = wMSum
      end if
      rwp => getNext(sample)
    end do

    !!print*,wPSum,wMSum,nPlus,nMinus
    call assertEqualAbsolute(wPSum,wMSum,1.d-6,"qmc_reconfNew: plus sum unequal minus sum")

    Nreconf = floor(wPSum+myran())

    !!print*,Nreconf
    do i=1,Nreconf
      ! pick plus walker to duplicate
      xi = wPsum*myran()
      do k=1,nPlus
        if (xi <= mWSumPlus(k)) exit
      end do
      kPlus = mIdxPlus(k)

      ! pick minus walker to replace
      xi = wMsum*myran()
      do k=1,nMinus
        if (xi <= mWSumMinus(k)) exit
      end do
      kMinus = mIdxMinus(k)
      !!print*,"replace:",kPlus,kMinus
      call replaceWalker(sample,kMinus,kPlus)
    end do

    if (mMode==1) then
      call resettoMean(sample)
    else
      call resetWeights(sample)
    end if

  end subroutine qmc_reconfNew



  subroutine qmc_reconf1(sample)
  !----------------------------!
   
  ! Our own implementation?! 
   
    type(RWSample), intent(inout) :: sample
    type(RandomWalker), pointer   :: rwp => null()
    type(Weighta),allocatable :: recfW(:)

    integer       :: count,count1,count2
    integer       :: Nrn,Nrp
    integer       :: i,j,k
    integer       :: siza,count3,count4

    logical       :: mCheck

    real(r8)        :: wgtsum1,wgtsum2
    real(r8)        :: twgtsum1,twgtsum2
    real(r8)        :: tmp1,tmp2,xi
    real(r8)        :: tmp
    real(r8)        :: allWgt
    real(r8)        :: W
    real(r8),allocatable        :: tmpWgt(:)

    siza = int(getMaxSampleSize(sample))
    allocate(recfW(siza),tmpWgt(siza))
    W = getMeanWeight(sample)

    count  = 0
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    wgtsum1 = 0
    wgtsum2 = 0
    
    rwp => getFirst(sample)
    do
      count = count + 1
      if ((wgt(rwp)/W) < 1 ) then   
        count1 = count1 + 1  
        recfW(count1)%in = count   
        recfW(count1)%wn = dabs((wgt(rwp)/W) - 1.0d0)
        wgtsum1 = wgtsum1 + dabs((wgt(rwp)/W) - 1.0d0)
      else if ((wgt(rwp)/W) >= 1) then
        count2 = count2 + 1
        recfW(count2)%ip = count   
        recfW(count2)%wp = dabs((wgt(rwp)/W) - 1.0d0) 
        wgtsum2 = wgtsum2 + dabs((wgt(rwp)/W) - 1.0d0)           
      end if    
      if (.not.isNext(sample)) exit
      rwp => getNext(sample)
    end do

    recfW%wp = recfW%wp/wgtsum2  

    tmpWgt = 0.0d0
    do i = 1, count2        
      tmp2 = 0.0d0
      do j = 1,i          
        tmp2 = tmp2 + recfW(j)%wp/wgtsum2 
      end do   
      tmpWgt(i) = tmp2    
    end do  
    recfW%wp  = tmpWgt

    do i = 1, count1
      xi  = myran()
      if (xi <= recfW(i)%wn) then
        count3 = count3 + 1 
        xi = myran()
        do j = 1, count2
            if (xi <= recfW(j)%wp) then
            count4 = count4 + 1
            call replaceWalker(sample, recfW(i)%in, recfW(j)%ip)
            exit
          end if
        end do          
      end if
    end do   
       
    if (count3 /= count4) call abortp("qmc_reconf:No.of del. rwps .ue. copied rwps")
    call resettoMean(sample)
    deallocate(recfW,tmpWgt)
   
  end subroutine qmc_reconf1
  
    
  !----------------------------!
  subroutine qmc_reconf2(sample)
  !----------------------------!
  

    type(RWSample), intent(inout) :: sample
    type(RandomWalker), pointer   :: rwp => null()
    type(Weighta),allocatable :: recfW(:)

    integer       :: count,count1,count2
    integer       :: Nrn,Nrp
    integer       :: i,j,k
    integer       :: siza,count3,count4

    logical       :: mCheck

    real(r8)        :: wgtsum1,wgtsum2
    real(r8)        :: twgtsum1,twgtsum2
    real(r8)        :: tmp1,tmp2,xi
    real(r8)        :: tmp
    real(r8)        :: allWgt
    real(r8)        :: W
    real(r8),allocatable        :: tmpWgt(:)


    siza = int(getMaxSampleSize(sample))
    allocate(recfW(siza),tmpWgt(siza))
    W = getMeanWeight(sample)

    count  = 0
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    wgtsum1 = 0.d0
    wgtsum2 = 0.d0

    rwp => getFirst(sample)
    Do
      count = count + 1  
      if ((wgt(rwp)/W) < 1d0 ) then   
         count1 = count1 + 1  
         recfW(count1)%in = count   
         recfW(count1)%wn = dabs((wgt(rwp)/W) - 1.0d0)
         wgtsum1 = wgtsum1 + dabs((wgt(rwp)/W) - 1.0d0)
      elseif ((wgt(rwp)/W) >= 1d0) then
         count2 = count2 + 1
         recfW(count2)%ip = count   
         recfW(count2)%wp = dabs((wgt(rwp)/W) - 1.0d0)
         wgtsum2 = wgtsum2 + dabs((wgt(rwp)/W) - 1.0d0) 
      endif    
      if (.not.isNext(sample)) exit
     rwp => getNext(sample)
    enddo
    
    Nrp = int(wgtsum2 + myran())

    tmpWgt = 0.0d0
    Do i = 1, count1
       tmp1 = 0.0d0      
      Do j = 1,i
        tmp1 = tmp1 + recfW(j)%wn/wgtsum1           
      enddo
      tmpWgt(i) = tmp1           
    enddo  
    recfW%wn  = tmpWgt 
    
     tmpWgt = 0.0d0
     Do i = 1, count2        
        tmp2 = 0.0d0
      Do j = 1,i          
        tmp2 = tmp2 + recfW(j)%wp/wgtsum2 
      enddo   
      tmpWgt(i) = tmp2    
    enddo  
    recfW%wp  = tmpWgt
    

!  Copy Np Walkers randomly and seperately
!-- aber dasselbe Feld soll nicht zweimal gelöscht/kopiert werden
 
if (Nrp>=1) then

    tmp=0
    k=1

    do while(k <= Nrp)
      xi = myran()
      do i = 1, count2
        if (xi <= recfW(i)%wp) then
         if (i /= tmp) then
             count3 = count3 + 1
             rwp => null()
             rwp => getWalker(sample,recfW(i)%ip)
             call appendWalker(sample,rwp)
             k = k+1
             tmp = i
             exit
         else 
            exit
         endif
        endif
      enddo
    enddo  

!   Delete Np Walkers randomly and seperately

   tmp = 0
   k = 1

   do while(k <= Nrp)
     xi = myran()
     do i = 1, count1
       if( xi <= recfW(i)%wn) then
         if(i /= tmp) then
           count4 = count4 + 1 
           call deleteWalker(sample, recfW(i)%in)
           k = k + 1
           tmp = i
           exit
         else
           exit
         endif
       endif
     enddo   
   enddo

endif ! Nrp>=1

if (count3 /= count4) call abortp("qmc_reconf: No. of deleted rwps /= copied rwps")

call resettoMean(sample)
deallocate(recfW, tmpWgt)
   
end subroutine qmc_reconf2


end module reconfg_m
