! Copyright (C) 1997-1999, 2012 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

! collection of machine-dependent utility routines
module machine_m
#ifdef NAG
  use, intrinsic :: f90_unix_env, only: getenv, iargc, getarg
#endif
  use global_m
  use error_m, only: error
  use mpiInterface_m, only: myMPIBcastInteger, myMPIBcastString

  implicit none
contains

!---------------------------!
subroutine myGetArg(n,s,ierr)
!---------------------------!

! master reads arg broadcasts to other procs 
  integer, intent(in)              :: n    ! get n-th command line argument
  character(len=*), intent(out)    :: s    ! return argument
  integer, intent(out)             :: ierr ! error code: 1=no argument


  ierr = 0
  if (mytid==0) then
     if (iargc() < n) then
        ierr = 1
     else
        call getarg(n,s)
     endif
  endif
  call myMPIBcastInteger(ierr,1)
  if (ierr == 0) call myMPIBcastString(s,len(s))

end subroutine myGetArg
  
!--------------------------------!
subroutine myGetArgLocal(n,s,ierr)
!--------------------------------!

  integer, intent(in)              :: n    ! get n-th command line argument
  character(len=*), intent(out)    :: s    ! return argument
  integer, intent(out)             :: ierr ! error code: 1=no argument

  ierr = 0
  if (iargc() < n) then
     ierr = 1
  else
     call getarg(n,s)
  endif

end subroutine myGetArgLocal
  

!--------------------------
subroutine myGetHost(hostn)
!--------------------------
  ! returns the hostname in 'hostn'
  character(len=*) :: hostn

  call getenv('HOSTNAME',hostn)
  !call hostnm(hostn)    ! SGI, Alpha
  !ihost = hostnm_(hostn)   ! AIX
end subroutine myGetHost


!-----------------------------
subroutine myGetEnv(name,value)
!-----------------------------
  ! returns content of environment variable name
  character(len=*) name, value

  if (mytid==0) then
    call getenv(name,value)
    if (value == "") then
      call error('environment variable '//trim(name)//' is not set/empty.')
    end if
  end if
  call myMPIBcastString(value,len(value))

end subroutine myGetEnv

!--------------------------
subroutine myGetDate(full_date)
!--------------------------
  ! returns current time and date in 'date'
  ! follows ISO8601 : https://www.w3.org/TR/NOTE-datetime

  character(len=*), intent(out) :: full_date
  character(len=8)              :: date
  character(len=10)             :: time
  character(len=5)              :: zone
  character(len=2)              :: hour, minute, month, day, zone2
  character(len=3)              :: zone1
  character(len=4)              :: year
  character(len=6)              :: second

  call date_and_time(date=date, time=time, zone=zone)

  year  = date(1:4)
  month = date(5:6)
  day   = date(7:8)
  hour   = time(1:2)
  minute = time(3:4)
  second = time(5:10)
  zone1 = zone(1:3)
  zone2 = zone(4:5)

  full_date = year//'-'//month//'-'//day//'T'//hour//':'//minute//':'//second//zone1//':'//zone2
  !call fdate_(full_date)     ! AIX
end subroutine myGetDate

end module machine_m
