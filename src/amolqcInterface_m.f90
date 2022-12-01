! Copyright (C) 2017 Arne Luechow
! Copyright (C) 2017 Michael A. Heuer
!
! SPDX-License-Identifier: GPL-3.0-or-later

! amolqcInterface_m provides a C interface to amolqc

module amolqcInterface_m
   use kinds_m, only: r8
   use iso_c_binding!, only: c_double, c_int, c_ptr, c_char, c_f_pointer
   use global_m
   use init_m
   use elocData_m
   use eloc_m, only: wf_init, eloc
   use eConfigs_m
   use initialPositions_m

   implicit none

   ! enum is not really used in this example
   enum, bind(c)
      enumerator :: GAUSSIAN, DENSITY, LMO
   end enum

contains

   function c_to_f_string(s) result(str)
      character(kind=c_char,len=1), intent(in) :: s(*)
      character(len=:), allocatable :: str
      integer i, nchars
      i = 1
      do
         if (s(i) == c_null_char) exit
         i = i + 1
      end do
      nchars = i - 1  ! Exclude null character from Fortran string
      allocate(character(len=nchars) :: str)
      str = transfer(s(1:nchars), str)
   end function c_to_f_string

   subroutine amolqc_init() bind(c)
      integer seed
      real(r8) :: dummy

      call initAmolqc()

      ! init rng and logmode (replacing init.f90:initGen)
      seed = 101 + mytid
      dummy = init_ran(seed)
      ! try supressing output
      logmode = 0
   end subroutine amolqc_init

   subroutine amolqc_set_wf(nelecs, natoms, filename) bind(c, name='amolqc_set_wf')
      integer(c_int) :: nelecs   ! # of elecs
      integer(c_int) :: natoms   ! # of atoms
      character(kind=c_char), dimension(*), intent(in) :: filename
      character(len=80) :: lines(1)
      integer :: nl

      nl = 1
      lines(1) = "$wf(read, file='"//trim(c_to_f_string(filename))//"')"
      
      print*, lines(1)
      call wf_init(lines,nl)
      !call setNumberOfElectrons(ne)    ! for randomwalker initialization
      !call setNumberOfCenters(ncenter)
      nelecs = ne
      natoms = ncenter
   end subroutine amolqc_set_wf

   subroutine amolqc_initial_positions(mode, n, x) bind(c)
      integer(c_int), value :: mode   ! enum above
      integer(c_int), value :: n   
      real(c_double)  :: x(3*n)
      real(r8)  xx(n), yy(n), zz(n)
      integer i

      if (n /= ne) call abortp("amolqc_initial_positions: n /= ne")

      do i = 1, n
         xx(i) = x(3*i-2)
         yy(i) = x(3*i-1)
         zz(i) = x(3*i)
      enddo

      call createRandomElectronPositions(mode,xx,yy,zz)

      do i = 1, n
         x(3*i-2) = xx(i) 
         x(3*i-1) = yy(i) 
         x(3*i)   = zz(i) 
      enddo
   end subroutine amolqc_initial_positions

   subroutine amolqc_eloc(x, n, phi, u, drift, elocal) bind(c)
      integer(c_int), value :: n
      real(c_double)        :: x(3*n)
      real(c_double)        :: phi
      real(c_double)        :: u
      real(c_double)        :: elocal
      real(c_double)        :: drift(3*n)
      real(r8)  xx(n), yy(n), zz(n)
      integer i
      type(eConfigArray) :: ec

      if (n /= ne) call abortp("amolqc_eloc: n /= ne")

      do i = 1, n
         xx(i) = x(3*i-2)
         yy(i) = x(3*i-1)
         zz(i) = x(3*i)
      enddo

      call eConfigArray_new(ec,n,1)
      call eConfigArray_set(ec,1,xx,yy,zz)
      call eloc(0,ec,'none')

      phi = elPhi(1)
      u = elU(1)
      elocal = elEloc(1)

      do i = 1, n
         drift(3*i-2) = elxDrift(i, 1)
         drift(3*i-1) = elyDrift(i, 1)
         drift(3*i)   = elzDrift(i, 1)
      enddo

      call eConfigArray_destroy(ec)
   end subroutine amolqc_eloc

end module amolqcInterface_m
