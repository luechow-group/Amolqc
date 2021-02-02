! Copyright (C) 2013 Kaveh Haghighi Mood
!
! SPDX-License-Identifier: GPL-3.0-or-later

program cusptest
    use kinds_m, only: r8
	use aosData_m, only: aoinputabs, aoinputex
	use global_m
	use init_m, only: initAmolqc
	use jastrow_m
	use random_m, only: init_ran, myran
	use wfData_m
	use eConfigs_m
	use eloc_m
	implicit none

	character(len=10) :: jt,tmp
	character(len=120) :: lines(1)
	integer :: argc, iflag, optmode = 1
	real(r8) :: dummy

	argc = iargc()

	if (argc < 1) then
		print*, "jastype required"
		stop
	else if (argc > 1) then
		call getarg(2, tmp)
		read(tmp,*) optmode
	endif
	call getarg(1, jt)

	print*,"Testing cusp for ", jt, " with optmode ", optmode

	! rng
	call initAmolqc("cusptest")
	dummy = init_ran(4) ! http://xkcd.com/221/

	! set up atoms
	ncenter = 1
	allocate(atoms(1))
	call atom_init(atoms(1),0d0,0d0,0d0,symbol="He")
	ne = atoms_countElectrons(atoms)
	nalpha = ne / 2
	nbeta = nalpha
  nscenter = atoms_getNSCenter(atoms)
  

	! set up basis
	basis = "TZPAE"
	call aoinputabs()

	

	!call jasCalcWDerivs
	call calcJastrow(ncenter, ne, optmode)
contains
	subroutine calcJastrow(ncenter, ne, optmode)
		integer, intent(in) :: ncenter, ne, optmode
		real(r8) :: x(ne), y(ne), z(ne)
		real(r8) :: rai(ncenter, ne), rij(ne, ne)
		real(r8) :: u, ugrad(3*ne), ulapl, ulapli(ne)
		integer :: a, i, j, np, iu = 0
		real(r8), allocatable :: params(:)
		type(eConfigArray) :: ec
 		real(r8) :: rrai(ncenter,ne,1)

		do i = 1, ne
			x(i) = i*myran()*1d-20; y(i) = 0d0; z(i) = 0d0
		enddo
  	
  	call eConfigArray_new(ec,ne,1)
  	call eConfigArray_set(ec,1,x,y,z)
  	call eloc_initialize(1)
  	call aos_initialize(1)
  	norb = 1
  	call mos_initialize(1)
  	call aocalc(0, ec, rrai)

  	! set up jastrow
		lines(1) = "$change_jastrow(new_jastrow=" // trim(jt) // ")"
		call jasChangeType(lines, 1)
		
		do i = 1, ne
			do j = i +1, ne
				rij(i,j) = x(i) - x(j)
			enddo
			do a = 1, ncenter
				rai(a,i) = x(i)
			enddo
		enddo

		np = getJastrowParamCnt(optmode)
		allocate(params(np))
		do i = 1, np
			params(i) = 5 * myran()
		enddo
		call putJastrowParamVector(optmode, params)
		print*,"Set Jastrow parameters"
		print*,"======================"
		call jasoutput(iu)
		print*,"======================"

		call jasCalcWDerivs(0,x,y,z,rai,rij,"none",u,ugrad,ulapl,ulapli,1)
		print*,u,0.5d0*rij(1,2)
	end subroutine
end program
