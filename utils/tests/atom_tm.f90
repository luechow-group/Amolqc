! TODO UNKNOWN ORIGIN AND AUTHOR (OLDER THAN 2011?)
!
! SPDX-License-Identifier: GPL-3.0-or-later

module atom_tm
   use atom_m
   implicit none
contains
   subroutine atom_test
      type(atom) :: a1,a2,a3
      type(atom) :: mol(5) !,mol1(9)
      ! integer i

      call atoms_configure(angstrom=.true.,withSA=.false.,withBasis=.false.,withCharge=.false.)

      call atom_init(a1,1.1d0,2.2d0,3.3d0,symbol='C ')
      a2 = atom(0.0d0,0.1d0,0.2d0,1.0d0,1,0,0,'H ','')
      call atom_init(mol(1),1.0d0,2.0d0,3.0d0,atomIdx=8)
      a3 = a2
      mol(2) = a1
      mol(3) = a2
      mol(4) = a3
      mol(5) = a1
      call atoms_write(mol,6)
      write(*,*)

      ! open(10,file='c2h5oh.mol')
      ! call atoms_readFromFile(mol1,10,.true.)
      ! call atoms_write(mol1,6)
      ! do i=1,size(mol1)
      ! write(*,'(i4)',advance='no') atom_getSA(mol1(i))
      ! enddo
      ! write(*,*)
      ! write(*,*) ' nscenter = ',atoms_getNSCenter(mol1)
      ! write(*,*) ' nclast = ',atoms_getNCLast(mol1)
      ! call atoms_ignoreHAtoms()
      ! write(*,*) ' nscenter = ',atoms_getNSCenter(mol1)
      ! write(*,*) ' nclast = ',atoms_getNCLast(mol1)
      ! write(*,*) ' ne = ',atoms_countElectrons(mol1)
      ! write(*,*) getPSEIdx('O ')
   end subroutine atom_test
end module atom_tm