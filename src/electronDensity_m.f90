! Copyright (C) 2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module electronDensity_m
   use kinds_m, only: r8
   use parsing_m, only: getdbla, getinta
   use global_m, only: iul, getNElec
   use eConfigs_m, only: eConfigArray, eConfigArray_new, eConfigArray_set, eConfigArray_destroy
   use error_m, only: debug, abortp
   use aos_m, only: aocalc, ao1calc
   use wfData_m, only: getNNuc, atoms
   use error_m, only: assert
   use mos_m, only: mocalc, mo1calc, mat, mat1x, mat1y, mat1z
   use atom_m, only: atom_getPosition, atom_getDistance
   use fctn_module, only: Function_t
   use wfData_m, only: mult

   implicit none
   private
   public :: calculateDensity, scanBondDensity, fctn_rho, fctn_log_rho, rho

   type, extends(Function_t) :: fctn_rho
   contains
      procedure :: eval_fg => rho_fg
      procedure :: eval => rho_f
   end type fctn_rho

   type, extends(Function_t) :: fctn_log_rho
   contains
      procedure :: eval_fg => log_rho_fg
      procedure :: eval => log_rho_f
   end type fctn_log_rho


contains
   subroutine scanBondDensity(lines, nl, verboseInput)
      ! atm only calculates the density path for all orbs
      integer, intent(in) :: nl
      character(len=120), intent(in) :: lines(nl)
      logical, intent(in), optional :: verboseInput
      logical :: verbose = .true.
      integer :: A, B
      integer :: steps
      real(r8) :: stepSize = 0._r8
      real(r8) :: x(3)
      integer :: iflag, i
      real(r8) :: value = 0._r8
      real(r8) :: gradient(3) = 0._r8
      real(r8) :: bondVector(3) = 0._r8

      if (PRESENT(verboseInput)) verbose = verboseInput

      call getinta(lines,nl,'A=',A,iflag)
      call getinta(lines,nl,'B=',B,iflag)
      call getinta(lines,nl,'steps=',steps,iflag)

      write(iul,'(A,2i3)') "Scanning bond density between atoms", A, B

      call assert(A > 0 .and. B > 0 .and. A <= SIZE(atoms) .and. B <= SIZE(atoms), 'scanBondDensity: invalid atom index')
      write(iul,*) "Steps: ", steps
      write(iul,'(A)') ""

      stepSize = atom_getDistance(atoms(A), atoms(B)) / (steps + 1)

      write(iul,'(A)') "Rho  ScalarProduct Gradient"

      bondVector = (atom_getPosition(atoms(B)) - atom_getPosition(atoms(A))) / atom_getDistance(atoms(A), atoms(B))

      x = atom_getPosition(atoms(A))
      do i = 1, steps
         x = x + stepSize * bondVector
         call rho(x, value, gradient)
         write(iul,*) value, DOT_PRODUCT(gradient, atom_getPosition(atoms(B)) - atom_getPosition(atoms(A))), NORM2(gradient)
      end do

   end subroutine scanBondDensity

   subroutine calculateDensity(lines, nl)
      ! atm only calculates the density for all orbs
      integer, intent(in) :: nl
      character(len=120), intent(in) :: lines(nl)
      real(r8) :: r(3) = 0._r8
      integer :: iflag
      real(r8) :: value = 0._r8
      real(r8) :: gradient(3) = 0._r8

      call getdbla(lines,nl,'x=',r(1),iflag)
      call getdbla(lines,nl,'y=',r(2),iflag)
      call getdbla(lines,nl,'z=',r(3),iflag)

      write(iul,'(A,3f10.5)') "Calculating determinants density at x,y,z = ", r(1), r(2), r(3)

      call rho(r, value, gradient)

      write(iul,'(A,f10.5)') "Density is ", value
      write(iul,'(A,3f10.5)') "Gradient is ", gradient
   end subroutine calculateDensity

   function rho_f(this, x) result(f)
      class(fctn_rho), intent(in)  :: this
      real(r8), intent(in) :: x(:)
      real(r8) :: f

      call assert(SIZE(x) == 3, "rho_f: illegal size")

      call rho(x, f)
      f = -1 * f
   end function rho_f

   subroutine rho_fg(this, x, f, g, mask)
      class(fctn_rho), intent(in)  :: this
      real(r8), intent(in) :: x(:)
      real(r8), intent(out) :: f
      real(r8), intent(out) :: g(:)
      logical, intent(in), optional :: mask(SIZE(x))

      call assert(SIZE(x) == 3 .and. SIZE(g) == 3, "rho_fg: illegal size")

      call rho(x, f, g)
      f = -1 * f
      g = -1 * g
      if (PRESENT(mask)) where (mask) g = 0._r8
   end subroutine rho_fg

   function log_rho_f(this, x) result(f)
      class(fctn_log_rho), intent(in)  :: this
      real(r8), intent(in) :: x(:)
      real(r8) :: f

      call assert(SIZE(x) == 3, "rho_f: illegal size")

      call rho(x, f)
      f = -LOG(f)
   end function log_rho_f

   subroutine log_rho_fg(this, x, f, g, mask)
      class(fctn_log_rho), intent(in)  :: this
      real(r8), intent(in) :: x(:)
      real(r8), intent(out) :: f
      real(r8), intent(out) :: g(:)
      logical, intent(in), optional :: mask(SIZE(x))

      call assert(SIZE(x) == 3 .and. SIZE(g) == 3, "rho_fg: illegal size")

      call rho(x, f, g)
      g = -1._r8/f * g
      f = -LOG(f)
      if (PRESENT(mask)) where (mask) g = 0._r8
   end subroutine log_rho_fg

   subroutine rho(position, value, gradient, orbitals)
      real(r8), intent(in) :: position(3)
      real(r8), intent(out) :: value
      real(r8), optional, intent(out) :: gradient(3)
      integer, optional, intent(in) :: orbitals
      real(r8) :: rrai(getNNuc(), getNElec(), 1)
      type(eConfigArray)  :: eca
      integer :: closed_orbitals, a

      real(r8) :: x(getNElec()), y(getNElec()), z(getNElec())

      if (PRESENT(orbitals)) then
         closed_orbitals = orbitals
      else
         call assert(mult == 1, "rho: implementation atm only works for singlet systems")
         closed_orbitals = getNElec() / 2
      end if

      ! filling rrai (electron-nucleus distances) only for electron 1
      rrai = 0._r8
      do a = 1, getNNuc()
         rrai(a, 1, 1) = NORM2(position - atom_getPosition(atoms(a)))
      enddo

      ! filling x, y, and z only for electron 1
      x = 0._r8; y = 0._r8; z = 0._r8;
      x(1) = position(1)
      y(1) = position(2)
      z(1) = position(3)

      ! calculating aos and mos for electron 1
      if (PRESENT(gradient)) then
         call eConfigArray_new(eca, getNElec(), 1)
         call eConfigArray_set(eca, 1, x, y, z)

         call aocalc(1, eca, rrai)
         call mocalc(1) ! calculates mat, mat1x, mat1y, mat1z

         gradient(1) = SUM( 2 * mat(1:closed_orbitals,1,1) * mat1x(1:closed_orbitals,1,1) )
         gradient(2) = SUM( 2 * mat(1:closed_orbitals,1,1) * mat1y(1:closed_orbitals,1,1) )
         gradient(3) = SUM( 2 * mat(1:closed_orbitals,1,1) * mat1z(1:closed_orbitals,1,1) )

         gradient = 2 * gradient
         call eConfigArray_destroy(eca)
      else
         call ao1calc(1, x, y, z, rrai(:,:,1))
         call mo1calc(1) ! calculates mat
      end if

      value = 2 * SUM(mat(1:closed_orbitals,1,1) ** 2)
   end subroutine rho
end module electronDensity_m
