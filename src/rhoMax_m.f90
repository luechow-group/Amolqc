! Copyright (C) 2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module rhoMax_m
   use kinds_m, only: r8
   use parsing_m, only: getinta, getdbla, getintarra, finda
   use error_m, only: asserts, assert
   use fctn_module, only: Function_t
   use electronDensity_m, only: fctn_rho, fctn_log_rho
   use minimizer_w_sing_module, only: minimizer_w_sing
   use minimizer_ws_factory_module, only: create_ws_minimizer
   use randomWalker_m, only: RandomWalker, pos, resetTo
   use wfData_m, only: atoms, getNNuc
   use global_m, only: iul, iull, getNElec, logmode
   use rhoData_m, only: rhoData_t
   use atom_m, only: atom_getPosition, atoms_getPositionMatrix
   use wfData_m, only: atoms

   implicit none
   private
   public :: rhoMax_t

   type rhoMax_t
      private
      logical :: initialized = .false.
      class(minimizer_w_sing), pointer :: minimizer_p => null()
      type(rhoData_t) :: rhoData
      class(Function_t), allocatable :: fg
      integer :: verbose_ = 0
      logical :: use_log = .false.
   contains
      procedure :: init => rhoMax_init
      procedure :: destroy => rhoMax_destroy
      procedure :: writeParams => rhoMax_writeParams
      procedure :: opt => rhoMax_opt
      procedure :: writeResults => rhoMax_writeResults
      procedure :: isConverged => rhoMax_isConverged
      procedure :: isInitialized => rhoMax_isInitialized
      procedure :: getBasin => rhoMax_getBasin
   end type rhoMax_t

contains
   subroutine rhoMax_init(this, lines, nl)
      class(rhoMax_t), intent(inout) :: this
      character(len=120), intent(in) :: lines(:)
      integer, intent(in) :: nl
      real(r8) :: assignThresh, printThresh, preAssignThresh
      integer, allocatable :: fragments(:)
      integer :: i, iflag

      this%minimizer_p => create_ws_minimizer(lines)

      call getinta(lines, nl, "verbose=", this%verbose_, iflag)
      if (iflag /= 0) this%verbose_ = 0

      call this%minimizer_p%set_verbose(this%verbose_)
      call this%minimizer_p%set_verbose_unit(iull)

      call this%minimizer_p%set_singularities(atoms_getPositionMatrix(atoms), &
                                              scalings=REAL(atoms%Get_atomic_number(), r8))

      call getdbla(lines, nl, "assign_thresh=", assignThresh, iflag)
      if (iflag /= 0) assignThresh = 0.1

      call getdbla(lines, nl, "print_thresh=", printThresh, iflag)
      if (iflag /= 0) printThresh = 0.05

      call getintarra(lines, nl, "fragments=", fragments, iflag)
      if (iflag /= 0) fragments = [(i, i=1, getNNuc())]

      if (finda(lines, nl, "use_log")) this%use_log = .true.

      call getdbla(lines, nl, "assign_pre_thresh=", preAssignThresh, iflag)
      if (iflag /= 0) preAssignThresh = 0.01

      if (asserts) call assert(SIZE(fragments) == getNNuc(), 'rhoMax_init: invalid fragments array size.')

      if (this%use_log) then
         allocate(this%fg, source=fctn_log_rho())
      else
         allocate(this%fg, source=fctn_rho())
      end if

      this%rhoData = rhoData_t(assignThresh, printThresh, fragments, preAssignThresh)

      this%initialized = .true.

      if (logmode >= 2) then
         call this%writeParams(iul)
      end if
   end subroutine rhoMax_init

   subroutine rhoMax_destroy(this)
      class(rhoMax_t), intent(inout)    :: this

      nullify(this%minimizer_p)
      deallocate(this%fg)
      this%verbose_ = 0
      this%initialized = .false.
   end subroutine rhoMax_destroy

   subroutine rhoMax_writeParams(this, iu)
      class(rhoMax_t), intent(in) :: this
      integer, intent(in) :: iu

      call this%minimizer_p%write_params(iu)
   end subroutine rhoMax_writeParams

   function rhoMax_isConverged(this) result(converged)
      class(rhoMax_t), intent(in) :: this
      logical :: converged

      converged = this%minimizer_p%is_converged()
   end function rhoMax_isConverged

   function rhoMax_isInitialized(this) result(initialized)
      class(rhoMax_t), intent(in) :: this
      logical :: initialized

      initialized = this%initialized
   end function rhoMax_isInitialized

   subroutine rhoMax_opt(this, rw, update_walker)
      class(rhoMax_t), intent(inout)      :: this
      type(RandomWalker), intent(inout) :: rw
      logical, intent(in)               :: update_walker
      real(r8)  :: r(3*getNElec()), rMax(3*getNElec())
      integer :: i, basins(getNElec())
      logical :: allConverged, converged

      call pos(rw, r)

      allConverged = .true.
      do i = 1, getNElec()
         call this%getBasin(r(3*i-2:3*i), basins(i), rMax=rMax(3*i-2:3*i), converged=converged)
         allConverged = allConverged .and. converged
      end do

      if (allConverged) then
         call this%rhoData%addData(basins)
      end if
   end subroutine rhoMax_opt

   subroutine rhoMax_getBasin(this, r, basin, rMax, converged)
      class(rhoMax_t), intent(inout)      :: this
      real(r8), intent(in)  :: r(3)
      integer, intent(out) :: basin
      real(r8), optional, intent(out) :: rMax(3)
      logical, optional, intent(out) :: converged

      real(r8) :: rTemp(3)

      basin = 0

      if (this%verbose_ >= 2) then
         write(iull,*) 'before minimize:'
         write(iull,'(3f15.6)') r
      end if

      basin = this%rhoData%getBasin(r, preAssign=.true.)
      if (basin /= 0) then  ! electron is already within %rhoData%assignThresh of a nuc
         rTemp = atom_getPosition(atoms(basin))
         converged = .true.
      else
         basin = 0

         rTemp = r
         call this%minimizer_p%reset()
         call this%minimizer_p%minimize(this%fg, rTemp)

         converged = this%minimizer_p%is_converged()

         if (converged) then
            basin = this%rhoData%getBasin(rTemp)
         end if
      end if

      if (this%verbose_ >= 2) then
         write(iull,*) 'after minimize:'
         write(iull,'(3f15.6)') rTemp
      end if

      if (PRESENT(rMax)) then
         rMax = rTemp
      end if
   end subroutine rhoMax_getBasin

   subroutine rhoMax_writeResults(this, iul)
      class(rhoMax_t), intent(inout) :: this
      integer, intent(in) :: iul

      call this%rhoData%writeResults(iul)
   end subroutine rhoMax_writeResults
end module rhoMax_m
