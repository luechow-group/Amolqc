! Copyright (C) 2019-2020 Vladimir Terzi

module rhoGrid_m
   use kinds_m, only: r8
   use parsing_m, only: getinta, getdbla, getintarra, getdblarra, finda
   use globalUtils_m, only: iul
   use global_m, only: getNElec
   use error_m, only: asserts, assert
   use wfData_m, only: getNNuc, atoms
   use rhoData_m, only: rhoData_t
   use rwSample_m, only: rwSample, getSampleSize, getFirst, getNext
   use randomWalker_m, only: RandomWalker, pos
   use atom_m, only: atom_getPosition
   use bins3D_m, only: Bins3D_t, GetBins3D
   use electronDensity_m, only: rho

   implicit none

   private
   public :: RhoGrid_t, InitRhoGrid

   real(r8), parameter :: defaultBinSize = .1_r8, defaultOffset = 5, defaultSearchRadius = 3

   type RhoGrid_t
      private
      integer, allocatable :: atomsIdx(:)
      real(r8) :: offset(3)
      integer :: searchRadius
      type(Bins3D_t) :: bins
      type(rhoData_t) :: rhoData
   contains
      procedure :: IsInitialized => RhoGrid_IsInitialized
      procedure :: GetCellSize => RhoGrid_GetCellSize
      procedure :: GetOffset => RhoGrid_GetOffset
      procedure :: AtomInGrid => RhoGrid_AtomInGrid
      procedure :: InBounds => RhoGrid_InBounds
      procedure :: AddSample => RhoGrid_AddSample
      procedure :: AddElecs => RhoGrid_AddElecs
      procedure :: GetValue => RhoGrid_GetValue
      procedure :: GetBasin => RhoGrid_GetBasin
      procedure :: CompareRho => RhoGrid_CompareRho
   end type RhoGrid_t

contains
   subroutine InitRhoGrid(lines, nl, smpl, rhoGrid)
      integer, intent(in) :: nl
      character(len = 120), intent(in) :: lines(nl)
      type(rwSample), intent(inout) :: smpl
      type(RhoGrid_t), intent(out) :: rhoGrid

      real(r8) :: binSize, offset(3), assignThresh
      integer, allocatable :: atomsIdx(:)
      real(r8), allocatable :: buffer(:)
      integer :: searchRadius, i, iflag

      call getdbla(lines, nl, 'bin_size=', binSize, iflag)
      if (iflag /= 0) binSize = defaultBinSize
      call getdblarra(lines, nl, 'offset=', buffer, iflag)
      if (iflag == 0) then
         if (SIZE(buffer) == 3) then
            offset = buffer
         else if (SIZE(buffer) == 1) then
            offset = buffer(1)
         else
            if (asserts) call assert(.false., 'InitRhoGrid: the offset should be either a 3D vector or a scalar.')
         end if
         if (asserts) call assert(ALL(offset >= 0._r8), 'InitRhoGrid: the offset should be positive.')
      else
         offset = defaultOffset
      end if
      call getintarra(lines, nl, 'atoms=', atomsIdx, iflag)
      if (iflag /= 0) atomsIdx = [(i, i = 1, SIZE(atoms))]
      call getinta(lines, nl, 'search_radius=', searchRadius, iflag)
      if (iflag /= 0) searchRadius = defaultSearchRadius
      if (asserts) call assert(searchRadius > 0, 'InitRhoGrid: the search radius should be positive.')
      assignThresh = binSize
      call getdbla(lines, nl, 'assign_thresh=', assignThresh, iflag)
      if (asserts) call assert(assignThresh > 0._r8, 'InitRhoGrid: the assign threshold should be positive.')
      rhoGrid = GetRhoGrid(binSize, offset, atomsIdx, smpl, searchRadius, assignThresh)
      if (finda(lines, nl, 'compare')) call rhoGrid%CompareRho()
   end subroutine InitRhoGrid


   function GetRhoGrid(binSize, offset, atomsIdx, smpl, searchRadius, assignThresh) result(rhoGrid)
      real(r8), intent(in) :: binSize, offset(3)
      integer, intent(in) :: atomsIdx(:)
      type(rwSample), intent(inout) :: smpl
      integer, intent(in) :: searchRadius
      real(r8), intent(in) :: assignThresh
      type(RhoGrid_t) :: rhoGrid

      real(r8) :: bounds(3, 2)
      real(r8), allocatable :: positions(:, :)
      integer :: i, j

      rhoGrid%offset = offset
      allocate(rhoGrid%atomsIdx(SIZE(atomsIdx)))
      rhoGrid%atomsIdx = atomsIdx
      allocate(positions(3, SIZE(atomsIdx)))
      do i = 1, SIZE(atomsIdx)
         positions(:, i) = atom_getPosition(atoms(atomsIdx(i)))
      end do
      bounds(:, 1) = positions(:, 1)
      bounds(:, 2) = positions(:, 1)
      do i = 2, SIZE(positions, 2)
         do j = 1, 3
            if (positions(j, i) < bounds(j, 1)) bounds(j, 1) = positions(j, i)
            if (positions(j, i) > bounds(j, 2)) bounds(j, 2) = positions(j, i)
         end do
      end do
      bounds(:, 1) = bounds(:, 1) - offset
      bounds(:, 2) = bounds(:, 2) + offset
      rhoGrid%bins = GetBins3D(bounds, binSize)

      rhoGrid%searchRadius = searchRadius
      call rhoGrid%rhoData%init(assignThresh, 0._r8)
      call rhoGrid%AddSample(smpl)
   end function GetRhoGrid


   function RhoGrid_IsInitialized(this) result(val)
      class(RhoGrid_t), intent(in) :: this
      logical :: val

      val = this%bins%IsInitialized()
   end function RhoGrid_IsInitialized


   function RhoGrid_GetCellSize(this) result(val)
      class(RhoGrid_t), intent(in) :: this
      real(r8) :: val

      val = this%bins%GetBinSize()
   end function RhoGrid_GetCellSize


   function RhoGrid_GetOffset(this) result(val)
      class(RhoGrid_t), intent(in) :: this
      real(r8) :: val(3)

      val = this%offset
   end function RhoGrid_GetOffset


   function RhoGrid_AtomInGrid(this, atomIdx) result(val)
      class(RhoGrid_t), intent(in) :: this
      integer, intent(in) :: atomIdx
      logical :: val

      val = ANY(this%atomsIdx == atomIdx)
   end function RhoGrid_AtomInGrid


   function RhoGrid_InBounds(this, coords) result(val)
      class(RhoGrid_t), intent(in) :: this
      real(r8), intent(in) :: coords(3)
      logical :: val

      val = this%bins%InBounds(coords(1), coords(2), coords(3))
   end function RhoGrid_InBounds


   subroutine RhoGrid_AddSample(this, smpl)
      class(RhoGrid_t), intent(inout) :: this
      type(rwSample), intent(inout) :: smpl

      integer :: smplSize
      real(r8) :: smplPoint(3 * getNElec())
      type(RandomWalker), pointer :: rwp
      integer :: x, i, iflag

      smplSize = getSampleSize(smpl)
      rwp => getFirst(smpl)
      do x = 1, smplSize
         call pos(rwp, smplPoint)
         do i = 1, 3 * getNElec(), 3
            call this%bins%IncValue(smplPoint(i), smplPoint(i + 1), smplPoint(i + 2), iflag)
         end do
         rwp => getNext(smpl)
      end do
   end subroutine RhoGrid_AddSample


   subroutine RhoGrid_AddElecs(this, elecs)
      class(RhoGrid_t), intent(inout) :: this
      real(r8), intent(in) :: elecs(:, :)

      integer :: i, iflag

      do i = 1, SIZE(elecs, 2)
         call this%bins%IncValue(elecs(1, i), elecs(2, i), elecs(3, i), iflag)
      end do
   end subroutine RhoGrid_AddElecs


   function RhoGrid_GetValue(this, coords) result(val)
      class(RhoGrid_t), intent(in) :: this
      real(r8), intent(in) :: coords(3)
      real(r8) :: val

      val = 1._r8 * this%bins%GetValue(coords(1), coords(2), coords(3)) / this%bins%GetTotNum()
   end function RhoGrid_GetValue


   subroutine RhoGrid_GetBasin(this, coords, basin, coordsMax, inBounds)
      class(RhoGrid_t), intent(in) :: this
      real(r8), intent(in) :: coords(3)
      integer, intent(out) :: basin
      real(r8), intent(out), optional :: coordsMax(3)
      logical, intent(out), optional :: inBounds

      integer :: iflag

      basin = 0
      inBounds = .false.
      call this%bins%FindMax(coords, this%searchRadius, coordsMax, iflag)
      if (iflag == 0) then
         basin = this%rhoData%getBasin(coordsMax)
         inBounds = .true.
      end if
   end subroutine RhoGrid_GetBasin


   subroutine RhoGrid_CompareRho(this)
      class(RhoGrid_t), intent(in) :: this

      integer :: gridShape(3), totNum, numCells
      real(r8) :: msd  ! mean squared deviation
      real(r8) :: cellVol, diff, maxDiff, maxValRef, maxValGrid
      real(r8) :: r(3), valRef, valGrid
      integer :: i, j, k

      gridShape = this%bins%GetShape()
      numCells = PRODUCT(gridShape)
      totNum = this%bins%GetTotNum()
      cellVol = this%bins%GetBinSize() ** 3
      msd = 0
      maxDiff = 0

      maxValRef = 0
      maxValGrid = 0
      do k = 1, gridShape(3)
         do j = 1, gridShape(2)
            do i = 1, gridShape(1)
               call this%bins%CellToCoords(i, j, k, r(1), r(2), r(3))
               call rho(r, valRef)
               valRef = valRef * cellVol
               valGrid = 1._r8 * this%bins%GetValue(i, j, k) / totNum
               if (maxValRef < valRef) maxValRef = valRef
               if (maxValGrid < valGrid) maxValGrid = valGrid
               diff = ABS(valGrid - valRef)
               if (maxDiff < diff) maxDiff = diff
               msd = msd + diff ** 2
            end do
         end do
      end do
      msd = msd / numCells
      write(iul, '(a, f0.16)') 'root mean squared deviation (RMSD) of electronic density from grid = ', SQRT(msd)
      write(iul, '(a, f0.16)') 'maximal deviation = ', maxDiff
      write(iul, '(a, f0.16)') 'maximal calculated electronic density in a cell = ', maxValRef
      write(iul, '(a, f0.16)') 'maximal grid electronic density in a cell = ', maxValGrid
   end subroutine RhoGrid_CompareRho
end module rhoGrid_m
