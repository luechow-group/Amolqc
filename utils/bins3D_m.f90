! Copyright (C) 2012 Arne Luechow
! Copyright (C) 2019-2020 Vladimir Terzi
!
! SPDX-License-Identifier: GPL-3.0-or-later


module bins3D_m
   use kinds_m, only: r8
   use error_m, only: error
#ifdef MPI
   use MPI_F08, only: MPI_ALLREDUCE, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD
#endif

   implicit none

   private
   public :: Bins3D_t, GetBins3D

   type Bins3D_t
      private
      integer, allocatable :: data(:, :, :)
      real(r8) :: binSize, offset(3)
      integer :: num = 0
      integer :: totNum = 0
   contains
      procedure :: IsInitialized => Bins3D_IsInitialized
      procedure :: GetData => Bins3D_GetData
      procedure :: GetNum => Bins3D_GetNum
      procedure :: GetTotNum => Bins3D_GetTotNum
      procedure :: GetBinSize => Bins3D_GetBinSize
      procedure :: GetOffset => Bins3D_GetOffset
      procedure :: GetShape => Bins3D_GetShape
      procedure, private :: Bins3D_InBounds3I, Bins3D_InBounds3R
      generic :: InBounds => Bins3D_InBounds3I, Bins3D_InBounds3R
      procedure :: CoordsToCell => Bins3D_CoordsToCell
      procedure :: CellToCoords => Bins3D_CellToCoords
      procedure, private :: Bins3D_GetValue3I, Bins3D_GetValue3R
      generic :: GetValue => Bins3D_GetValue3I, Bins3D_GetValue3R
      procedure :: IncValue => Bins3D_IncValue
      procedure :: FindMax => Bins3D_FindMax
      procedure :: AllNodes => Bins3D_AllNodes
   end type Bins3D_t

contains
   function GetBins3D(bounds, binSize) result(bins)
      real(r8), intent(in) :: bounds(3, 2), binSize
      type(Bins3D_t) :: bins

      integer :: shape(3), stat

      if (ANY(bounds(:, 2) <= bounds(:, 1))) call error('GetBins3D: invalid bounds')
      shape = CEILING((bounds(:, 2) - bounds(:, 1)) / binSize)
      bins%binSize = binSize
      bins%offset = .5_r8 * (bounds(:, 2) + bounds(:, 1) - shape * binSize)
      allocate(bins%data(shape(1), shape(2), shape(3)), stat = stat)
      if (stat /= 0) call error('GetBins3D: allocation error')
      bins%data = 0
   end function GetBins3D


   function Bins3D_IsInitialized(this) result(val)
      class(Bins3D_t), intent(in) :: this
      logical :: val

      val = ALLOCATED(this%data)
   end function Bins3D_IsInitialized


   function Bins3D_GetData(this) result(val)
      class(Bins3D_t), intent(in) :: this
      integer :: val(SIZE(this%data, 1), SIZE(this%data, 2), SIZE(this%data, 3))

      val = this%data
   end function Bins3D_GetData


   function Bins3D_GetNum(this) result(val)
      class(Bins3D_t), intent(in) :: this
      integer :: val

      val = this%num
   end function Bins3D_GetNum


   function Bins3D_GetTotNum(this) result(val)
      class(Bins3D_t), intent(in) :: this
      integer :: val

      val = this%totNum
   end function Bins3D_GetTotNum


   function Bins3D_GetBinSize(this) result(val)
      class(Bins3D_t), intent(in) :: this
      real(r8) :: val

      val = this%binSize
   end function Bins3D_GetBinSize


   function Bins3D_GetOffset(this) result(val)
      class(Bins3D_t), intent(in) :: this
      real(r8) :: val(3)

      val = this%offset
   end function Bins3D_GetOffset


   function Bins3D_GetShape(this) result(val)
      class(Bins3D_t), intent(in) :: this
      integer :: val(3)

      val = SHAPE(this%data)
   end function Bins3D_GetShape


   function Bins3D_InBounds3I(this, i1, i2, i3) result(val)
      class(Bins3D_t), intent(in) :: this
      integer, intent(in) :: i1, i2, i3
      logical :: val

      val = i1 >= 1 .and. i1 <= SIZE(this%data, 1) &
         .and. i2 >= 1 .and. i2 <= SIZE(this%data, 2) &
         .and. i3 >= 1 .and. i3 <= SIZE(this%data, 3)
   end function Bins3D_InBounds3I


   function Bins3D_InBounds3R(this, c1, c2, c3) result(val)
      class(Bins3D_t), intent(in) :: this
      real(r8), intent(in) :: c1, c2, c3
      logical :: val

      integer :: i1, i2, i3

      call this%CoordsToCell(c1, c2, c3, i1, i2, i3)
      val = this%InBounds(i1, i2, i3)
   end function Bins3D_InBounds3R


   subroutine Bins3D_CoordsToCell(this, c1, c2, c3, i1, i2, i3)
      class(Bins3D_t), intent(in) :: this
      real(r8), intent(in) :: c1, c2, c3
      integer, intent(out) :: i1, i2, i3

      i1 = FLOOR((c1 - this%offset(1)) / this%binSize) + 1
      i2 = FLOOR((c2 - this%offset(2)) / this%binSize) + 1
      i3 = FLOOR((c3 - this%offset(3)) / this%binSize) + 1
   end subroutine Bins3D_CoordsToCell


   subroutine Bins3D_CellToCoords(this, i1, i2, i3, c1, c2, c3)
      class(Bins3D_t), intent(in) :: this
      integer, intent(in) :: i1, i2, i3
      real(r8), intent(out) :: c1, c2, c3

      c1 = (i1 - .5_r8) * this%binSize + this%offset(1)
      c2 = (i2 - .5_r8) * this%binSize + this%offset(2)
      c3 = (i3 - .5_r8) * this%binSize + this%offset(3)
   end subroutine Bins3D_CellToCoords


   function Bins3D_GetValue3I(this, i1, i2, i3) result(val)
      class(Bins3D_t), intent(in) :: this
      integer, intent(in) :: i1, i2, i3
      integer :: val

      val = 0
      if (this%InBounds(i1, i2, i3)) val = this%data(i1, i2, i3)
   end function Bins3D_GetValue3I


   function Bins3D_GetValue3R(this, c1, c2, c3) result(val)
      class(Bins3D_t), intent(in) :: this
      real(r8), intent(in) :: c1, c2, c3
      integer :: val

      integer :: i1, i2, i3

      call this%CoordsToCell(c1, c2, c3, i1, i2, i3)
      val = 0
      if (this%InBounds(i1, i2, i3)) val = this%data(i1, i2, i3)
   end function Bins3D_GetValue3R


   subroutine Bins3D_IncValue(this, c1, c2, c3, iflag)
      class(Bins3D_t), intent(inout) :: this
      real(r8), intent(in) :: c1, c2, c3
      integer, intent(out) :: iflag

      integer :: i1, i2, i3

      this%totNum = this%totNum + 1
      call this%CoordsToCell(c1, c2, c3, i1, i2, i3)
      if (this%InBounds(i1, i2, i3)) then
         this%data(i1, i2, i3) = this%data(i1, i2, i3) + 1
         this%num = this%num + 1
         iflag = 0
      else
         iflag = 1
      end if
   end subroutine Bins3D_IncValue


   subroutine Bins3D_FindMax(this, coords, radius, coordsMax, iflag)
      class(Bins3D_t), intent(in) :: this
      real(r8), intent(in) :: coords(3)
      integer, intent(in) :: radius
      real(r8), intent(out) :: coordsMax(3)
      integer, intent(out) :: iflag

      integer :: i1, i2, i3, j1, j2, j3, k1, k2, k3, val, newVal
      logical :: maxNotFound

      call this%CoordsToCell(coords(1), coords(2), coords(3), i1, i2, i3)
      if (.not. this%InBounds(i1, i2, i3)) then
         coordsMax = coords
         iflag = 1
         return
      end if
      val = this%data(i1, i2, i3)
      maxNotFound = .true.
      do while(maxNotFound)
         maxNotFound = .false.
         k1 = i1
         k2 = i2
         k3 = i3
         do j3 = i3 - radius, i3 + radius
            do j2 = i2 - radius, i2 + radius
               do j1 = i1 - radius, i1 + radius
                  if (.not. this%InBounds(j1, j2, j3)) cycle
                  newVal = this%data(j1, j2, j3)
                  if (newVal > val) then
                     maxNotFound = .true.
                     k1 = j1
                     k2 = j2
                     k3 = j3
                     val = newVal
                  end if
               end do
            end do
         end do
         i1 = k1
         i2 = k2
         i3 = k3
      end do
      call this%CellToCoords(i1, i2, i3, coordsMax(1), coordsMax(2), coordsMax(3))
      iflag = 0
   end subroutine Bins3D_FindMax


   subroutine Bins3D_AllNodes(this)
      class(Bins3D_t), intent(inout) :: this
      integer, dimension(SIZE(this%data)) :: sendData, recvData
      integer :: sendNum, recvNum, sendTotNum, recvTotNum

      sendData = RESHAPE(this%data, SHAPE(sendData))
      sendNum = this%num
      sendTotNum = this%totNum
      recvData = 0
      recvNum = 0
      recvTotNum = 0
#ifdef MPI
      call MPI_ALLREDUCE(sendData, recvData, SIZE(this%data), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD)
      call MPI_ALLREDUCE(sendNum, recvNum, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD)
      call MPI_ALLREDUCE(sendTotNum, recvTotNum, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD)
#else
      recvData = sendData
      recvNum = sendNum
      recvTotNum = sendTotNum
#endif
      this%data = RESHAPE(recvData, SHAPE(this%data))
      this%num = recvNum
      this%totNum = recvTotNum
   end subroutine Bins3D_AllNodes
end module bins3D_m
