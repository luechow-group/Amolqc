! Copyright (C) 2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module rhoData_m
   use kinds_m, only: r8
   use error_m, only: asserts, assert
   use wfData_m, only: atoms, getNNuc
   use atom_m, only: atom_getPosition
   use partition_m, only: partition_t, assignment (=)
#ifdef MPI
   use MPI_F08
   use globalUtils_m, only: mytid, nproc, MASTER
   use mpiInterface_m, only: myMPI_Recv_DynamicInteger4
#else
   use globalUtils_m, only: MASTER
#endif
   implicit none
   private
   public :: rhoData_t, operator (+), assignment (=)

   type rhoData_t
      !private
      integer :: n !number of partitions
      integer :: nSamples ! number of total samples
      type(partition_t), allocatable :: partitions(:)
      real(r8) :: assignThresh, preAssignThresh
      real(r8) :: printThresh
      integer, allocatable :: fragments(:)
   contains
      procedure :: init => rhoData_init
      procedure :: addData => rhoData_addData
      procedure :: addPartition => rhoData_addPartition
      procedure :: addInstance => rhoData_addInstance
      procedure :: writeResults => rhoData_writeResults
      procedure :: getBasin => rhoData_getBasin
   end type rhoData_t

   interface assignment (=)
      module procedure rhoData_copy
   end interface

   interface operator (+)
      module procedure rhoData_sum
   end interface

   interface rhoData_t
      module procedure rhoData_constructor
   end interface rhoData_t

contains
   subroutine rhoData_init(this, assignThresh, printThresh, fragments, preAssignThresh)
      class(rhoData_t), intent(inout) :: this
      real(r8), intent(in) :: assignThresh, printThresh
      integer, optional, intent(in) :: fragments(:)
      real(r8), optional, intent(in) :: preAssignThresh
      integer :: i

      if (allocated(this%partitions)) deallocate(this%partitions)
      if (allocated(this%fragments)) deallocate(this%fragments)

      ! next lines use allocation on assignment
      if (PRESENT(fragments)) then
         this%fragments = fragments
      else
         this%fragments = [(i, i=1, getNNuc())]
      end if

      allocate(this%partitions(0))
      this%n = 0
      this%nSamples = 0
      this%assignThresh = assignThresh
      this%printThresh = printThresh
      if (PRESENT(preAssignThresh)) this%preAssignThresh = preAssignThresh
   end subroutine rhoData_init

   function rhoData_constructor(assignThresh, printThresh, fragments, preAssignThresh) result(rhoData)
      real(r8), intent(in) :: assignThresh, printThresh
      type(rhoData_t) :: rhoData
      integer, optional, intent(in) :: fragments(:)
      real(r8), optional, intent(in) :: preAssignThresh

      call rhoData%init(assignThresh, printThresh, fragments=fragments, preAssignThresh=preAssignThresh)

   end function rhoData_constructor

   subroutine rhoData_copy(this, other)
      class(rhoData_t), intent(inout) :: this
      class(rhoData_t), intent(in) :: other

      call this%init(other%assignThresh, other%printThresh)
      this%n = other%n
      this%nSamples = other%nSamples

      deallocate(this%partitions)
      ! using allocate on assign
      this%partitions = other%partitions
   end subroutine rhoData_copy

   function rhoData_getBasin(this, r, preAssign) result(basin)
      class(rhoData_t), intent(in) :: this
      real(r8), intent(in) :: r(3)
      logical, optional, intent(in) :: preAssign
      integer :: basin
      integer :: j
      logical :: found
      real(r8) :: minDist, dist

      basin = 0
      found = .false.
      if (PRESENT(preAssign) .and. preAssign) then
         minDist = this%preAssignThresh
      else
         minDist = this%assignThresh
      end if
      do j=1, getNNuc()
         dist = NORM2(r - atom_getPosition(atoms(j)))
         if (dist < minDist) then
            if (asserts) call assert(.not. found, 'rhoData_getBasin: assignThresh too large.')
            found = .true.
            basin = this%fragments(j)
            minDist = dist
         end if
      end do
   end function rhoData_getBasin

   subroutine rhoData_addData(this, basins)
      class(rhoData_t), intent(inout) :: this
      integer, intent(in) :: basins(:)
      integer :: electronsPerFragment(MAXVAL(this%fragments) + 1)
      integer :: i

      electronsPerFragment = 0
      do i=1,SIZE(basins)
         if (basins(i) /= 0) then
            electronsPerFragment(basins(i)) = electronsPerFragment(basins(i)) + 1
         else
            electronsPerFragment(MAXVAL(this%fragments) + 1) = electronsPerFragment(MAXVAL(this%fragments) + 1) + 1
         end if
      end do

      call this%addPartition(partition_t(electronsPerFragment, 1))

   end subroutine rhoData_addData

   subroutine rhoData_addPartition(this, partition)
      class(rhoData_t), intent(inout) :: this
      type(partition_t), intent(in) :: partition
      type(partition_t) :: partitions(this%n + 1) ! needed iff new partition
      integer :: i
      logical :: sorted

      sorted = .false.
      do i=1, this%n
         if (ALL(this%partitions(i)%electronsPerFragment == partition%electronsPerFragment)) then
            this%partitions(i)%count = this%partitions(i)%count + partition%count
            sorted = .true.
            exit
         end if
      end do

      if (.not. sorted) then

         if (this%n /= 0) then
            call assert(SUM(this%partitions(1)%electronsPerFragment)==SUM(partition%electronsPerFragment), &
            'rhoData_addPartition: invalid size of elecronsPerCore')
         end if

         ! copying partitions
         do i=1, this%n
            partitions(i) = this%partitions(i)
         end do

         ! setting new partition
         partitions(this%n + 1) = partition
         this%partitions = partitions
         this%n = this%n + 1
      end if

      this%nSamples = this%nSamples + partition%count
   end subroutine rhoData_addPartition

   subroutine rhoData_addInstance(this, other)
      class(rhoData_t), intent(inout) :: this
      class(rhoData_t), intent(in) :: other
      integer :: i

      do i=1, other%n
         call this%addPartition(other%partitions(i))
      end do
   end subroutine rhoData_addInstance

   function rhoData_sum(this, other) result(sum)
      class(rhoData_t), intent(in) :: this, other
      type(rhoData_t) :: sum

      sum = this
      call sum%addInstance(other)
   end function rhoData_sum

   subroutine rhoData_writeResults(this, iul)
      class(rhoData_t), intent(inout) :: this
      integer, intent(in) :: iul
      real(r8) :: integral, sum
      integer :: i, j

#ifdef MPI
      integer :: p, occLen
      integer, allocatable :: all_counts(:)
      integer, allocatable :: all_electronsPerFragment(:)

      occLen = MAXVAL(this%fragments) + 1

      if (.not.MASTER) then
         allocate(all_counts(this%n))
         allocate(all_electronsPerFragment(this%n*occLen))
         do i = 1,this%n
            all_counts(i) = this%partitions(i)%count
            do j = 1, occLen
               all_electronsPerFragment((i-1)*occLen+j) = this%partitions(i)%electronsPerFragment(j)
            end do
         end do
         call MPI_SEND(all_counts, this%n, MPI_INTEGER4, 0, 0, MPI_COMM_WORLD)
         call MPI_SEND(all_electronsPerFragment, this%n*occLen, MPI_INTEGER4, 0, 1, MPI_COMM_WORLD)
      else
         do p = 1,nproc-1
            call myMPI_Recv_DynamicInteger4(all_counts,p,0)
            call myMPI_Recv_DynamicInteger4(all_electronsPerFragment,p,1)
            do i = 1, SIZE(all_counts)
               call this%addPartition(partition_t(all_electronsPerFragment((i - 1)*occLen + 1:i*occLen), all_counts(i)))
            end do
         end do
      end if
#endif

      if (MASTER) then
         write(iul,*) 'Converged minimizers:', this%nSamples
         write(iul,*) ''

         write(iul,'(A)') ' Fragments'
         write(iul,'(99i3)') this%fragments
         write(iul,*) ''

         write(iul,'(A)') ' Weight     Partition'
         do i = 1, this%n
            if (1._r8 * this%partitions(i)%count/this%nSamples > this%printThresh) then
               call this%partitions(i)%write(iul, this%nSamples)
            end if
         end do

         write(iul,*) ''
         write(iul,'(A)') ' Integrals of Basins'
         sum = 0._r8
         do i = 1, MAXVAL(this%fragments) + 1
            integral = 0._r8
            do j = 1, this%n
               integral = integral + this%partitions(j)%electronsPerFragment(i) * this%partitions(j)%count
            end do
            sum = sum + integral
            write(iul,'(i4,A,f9.4)') i, ':', integral / this%nSamples
         end do
         write(iul,'(A,f9.4)') ' sum:', sum / this%nSamples
      end if
   end subroutine rhoData_writeResults

end module rhoData_m
