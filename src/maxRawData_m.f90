! Copyright (C) 2018 Michael A. Heuer
!
! SPDX-License-Identifier: GPL-3.0-or-later

module maxRawData_m
    use kinds_m, only : r8, i4
#ifdef MPI
    use MPI_F08
#endif
    use global_m, only : getNElec, iul, baseName
    use wfData_m, only : getNAlpha, getNBeta, atoms
    use globalUtils_m, only : mytid, nproc, MASTER
    use parsing_m, only: getinta, getstra
    use mpiInterface_m, only : myMPIBcastInteger, myMPI_Recv_DynamicReal8
    use blockAllocator_m, only: BlockAllocator_t

    implicit none
    private
    public :: maxraw_init, maxraw_destroy, maxraw_add, maxraw_writeResults, maxraw_writeParams, &
            maxraw_isInitialized

    integer(i4) :: mVerbose = 0, mBinaryFileUnit = 800
    logical :: mIsInitialized = .false.
    character(len = 120) :: mFileBasename
    type(BlockAllocator_t) :: serializedDataContainer
    integer(i4) :: serializedDataLength
    integer(i4) :: numberOfRecords = 0
    integer(i4) :: maxNumberOfRecords
    integer(i4) :: numberOfFiles = 0

contains

    subroutine maxraw_init(lines, nl)
        character(len = 120), intent(in) :: lines(:)
        integer(i4), intent(in) :: nl
        integer(i4) :: initSuccess = 0, fileUnit, iflag
        character(len = 120) :: filename
        
        call getinta(lines, nl, 'verbose=', mVerbose, iflag)
        call getstra(lines, nl, 'basename=', mFileBasename, iflag)
        if (iflag == 1) mFileBaseName = baseName
        call getinta(lines,nl,'max_records=',maxNumberOfRecords,iflag)
        if (iflag == 1) maxNumberOfRecords = 1000000

        serializedDataLength = 7*getNElec() + 1
        call serializedDataContainer%create(serializedDataLength*100)

        if (MASTER) then
            write(filename, '(a,a,I2.2,a)') trim(mFileBasename),'-',numberOfFiles,'.bin'
            open(mBinaryFileUnit, file = filename, form = 'unformatted')
            numberOfFiles = numberOfFiles + 1

            inquire(file = filename, number = fileUnit)

            if(fileUnit == mBinaryFileUnit) initSuccess = 1

            write(mBinaryFileUnit) size(atoms)
            numberOfRecords = numberOfRecords +1

            call writeAtoms()

            write(mBinaryFileUnit) getNElec()
            write(mBinaryFileUnit) getNAlpha()
            numberOfRecords = numberOfRecords +2
        end if

        call myMPIBcastInteger(initSuccess, 1)
        mIsInitialized = (initSuccess == 1)

        if (MASTER .and. mVerbose >= 0) call maxraw_writeParams(iul)
    end subroutine maxraw_init


    subroutine openNewFile()
        character(len = 120):: filename

        if (MASTER) then
            !close old one
            close(mBinaryFileUnit)

            write(filename, '(a,a,I2.2,a)') trim(mFileBasename),'-',numberOfFiles,'.bin'
            open(mBinaryFileUnit, file = filename, form = 'unformatted')
            numberOfFiles = numberOfFiles + 1

            inquire(file = filename, number = mBinaryFileUnit)
            numberOfRecords = 0
        end if
    end subroutine

    subroutine writeAtoms()
        real(r8) :: atomsVec(size(atoms)*3)
        integer :: atomicNumbersVec(size(atoms))
        integer :: i

        do i = 1, size(atoms)
            atomicNumbersVec(i) = atoms(i)%elemIdx
            atomsVec(i*3-2) = atoms(i)%cx
            atomsVec(i*3-1) = atoms(i)%cy
            atomsVec(i*3-0) = atoms(i)%cz
        end do

        write(mBinaryFileUnit) atomicNumbersVec
        write(mBinaryFileUnit) atomsVec
        numberOfRecords = numberOfRecords+2
    end subroutine

    subroutine maxraw_destroy()
        close(mBinaryFileUnit)
    end subroutine maxraw_destroy

    function maxraw_isInitialized() result(res)
        logical res
        res = mIsInitialized
    end function maxraw_isInitialized
    
    subroutine maxraw_writeParams(iu)
        integer, intent(in) :: iu
        write(iu, *) 'Raw data is written in binary format into "', trim(mFileBasename), '"'
        write(iu, *) 'in the following order, dimensionality, and precision:'
        write(iu, *) '---------------Header----------------'
        write(iu, *) '              number of atoms ( 1*i4)'
        write(iu, *) '        atomic numbers vector ( N*i4)'
        write(iu, *) '        atom positions vector (3N*r8)'
        write(iu, *) '          number of electrons ( 1*i4)'
        write(iu, *) '    number of alpha electrons ( 1*i4)'
        write(iu, *) '----------------Body-----------------'
        write(iu, *) '  1.  sample positions vector (3N*r8)'
        write(iu, *) '      kinetic energies vector ( N*r8)'
        write(iu, *) '     maximum positions vector (3N*r8)'
        write(iu, *) '             -ln|Psi|^2 value ( 1*r8)'
        write(iu, *) '  2. ...                             '
    end subroutine maxraw_writeParams

    subroutine maxraw_add(sample, kinetic_energies, maximum, value)
        real(r8), intent(in) :: sample(:)
        real(r8), intent(in) :: kinetic_energies(:)
        real(r8), intent(in) :: maximum(:)
        real(r8), intent(in) :: value
        real(r8), allocatable :: serializedData(:)
        integer(i4) :: ne

        ne = getNElec()
        allocate(serializedData(serializedDataLength))

        serializedData(1 : 3*ne) = sample
        serializedData(3*ne+1 : 4*ne) = kinetic_energies
        serializedData(4*ne+1 : 7*ne) = maximum
        serializedData(7*ne+1) = value

        call serializedDataContainer%add(serializedData)

    end subroutine maxraw_add

    subroutine writeSerializedDataToBinary(allSerializedData)
        real(r8), intent(in) :: allSerializedData(:)
        integer :: ne, l, i

        ne = getNElec()
        l = serializedDataLength

        do i = 0, size(allSerializedData)/l -1
            if(numberOfRecords + 1 > maxNumberOfRecords) call openNewFile()

            write(mBinaryFileUnit) allSerializedData(i*l+1 : i*l+l)
            numberOfRecords = numberOfRecords+1

        end do
    end subroutine writeSerializedDataToBinary

    subroutine maxraw_writeResults()
        real(r8), allocatable :: allSerializedData(:)
        
#ifdef MPI
        integer(i4) :: p
#endif

        if(MASTER) then
            ! write own block
            allSerializedData =  serializedDataContainer%get()
            call serializedDataContainer%destroy()
            call writeSerializedDataToBinary(allSerializedData)
        end if
        
#ifdef MPI
        if (.not.MASTER) then 
            ! send block
            allSerializedData = serializedDataContainer%get()
            call serializedDataContainer%destroy()
            call MPI_SEND(allSerializedData, size(allSerializedData), MPI_REAL8, 0, 0, MPI_COMM_WORLD)
        else
            ! recieve other blocks
            do p = 1,nproc-1
                call myMPI_Recv_DynamicReal8(allSerializedData,MPI_ANY_SOURCE,0)
                call writeSerializedDataToBinary(allSerializedData)
            end do
        end if
#endif

        if (MASTER) then
            write(iul, *) 'Raw data was written in file ', trim(mFileBasename), '.'
        end if
    end subroutine maxraw_writeResults

end module maxRawData_m
