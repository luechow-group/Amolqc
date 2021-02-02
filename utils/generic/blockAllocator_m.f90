! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module blockAllocator_m
    use kinds_m, only : r8, i4
    implicit none
    private
    public :: BlockAllocator_t

    type :: BlockAllocator_t
        integer(i4) :: size = 0
        integer(i4) :: block_size = 100
        real(r8), allocatable :: vector(:)
    contains
        procedure, pass :: Add => BlockAllocator_t_add
        procedure, pass :: Get => BlockAllocator_t_get
        procedure, pass :: Create => BlockAllocator_t_create
        procedure, pass :: Destroy => BlockAllocator_t_destroy
    end type BlockAllocator_t

contains
    pure subroutine BlockAllocator_t_create(this,block_size)
        class(BlockAllocator_t), intent(inout) :: this
        integer(i4), intent(in), optional :: block_size

        if (present(block_size)) then
            this%block_size = block_size
        end if

        allocate(this%vector(block_size))
        this%vector = 0._r8
    end subroutine BlockAllocator_t_create

    pure subroutine BlockAllocator_t_destroy(this)
        class(BlockAllocator_t), intent(inout) :: this

        deallocate(this%vector)
        this%size = 0
    end subroutine BlockAllocator_t_destroy

    pure function BlockAllocator_t_get(this) result(vector)
        class(BlockAllocator_t), intent(in) :: this
        real(r8) :: vector(this%size)

        vector = this%vector(1:this%size)
    end function BlockAllocator_t_get

    pure subroutine BlockAllocator_t_add(this, r)
        class(BlockAllocator_t), intent(inout) :: this
        real(r8), intent(in) :: r(:)
        real(r8), allocatable :: new_vector(:)

        do while ( size(r) + this%size > size(this%vector) )
            ! enlarges this%vector by one block_size if r does not fit
            allocate( new_vector( size(this%vector) + this%block_size ) )
            new_vector = 0._r8
            new_vector(1:this%size) = this%vector(1:this%size)
            ! MOVE_ALLOC automatically deallocates new_vector
            call MOVE_ALLOC(new_vector, this%vector)
        end do

        this%vector(this%size + 1:this%size + size(r)) = r
        this%size = this%size + size(r)
    end subroutine BlockAllocator_t_add

end module blockAllocator_m