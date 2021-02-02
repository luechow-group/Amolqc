! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

! this module is mostly an implementation of pseudo-code from https://doi.org/10.1002/jcc.20746
module statistic_m
    use kinds_m, only: i8, r8
    use, intrinsic :: iso_c_binding, only: c_int64_t, c_double
    
    implicit none
    private
    public :: Statistic_t, StatisticData_t, &
              operator (+), operator (==), operator (.similar.), assignment (=)

    type, bind(C) :: StatisticData_t
        integer(c_int64_t) :: nSamp = 0_c_int64_t
        real(c_double) :: sum   = 0._c_double
        real(c_double) :: sqSum = 0._c_double
    end type

    type :: Statistic_t
        integer(i8), pointer :: nSamp => null()
        real(r8), pointer :: sum => null()
        real(r8), pointer :: sqSum => null()
    contains
        procedure, pass :: Add_data           => Statistic_t_add_data
        procedure, pass :: Write              => Statistic_t_write
        procedure, pass :: Get_average        => Statistic_t_get_average
        procedure, pass :: Get_error_estimate => Statistic_t_get_error_estimate
        procedure, pass :: Get_errors_error   => Statistic_t_get_errors_error
        procedure, pass :: Create             => Statistic_t_create
        procedure, pass :: Destroy            => Statistic_t_destroy
        procedure, pass :: Nullify            => Statistic_t_nullify
    end type Statistic_t

    interface operator (+)
        module procedure StatisticData_t_add
    end interface operator (+)

    interface operator (==)
        module procedure Statistic_t_is_equal
    end interface operator (==)

    interface operator (.similar.)
        module procedure Statistic_t_is_similar
    end interface operator (.similar.)

    interface assignment (=)
        module procedure Statistic_t_assign
    end interface assignment (=)

contains
    ! interfaces
    function Statistic_t_is_equal(this, other) result(bool)
        type(Statistic_t), intent(in) :: this, other
        logical :: bool

        bool = .true.
        if (.not. this%nSamp == other%nSamp) then
            bool = .false.
            return
        else if (.not. this%sum == other%sum) then
            bool = .false.
            return
        else if (.not. this%sqSum == other%sqSum) then
            bool = .false.
            return
        end if
    end function Statistic_t_is_equal

    function Statistic_t_is_similar(this, other) result(bool)
        type(Statistic_t), intent(in) :: this, other
        logical                     :: bool

        bool = .true.
        if (.not. this%nSamp == other%nSamp) then
            bool = .false.
            return
        else if (.not. this%sum == other%sum) then
            bool = .false.
            return
        end if
    end function Statistic_t_is_similar

    function StatisticData_t_add(this, other) result(sum)
        type(StatisticData_t), intent(in) :: this, other
        type(StatisticData_t)             :: sum

        sum%nSamp = this%nSamp + other%nSamp
        sum%sum   = this%sum   + other%sum
        sum%sqSum = this%sqSum + other%sqSum
    end function StatisticData_t_add

    subroutine Statistic_t_assign(this, other)
        type(Statistic_t), intent(inout) :: this
        type(Statistic_t), intent(in)    :: other

        this%nSamp = other%nSamp
        this%sum   = other%sum
        this%sqSum = other%sqSum
    end subroutine Statistic_t_assign

    ! create and destroy
    subroutine Statistic_t_create(this, data)
        class(Statistic_t),            intent(inout) :: this
        type(StatisticData_t), target, intent(in)    :: data
        this%nSamp => data%nSamp
        this%sum   => data%sum
        this%sqSum => data%sqSum
    end subroutine Statistic_t_create

    subroutine Statistic_t_destroy(this)
        class(Statistic_t), intent(inout) :: this
        this%nSamp = 0
        this%sum = 0._r8
        this%sqSum = 0._r8

        call this%Nullify()
    end subroutine Statistic_t_destroy

    subroutine Statistic_t_nullify(this)
        class(Statistic_t), intent(inout) :: this

        nullify(this%nSamp)
        nullify(this%sum)
        nullify(this%sqSum)
    end subroutine Statistic_t_nullify

    ! other type procedures
    subroutine Statistic_t_add_data(this, newSample)
        class(Statistic_t), intent(inout) :: this
        real(r8),     intent(in)    :: newSample

        this%nSamp = this%nSamp + 1
        this%sum   = this%sum   + newSample
        this%sqSum = this%sqSum + newSample**2
    end subroutine Statistic_t_add_data

    pure function Statistic_t_get_average(this) result(average)
        class(Statistic_t), intent(in) :: this
        real(r8)                 :: average

        average = this%sum / this%nSamp
    end function Statistic_t_get_average

    pure function Statistic_t_get_error_estimate(this) result(error)
        class(Statistic_t), intent(in) :: this
        real(r8)                 :: error

        error = sqrt( ( this%sqSum / this%nSamp - this%sum**2 / this%nSamp**2 ) / (this%nSamp - 1) )
    end function Statistic_t_get_error_estimate

    pure function Statistic_t_get_errors_error(this) result(errorsError)
        class(Statistic_t), intent(in) :: this
        real(r8)                 :: errorsError

        errorsError = this%Get_error_estimate() * sqrt(2._r8 / ( this%nSamp - 1 ) )
    end function Statistic_t_get_errors_error

    subroutine Statistic_t_write(this, iu)
        class(Statistic_t), intent(in) :: this
        integer,          intent(in) :: iu
        write(iu,*) this%nSamp, this%sum, this%sqSum
    end subroutine Statistic_t_write

end module statistic_m
