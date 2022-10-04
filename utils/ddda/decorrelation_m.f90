! Copyright (C) 2018-2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

! this module is mostly an implementation of pseudo-code from https://doi.org/10.1002/jcc.20746
module decorrelation_m
    use, intrinsic :: iso_c_binding, only: c_bool, c_int64_t, c_double
    use error_m, only: assert
    use kinds_m, only: i8, r8, i4
    use boltzmann_m, only: Boltzmann_fit
    use statistic_m
    use plateau_m, only: Find_plateau

    implicit none

    public :: DecorrelationData_t, Decorrelation_t, &
            operator (+), operator (==), operator (.similar.)

    type, bind(C) :: DecorrelationData_t
        integer(c_int64_t)                 :: size     = 0_c_int64_t
        integer(c_int64_t)                 :: nSamp    = 0_c_int64_t
        real(c_double),      dimension(64) :: wSamps   = 0._c_double
        logical(c_bool),     dimension(64) :: wSampsEx = .false.
        type(StatisticData_t), dimension(64) :: stats
    end type DecorrelationData_t

    type :: Decorrelation_t
        integer(i8),                pointer :: size => null()
        integer(i8),                pointer :: nSamp => null()
        real(r8),    dimension(:), pointer :: wSamps => null()
        logical(kind=1), dimension(:), pointer :: wSampsEx => null()
        type(Statistic_t), dimension(64)         :: stats
    contains
        procedure, pass :: Add_data     => Decorrelation_t_add_data
        procedure, pass :: Write        => Decorrelation_t_write
        procedure, pass :: Write_results => Decorrelation_t_write_results
        procedure, pass :: Get_error     => Decorrelation_t_get_error
        procedure, pass :: Create       => Decorrelation_t_create
        procedure, pass :: Destroy      => Decorrelation_t_destroy
        procedure, pass :: Nullify      => Decorrelation_t_nullify
    end type Decorrelation_t

    interface operator (+)
        module procedure DecorrelationData_t_add
    end interface operator (+)

    interface operator (.similar.)
        module procedure Decorrelation_t_is_similar
    end interface operator (.similar.)

    interface operator (==)
        module procedure Decorrelation_t_is_equal
    end interface operator (==)

    contains
        !interfaces
        function DecorrelationData_t_add(AD, BD) result(CD)
            type(DecorrelationData_t), target, intent(in) :: AD, BD
            type(DecorrelationData_t), target     :: CD

            type(Decorrelation_t)                 :: A, B, C
            type(Statistic_t)                     :: statA, statB, newStat
            type(StatisticData_t), target         :: statAD, statBD, newStatD
            real(r8)                            :: carry, wSampleA, wSampleB
            logical                             :: carryEx, wSampleAEx, wSampleBEx
            integer                             :: i

            call A%Create(AD)
            call B%Create(BD)
            call C%Create(CD)
            call statA%Create(statAD)
            call statB%Create(statBD)
            call newStat%Create(newStatD)
            
            C%nSamp = A%nSamp + B%nSamp

            carryEx = .false.
            carry = 0d0

            do while(C%nSamp >= 2**C%size)
                C%size = C%size + 1
            end do

            do i=1,C%size

                if (i <= A%size) then
                    statA = A%stats(i)
                    wSampleA = A%wSamps(i)
                    wSampleAEx = A%wSampsEx(i)
                else
                    statA = newStat
                    wSampleA = 0d0
                    wSampleAEx = .false.
                end if

                if (i <= B%size) then
                    statB = B%stats(i)
                    wSampleB = B%wSamps(i)
                    wSampleBEx = B%wSampsEx(i)
                else
                    statB = newStat
                    wSampleB = 0d0
                    wSampleBEx = .false.
                end if

                CD%stats(i) = CD%stats(i) + statAD + statBD

                if (carryEx .and. wSampleAEx .and. wSampleBEx) then
                    call C%stats(i+1)%Add_data((wSampleA + wSampleB)/2)
                    C%wSamps(i) = carry
                    C%wSampsEx(i) = .true.
                    carryEx = .true.
                    carry = (wSampleA + wSampleB)/2
                else if ((.not. carryEx) .and. wSampleAEx .and. wSampleBEx) then
                    call C%stats(i+1)%Add_data((wSampleA + wSampleB)/2)
                    C%wSamps(i) = 0d0
                    C%wSampsEx(i) = .false.
                    carryEx = .true.
                    carry = (wSampleA + wSampleB)/2
                else if (carryEx .and. (.not. wSampleAEx) .and. wSampleBEx) then
                    call C%stats(i+1)%Add_data((carry + wSampleB)/2)
                    C%wSamps(i) = 0d0
                    C%wSampsEx(i) = .false.
                    carryEx = .true.
                    carry = (wSampleB + carry)/2
                else if (carryEx .and. wSampleAEx .and. (.not. wSampleBEx)) then
                    call C%stats(i+1)%Add_data((wSampleA + carry)/2)
                    C%wSamps(i) = 0d0
                    C%wSampsEx(i) = .false.
                    carryEx = .true.
                    carry = (wSampleA + carry)/2
                else if (carryEx .or. wSampleAEx .or. wSampleBEx) then
                    C%wSamps(i) = carry + wSampleA + wSampleB
                    C%wSampsEx(i) = .true.
                    carryEx = .false.
                    carry = 0d0
                else
                    C%wSamps(i) = 0d0
                    C%wSampsEx(i) = .false.
                    carryEx = .false.
                    carry = 0d0
                end if
            end do

            call A%Nullify()
            call B%Nullify()
            call C%Nullify()
            call statA%Nullify()
            call statB%Nullify()
            call newStat%Nullify()
        end function DecorrelationData_t_add

        function Decorrelation_t_is_similar(this, other) result(bool)
            type(Decorrelation_t), intent(in) :: this, other
            logical                         :: bool
            integer                         :: i

            bool = .true.
            if (.not. this%size == other%size) then
                bool = .false.
                return
            else if (.not. this%nSamp == other%nSamp) then
                bool = .false.
                return
            end if

            do i = 1, this%size
                if (.not. this%wSamps(i) == other%wSamps(i)) then
                    bool = .false.
                    return
                else if (.not. this%wSampsEx(i) .eqv. other%wSampsEx(i)) then
                    bool = .false.
                    return
                else if ((i == 1) .and. (.not. this%stats(i) == other%stats(i))) then
                    bool = .false.
                    return
                ! .similar. denotes: nearly equal (equal sum and nSamp)
                else if ((i > 1) .and. (.not. (this%stats(i) .similar. other%stats(i)))) then
                    bool = .false.
                    return
                end if
            end do

        end function Decorrelation_t_is_similar

        function Decorrelation_t_is_equal(this, other) result(bool)
            type(Decorrelation_t), intent(in) :: this, other
            logical                         :: bool
            integer                         :: i

            bool = .true.
            if (.not. this%size == other%size) then
                bool = .false.
                return
            else if (.not. this%nSamp == other%nSamp) then
                bool = .false.
                return
            end if

            do i = 1, this%size
                if (.not. this%wSamps(i) == other%wSamps(i)) then
                    bool = .false.
                    return
                else if (.not. this%wSampsEx(i) .eqv. other%wSampsEx(i)) then
                    bool = .false.
                    return
                else if ((i == 1) .and. (.not. this%stats(i) == other%stats(i))) then
                    bool = .false.
                    return
                else if ((i > 1) .and. (.not. this%stats(i) == other%stats(i))) then
                    bool = .false.
                    return
                end if
            end do

        end function Decorrelation_t_is_equal

        subroutine Decorrelation_t_write(this, iu)
            class(Decorrelation_t), intent(in) :: this
            integer,              intent(in) :: iu
            integer                          :: i

            write(iu,*) this%size, this%nSamp
            do i = 1, this%size
                write(iu,*) i, this%wSampsEx(i), this%wSamps(i)
                call this%stats(i)%Write(iu)
            end do
        end subroutine Decorrelation_t_write

        subroutine Decorrelation_t_write_results(this, iu)
            class(Decorrelation_t), intent(in) :: this
            integer, intent(in) :: iu
            real(r8) :: error, errorsError
            integer(i8) :: blockLength, k
            logical :: converged

            call this%get_error(error, errorsError, blockLength, converged)

            write(iu,'(a,f14.5,a)')  ' variance                     =', &
                    this%stats(1)%get_error_estimate()**2 * this%stats(1)%nSamp, ' E_h^2'

            if (.not. converged) then
                write(iu,'(a)')  ' WARNING: Plateau has not been reached, following values are not reliable!'
            end if

            write(iu,'(a,f14.5,a)')  ' error                        =', &
                    error, ' E_h'
            ! write(iu,'(a,f9.7,a)')  ' errors error                 =       ', errorsError, ' E_h'
            write(iu,'(a,i8)')       ' block length                 =', blockLength
            write(iu,'(a,f11.2)')    ' n_corr                       =',&
                    (error / this%stats(1)%get_error_estimate())**2
            write(iu,'(a)')  ''
            write(iu,'(a)')  ' log2(blen)   error estimate   errors error'
            do k = 1,this%size - 1
                write(iu,'(a,i10,f17.7,f15.7)') ' ',k-1, this%stats(k)%get_error_estimate(), this%stats(k)%get_errors_error()
            end do

        end subroutine Decorrelation_t_write_results

        subroutine Decorrelation_t_create(this, data)
            class(Decorrelation_t),            intent(inout) :: this
            type(DecorrelationData_t), target, intent(in)    :: data
            integer                                        :: i

            this%size     => data%size
            this%nSamp    => data%nSamp
            this%wSamps   => data%wSamps
            this%wSampsEx => data%wSampsEx

            do i = 1, 64
                call this%stats(i)%Create(data%stats(i))
            end do

        end subroutine Decorrelation_t_create

        subroutine Decorrelation_t_destroy(this)
            class(Decorrelation_t), intent(inout) :: this
            integer                             :: i

            this%size = 0
            this%nSamp = 0
            this%wSamps = 0._r8
            this%wSampsEx = .false.

            nullify(this%size)
            nullify(this%nSamp)
            nullify(this%wSamps)
            nullify(this%wSampsEx)

            do i = 1, 64
                call this%stats(i)%Destroy
            end do
        end subroutine Decorrelation_t_destroy

        subroutine Decorrelation_t_nullify(this)
            class(Decorrelation_t), intent(inout) :: this
            integer                             :: i

            nullify(this%size)
            nullify(this%nSamp)
            nullify(this%wSamps)
            nullify(this%wSampsEx)

            do i = 1, 64
                call this%stats(i)%Nullify
            end do
        end subroutine Decorrelation_t_nullify

        subroutine Decorrelation_t_add_data(this, newSample)
            class(Decorrelation_t), intent(inout) :: this
            real(r8),         intent(in)    :: newSample
            real(r8)                        :: carry
            integer                             :: i

            call assert(associated(this%nSamp), "Decorrelation_t_add_data: not associated") 
            this%nSamp = this%nSamp + 1

            ! Lengthen the vectors, when necessary, to accomodate
            if (this%nSamp >= 2**this%size) then
                this%size = this%size + 1
            end if

            call this%stats(1)%Add_data(newSample)

            carry = newSample
            ! Propagate the new sample up through the data structure
            do i=1,this%size
                if (this%wSampsEx(i)) then
                    carry = (this%wSamps(i) + carry) / 2
                    call this%stats(i+1)%Add_data(carry)
                    this%wSampsEx(i) = .false.
                    !next line is unnecessary, but makes things more clean
                    this%wSamps(i) = 0d0
                else
                    this%wSampsEx(i) = .true.
                    this%wSamps(i) = carry
                    exit
                end if
            end do

        end subroutine Decorrelation_t_add_data

        subroutine Decorrelation_t_get_error(this, error, errorsError, blockLength, converged)
            class(Decorrelation_t), intent(in)  :: this
            real(r8), intent(out) :: error, errorsError
            integer(i8), intent(out) :: blockLength
            logical, intent(out) :: converged

            logical :: plateau(this%size - 1)
            integer :: i
            real(r8) :: errorEstimates(this%size - 1), errorsErrors(this%size - 1)

            converged = .false.

            do i = 1, this%size - 1
                errorEstimates(i) = this%stats(i)%get_error_estimate()
                errorsErrors(i) = this%stats(i)%get_errors_error()
            end do

            plateau = Find_plateau(errorEstimates)
            if (ANY(plateau)) then
                converged = .true.
                error = SUM(errorEstimates, mask = plateau) / COUNT(plateau)
                errorsError = (MAXVAL(errorEstimates, mask = plateau) - MINVAL(errorEstimates, mask = plateau)) &
                        / 2._r8 + MINVAL(errorsErrors, mask = plateau)
                do i = this%size - 1, 1, -1
                    if (plateau(i)) blockLength = 2 ** (i - 1)
                end do
            else
                error = errorEstimates(MAX(this%size - 1,1_c_int64_t))
                errorsError = errorsErrors(MAX(this%size - 1,1_c_int64_t))
                blockLength = 2 ** (MAX(this%size - 2,0_c_int64_t))
            end if
        end subroutine Decorrelation_t_get_error

end module decorrelation_m
