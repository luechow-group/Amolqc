! Copyright (C) 2012 Alexander Sturm
!
! SPDX-License-Identifier: GPL-3.0-or-later

program randomTest_p
    use kinds_m, only: r8, i8
    use mrg_m, only: init_mrgran, mrg_ran, mrg_gran, mrg_intran
    use mt_m, only: mt_ran, mt_gran, init_mtran
    implicit none

    character(len=5) :: rng
    character(len=10) :: tmp
    integer :: i
    real(r8) :: cpuStart, cpuEnd, wallStart, wallEnd, s
    integer :: seed = 1633837925
    integer,parameter :: ops = 10000000
    rng = "mrg"

    if(COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1, rng)
        if(COMMAND_ARGUMENT_COUNT() > 1) then
            call GET_COMMAND_ARGUMENT(2, tmp)
            read(tmp, "(I10)") seed
            write(0, *) "Using seed", seed
        endif
    endif

    if(rng == "mrg") then
        write(0, *) "Testing MRG8 (integer)"
        call test_mrg_int()
    elseif(rng == "mrgr") then
        write(0, *) "Testing MRG8 (real)"
        call test_mrg_real()
    elseif(rng == "mt") then
        write(0, *) "Testing MT19937"
        call test_mt()
    elseif(rng == "test") then
        write(0, *) "MRG8 rands"
        write(*,*) init_mrgran(seed)
        !write(*,*) mrg_intran(), mrg_intran(), mrg_intran(), mrg_intran(), mrg_intran()
        write(*,*) mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), &
                   mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran()
    elseif(rng == "speed") then
        write(*,*) "Testing MRG8 speed"
        s = init_mrgran(seed)
        call startTime()
        do i = 1, ops
            s = s + mrg_ran()
        enddo
        call endTime()
        write(*,*) "Result:", s


        write(*,*) "Testing MT19937 speed"
        s = init_mtran(seed)
        call startTime()
        do i = 1, ops
            s = s + mt_ran()
        enddo
        call endTime()
        write(*,*) "Result:", s
    endif


contains
    subroutine print_int(foo, bytes)
        integer(i8), intent(in) :: foo
        integer, intent(in) :: bytes
        integer(i8) :: i
        do i = 0, bytes - 1
#ifndef NVHPC
            write(*,"(A)",advance="no") char(iand(shiftr(foo, i*8), 255_i8))
#else
            write(*,"(A)",advance="no") char(iand(rshift(foo, i*8), 255_i8))
#endif
        enddo

    end subroutine

    subroutine print_reals(m, n)
        real(r8), intent(in) :: m, n
        integer(i8) :: a, b, rand
        integer(i8), parameter :: k = 4503599627370495_i8 ! 2**52 - 1
        a = transfer(m, a)
        a = iand(a, k)
        b = transfer(n, b)
        b = iand(b, k)
        call print_int(a, 6)
        call print_int(b, 6)
#ifndef NVHPC
        rand = shiftl(shiftr(a, 48), 4) + shiftr(b, 48)
#else
        rand = lshift(rshift(a, 48), 4) + rshift(b, 48)
#endif
        call print_int(rand, 1)
    end subroutine
    
    subroutine startTime()
        integer(i8) :: cnt, rate
        call cpu_time(cpuStart)

        call system_clock(cnt,rate)
        wallStart = dble(cnt)/dble(rate)
    end subroutine

    subroutine endTime()
        integer(i8) :: cnt, rate
        call cpu_time(cpuEnd)

        call system_clock(cnt,rate)
        wallEnd = dble(cnt)/dble(rate)
        write(*,*) "CPU time:", (cpuEnd - cpuStart)
        write(*,*) "Walltime:", (wallEnd - wallStart)
        write(*,*) "Ops/sec:", ops/(wallEnd - wallStart)

    end subroutine

    subroutine test_mrg_int()
        real(r8) :: dummy
        integer(i8) :: a, b, c, d, rand
        integer, parameter :: m = 1073741823 ! 2**30 - 1
        integer :: i
        dummy = init_mrgran(seed)
        do i=1,1
            a = iand(mrg_intran(), m)
            b = iand(mrg_intran(), m)
#ifndef NVHPC
            rand = shiftl(a, 30) + b
#else
            rand = lshift(a, 30) + b
#endif
            call print_int(rand, 7)
            c = iand(mrg_intran(), m)
            d = iand(mrg_intran(), m)
#ifndef NVHPC
            rand = shiftl(c, 30) + d
#else
            rand = lshift(c, 30) + d
#endif
            call print_int(rand, 7)
#ifndef NVHPC
            rand = shiftl(shiftr(c, 26), 4) + shiftr(a, 26)
#else
            rand = lshift(rshift(c, 26), 4) + rshift(a, 26)
#endif
            call print_int(rand, 1)
        enddo
    end subroutine

    subroutine test_mrg_real()
        real(r8) :: dummy
        dummy = init_mrgran(seed)
        do
            call print_reals(mrg_ran(), mrg_ran())
        enddo
    end subroutine

    subroutine test_mt()
        real(r8) :: dummy
        dummy = init_mtran(seed)
        do
            call print_reals(mt_ran(), mt_ran())
        enddo
    end subroutine

end program

