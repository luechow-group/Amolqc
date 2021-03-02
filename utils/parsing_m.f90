! Copyright (C) 1998-1999, 2006-2008, 2015, 2018 Arne Luechow
! Copyright (C) 2003 Christian Diedrich
!
! SPDX-License-Identifier: GPL-3.0-or-later

!     -------------------------------------------------
module parsing_m
    use error_m, only: error
    use kinds_m, only : i8, r8
    use string_utility_module, only: str_parse_next_token, str_is_real, str_is_integer
    implicit none
contains

function Py_split(string) result(words)
    character(len=*), intent(in) :: string
    character(len=:), allocatable :: words(:)
    integer :: words_count

    character(len=:), allocatable :: word
    integer :: i, i_old
    integer :: max_length = 0

    words_count = 0
    i = 1
    i_old = 1
    call str_parse_next_token(string, i, word)
    do while (abs(i - i_old) > 1)
        if (len(word) > max_length) max_length = len(word)
        words_count = words_count + 1
        i_old = i
        call str_parse_next_token(string, i, word)
    end do

    if (ALLOCATED(words)) deallocate(words)
    allocate(character(len=max_length) :: words(1:words_count * max_length))

    words_count = 0
    i = 1
    i_old = 1
    call str_parse_next_token(string, i, word)
    do while (abs(i - i_old) > 1)
        words_count = words_count + 1
        words(words_count) = word
        i_old = i
        call str_parse_next_token(string, i, word)
    end do

end function Py_split

subroutine getblk(iu, itoken, ftoken, ldim, lines, nl)
    !     -------------------------------------------------

    ! getblk: get string array lines(1:nl) from open file with unit iu
    ! _between_ the line with the initial token itoken and line with the
    ! final token ftoken.
    ! Returns nl=0 if either the tokens have not been found or
    ! no lines are in between.
    !
    integer, intent(in) :: iu
    character(len = *), intent(in) :: itoken, ftoken
    integer, intent(in) :: ldim
    character(len = *), intent(inout) :: lines(ldim)
    integer, intent(out) :: nl
    integer k, io
    character*120 line

    nl = 0
    rewind(iu)
    do
        read(iu, '(A)', iostat = io) line
        k = index(line, itoken)
        if (k .gt. 0 .or. io .ne. 0) goto 101
    enddo
    101  continue
    if (io .eq. 0) then
        do
            read(iu, '(A)', iostat = io) line
            k = index(line, ftoken)
            if (k .gt. 0 .or. io .ne. 0) goto 201
            if (nl == ldim)&
                    call error("parselib:getblk: ldim too small")
            nl = nl + 1
            lines(nl) = line
        enddo
        201     continue
    endif

    return
end

!====================================================

!     --------------------------------------------------------------
subroutine getNextBlock(allLines, nla, idx, itoken, ftoken, ctoken, &
        ldim, lines, nl)
    !     --------------------------------------------------------------

    ! getNextBlock: get string array lines(1:nl) from AllLines(1:nla)
    ! _starting_ from the line with the initial token itoken and up to the
    ! line with the final token ftoken.
    ! lines starting in column 1 with the comment token ctoken are ignored
    ! start search from idx. change idx to line of final token ftoken +1
    ! Returns nl=0 if either the tokens have not been found or
    ! number of lines in block 'lines'
    !
    integer, intent(in) :: nla
    character(len = *), intent(in) :: allLines(nla)
    character(len = *), intent(in) :: itoken, ftoken, ctoken
    integer, intent(inout) :: idx
    integer, intent(in) :: ldim
    character(len = *), intent(inout) :: lines(ldim)
    integer, intent(out) :: nl

    integer i, k, n, idx0

    nl = 0
    idx0 = idx
    do i = idx, nla
        if (allLines(i)(1:1)==ctoken(1:1)) goto 100
        k = index(allLines(i), itoken)
        if (k .gt. 0) goto 101
        100     continue
    enddo
    101  continue

    if (k .gt. 0) then
        n = i
        nl = 1
        do i = n, nla
            k = index(allLines(i), ftoken)
            if (nl > ldim) call error("getNextBlock: wrong dimension")
            lines(nl) = allLines(i)
            idx = i + 1
            if (k .gt. 0) goto 201
            nl = nl + 1
        enddo
        nl = 0
        idx = idx0
        201     continue
    endif

    return
end

!====================================================

!     ----------------------------------------
subroutine getdblf(iu, target, value, iflag)
    !     ----------------------------------------

    ! getdblf: get double precision value after string target in open file
    ! with unit iu
    ! iflag: 0 if target found, 1 if not
    !
    integer i, iu, k, kf, io, iflag
    real(r8) value
    character line*120, target*(*)

    rewind(iu)
    do i = 1, 1000
        read(iu, '(A)', iostat = io) line
        k = index(line, target)
        if (k .gt. 0 .or. io .ne. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                      ! target not found
    else
        iflag = 0                      ! target found
        k = k + len(target)
        kf = index(line(k:), ')')
        if (kf == 0) then
            kf = len(line)
        else
            kf = k + kf - 2
        endif
        read(line(k:kf), *) value         ! internal file used for conversion
    endif

    return
end

!============================================

!     ----------------------------------------------
subroutine getdbla(lines, nl, target, value, iflag)
    !     ----------------------------------------------

    ! getdbla: get double precision value after string target in string array
    ! lines(1:nl)
    ! iflag: 0 if target found, 1 if not
    !
    integer i, k, kf, nl, iflag
    real(r8) value
    character lines(nl)*(*), target*(*)

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                          ! target not found
    else
        iflag = 0                          ! target found
        k = k + len(target)
        kf = index(lines(i)(k:), ')')
        if (kf == 0) then
            kf = len(lines(i))
        else
            kf = k + kf - 2
        endif
        read(lines(i)(k:kf), *) value         ! internal file used for conversion
    endif

    return
end

!====================================================

!     ----------------------------------------------
subroutine getdblarra(lines, nl, target, value, iflag)
    !     ----------------------------------------------

    ! getdbla: get double precision array after string target in string array
    ! lines(1:nl)
    ! iflag: 0 if target found, 1 if not
    !
    integer i, k, kf, nl, iflag, n, io
    character str*80
    character lines(nl)*(*), target*(*)
    character(len=:), allocatable :: token
    real(r8), allocatable :: value(:)

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                          ! target not found
    else
        iflag = 0                          ! target found
        k = k + len(target)
        kf = index(lines(i)(k:), ',')
        if (kf == 0) kf = index(lines(i)(k:), ')')
        if (kf == 0) then
            kf = len(lines(i))
        else
            kf = k + kf - 2
        endif
        read(lines(i)(k:kf), '(a)', iostat = io) str
        if (io /= 0) then
            read(lines(i)(k:), '(a)', iostat = io) str
            if (io /= 0) call error("cannot parse string")
        endif
        i = 1
        n = 0
        do
            call str_parse_next_token(str, i, token)
            if (token == '') exit
            if (.not. str_is_real(token)) call error("cannot parse string")
            n = n + 1
        end do
        do i = 1, LEN_TRIM(str)
            if (str(i : i) == ';') str = str(: i - 1) // ' ' // str(i + 1 :)
        end do
        allocate(value(n))
        read(str, *) value         ! internal file used for conversion
    endif

    return
end

!====================================================

!     -----------------------------------------
subroutine getintf(iu, target, value, iflag)
    !     -----------------------------------------

    ! getintf: get integer value after string target in open file
    ! with unit iu
    ! iflag: 0 if target found, 1 if not
    !
    integer i, iu, k, kf, io, value, iflag
    character line*120, target*(*)

    rewind(iu)
    do i = 1, 1000
        read(iu, '(A)', iostat = io) line
        k = index(line, target)
        if (k .gt. 0 .or. io .ne. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                      ! target not found
    else
        iflag = 0                      ! target found
        k = k + len(target)
        kf = index(line(k:), ')')
        if (kf == 0) then
            kf = len(line)
        else
            kf = k + kf - 2
        endif
        read(line(k:kf), *) value         ! internal file used for conversion
    endif

    return
end

!====================================================

!     ------------------------------------------
subroutine getint8f(iu, target, value, iflag)
    !     ------------------------------------------

    ! getintf: get integer value after string target in open file
    ! with unit iu
    ! iflag: 0 if target found, 1 if not
    !
    integer(i8) value
    integer i, iu, k, kf, io, iflag, n
    character line*120, target*(*)
    character str*80

    rewind(iu)
    do i = 1, 1000
        read(iu, '(A)', iostat = io) line
        k = index(line, target)
        if (k .gt. 0 .or. io .ne. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                      ! target not found
    else
        iflag = 0                      ! target found
        k = k + len(target)
        kf = index(line(k:), ')')
        if (kf == 0) then
            kf = len(line)
        else
            kf = k + kf - 2
        endif
        read(line(k:kf), *, iostat = io) str
        if (io /= 0) then
            read(line(k:), *, iostat = io) str
            if (io /= 0) call error("cannot parse string")
        endif
        n = len(trim(str))
        if (str(n:n)=='k') then
            read(str(:n - 1), *) value
            value = value * 1000
        else if (str(n:n)=='M') then
            read(str(:n - 1), *) value
            value = value * 1000000
        else
            read(str, *) value
        endif
    endif

    return
end

!============================================

!     -----------------------------------------------
subroutine getinta(lines, nl, target, value, iflag)
    !     -----------------------------------------------

    ! getinta: get integer value after string target in string array
    ! lines(1:nl)
    ! iflag: 0 if target found, 1 if not
    !
    integer i, k, kf, nl, value, iflag, n, io
    character lines(nl)*(*), target*(*)
    character str*80

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                        ! target not found
    else
        iflag = 0                        ! target found
        k = k + len(target)
        kf = index(lines(i)(k:), ')')
        if (kf == 0) then
            kf = len(lines(i))
        else
            kf = k + kf - 2
        endif
        read(lines(i)(k:kf), *, iostat = io) str
        if (io /= 0) then
            read(lines(i)(k:), *, iostat = io) str
            if (io /= 0) call error("cannot parse string")
        endif
        n = len(trim(str))
        if (str(n:n)=='k') then
            read(str(:n - 1), *) value
            value = value * 1000
        else if (str(n:n)=='M') then
            read(str(:n - 1), *) value
            value = value * 1000000
        else
            read(str, *) value
        endif
    endif

    return
end

!============================================

!     -----------------------------------------------
subroutine getintarra(lines, nl, target, value, iflag)
    !     -----------------------------------------------

    ! getinta: get integer array after string target in string array
    ! lines(1:nl)
    ! iflag: 0 if target found, 1 if not
    !
    integer i, k, kf, nl, iflag, n, io
    character lines(nl)*(*), target*(*)
    character str*80
    character(len=:), allocatable :: token
    integer, allocatable :: value(:)

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                        ! target not found
    else
        iflag = 0                        ! target found
        k = k + len(target)
        kf = index(lines(i)(k:), ',')
        if (kf == 0) kf = index(lines(i)(k:), ')')
        if (kf == 0) then
            kf = len(lines(i))
        else
            kf = k + kf - 2
        endif
        read(lines(i)(k:kf), '(a)', iostat = io) str
        if (io /= 0) then
            read(lines(i)(k:), '(a)', iostat = io) str
            if (io /= 0) call error("cannot parse string")
        endif
        i = 1
        n = 0
        do
            call str_parse_next_token(str, i, token)
            if (token == '') exit
            if (.not. str_is_integer(token)) call error("cannot parse string")
            n = n + 1
        end do
        do i = 1, LEN_TRIM(str)
            if (str(i : i) == ';') str = str(: i - 1) // ' ' // str(i + 1 :)
        end do
        allocate(value(n))
        read(str, *) value
    endif

    return
end

!============================================

!     ------------------------------------------------
subroutine getint8a(lines, nl, target, value, iflag)
    !     ------------------------------------------------

    ! getint8a: get integer value after string target in string array
    ! lines(1:nl)
    ! iflag: 0 if target found, 1 if not
    !
    integer(i8) value
    integer i, k, kf, n, nl, iflag, io
    character lines(nl)*(*), target*(*)
    character str*80

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                        ! target not found
    else
        iflag = 0                        ! target found
        k = k + len(target)
        kf = index(lines(i)(k:), ')')
        if (kf == 0) then
            kf = len(lines(i))
        else
            kf = k + kf - 2
        endif
        read(lines(i)(k:kf), *, iostat = io) str
        if (io /= 0) then
            read(lines(i)(k:), *, iostat = io) str
            if (io /= 0) call error("cannot parse string")
        endif
        n = len(trim(str))
        if (str(n:n)=='k') then
            read(str(:n - 1), *) value
            value = value * 1000
        else if (str(n:n)=='M') then
            read(str(:n - 1), *) value
            value = value * 1000000
        else
            read(str, *) value
        endif
    endif

    return
end

!====================================================

!     ---------------------------------------
subroutine getstrf(iu, target, str, iflag)
    !     ---------------------------------------

    ! getstrf: get string str (with apostroph!) after string target in open file
    ! with unit iu
    ! iflag: 0 if target found, 1 if not
    !
    integer i, iu, k, kf, io, iflag
    character line*120, target*(*), str*(*)

    rewind(iu)
    do i = 1, 1000
        read(iu, '(A)', iostat = io) line
        k = index(line, target)
        if (k .gt. 0 .or. io .ne. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                      ! target not found
    else
        iflag = 0                      ! target found
        k = k + len(target)
        read(line(k:), *, iostat = io) str         ! internal file used for conversion
        if (io /= 0) then
            kf = index(line(k:), ')')
            if (kf == 0) then
                kf = len(line)
            else
                kf = k + kf - 2
            endif
            read(line(k:kf), *) str         ! internal file used for conversion
        endif
    endif

    return
end

!============================================

!     ---------------------------------------------
subroutine getstra(lines, nl, target, str, iflag)
    !     ---------------------------------------------

    ! getstra: get string str (with apostrophs) after string target in string array
    ! lines(1:nl)
    ! iflag: 0 if target found, 1 if not
    !
    integer i, k, kf, nl, iflag, io
    character lines(nl)*(*), target*(*), str*(*)

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                        ! target not found
    else
        iflag = 0                        ! target found
        k = k + len(target)
        kf = index(lines(i)(k:), ')')
        if (kf == 0) then
            kf = len(lines(i))
        else
            kf = k + kf - 2
        endif
        read(lines(i)(k:kf), *, iostat = io) str
        if (io /= 0) then
            read(lines(i)(k:), *, iostat = io) str
            if (io /= 0) call error("cannot parse string")
        endif
    endif

    return
end

!====================================================

!     -----------------------------------------
subroutine getlogf(iu, target, value, iflag)
    !     -----------------------------------------

    ! getlogf: get logical value after string target in open file
    ! with unit iu
    ! iflag: 0 if target found, 1 if not
    !
    integer i, iu, k, kf, io, iflag
    character line*120, target*(*)
    logical value

    rewind(iu)
    do i = 1, 1000
        read(iu, '(A)', iostat = io) line
        k = index(line, target)
        if (k .gt. 0 .or. io .ne. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                      ! target not found
    else
        iflag = 0                      ! target found
        k = k + len(target)
        kf = index(line, ')')
        if (kf == 0) then
            kf = len(line)
        else
            kf = k + kf - 2
        endif
        read(line(k:kf), *) value         ! internal file used for conversion
    endif

    return
end

!============================================

!     -----------------------------------------------
subroutine getloga(lines, nl, target, value, iflag)
    !     -----------------------------------------------

    ! getloga: get logical value after string target in string array
    ! lines(1:nl)
    ! iflag: 0 if target found, 1 if not
    !
    integer i, k, kf, nl, iflag
    character lines(nl)*(*), target*(*)
    logical value

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        iflag = 1                        ! target not found
    else
        iflag = 0                        ! target found
        k = k + len(target)
        kf = index(lines(i)(k:), ')')
        if (kf == 0) then
            kf = len(lines(i))
        else
            kf = k + kf - 2
        endif
        read(lines(i)(k:kf), *) value       ! internal file used for conversion
    endif

    return
end


!====================================================

!     ---------------------------------
logical function findf(iu, target)
    !     ---------------------------------

    ! findf: find string target in open file
    ! with unit iu
    !
    integer iu, k, io
    character line*120, target*(*)

    rewind(iu)
    !      do i=1,1000
    !        read(iu,'(A)',iostat=io) line
    !        k = index(line,target)
    !        if ( k .gt. 0 .or. io .ne. 0) goto 101
    !      enddo
    ! 101  continue

    !cc modified 24.02.03 by CD

    io = 0
    do while (io.eq.0)
        read(iu, '(A)', iostat = io) line
        k = index(line, target)
        if (k .gt. 0) io = 1
    enddo

    !cc end CD

    if (k .eq. 0) then
        findf = .false.                ! target not found
    else
        findf = .true.                 ! target found
    endif

    return
end

!============================================

!     ---------------------------------------
logical function finda(lines, nl, target)
    !     ---------------------------------------

    ! finda: find string target in string array
    ! lines(1:nl)
    !
    integer i, k, nl
    character lines(nl)*(*), target*(*)

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        finda = .false.                ! target not found
    else
        finda = .true.                 ! target found
    endif

    return
end

!============================================

!     ----------------------------------------
integer function ifinda(lines, nl, target)
    !     ----------------------------------------

    ! ifinda: find string target in string array
    ! lines(1:nl). return line number if found,
    ! else zero
    !
    integer i, k, nl
    character lines(nl)*(*), target*(*)

    do i = 1, nl
        k = index(lines(i), target)
        if (k .gt. 0) goto 101
    enddo
    101  continue

    if (k .eq. 0) then
        ifinda = 0                ! target not found
    else
        ifinda = i                 ! target found
    endif

    return
end


!============================================

!     -----------------------------------------------
subroutine replaceEntry(line, entryIdx, newEntry)
    !     -----------------------------------------------

    integer entryIdx
    character line*(*), newEntry*(*)
    character str*40
    integer nnew, io, nold, offset

    nnew = len_trim(newEntry)
    read(line(entryIdx:), *, iostat = io) str
    nold = len_trim(str)
    if (str(nold:nold)==')') nold = nold - 1

    offset = nnew - nold
    call shiftTail(line, entryIdx, offset)
    line(entryIdx:entryIdx + nnew - 1) = newEntry(1:nnew)
    return
end


!============================================

!     -------------------------------------
subroutine shiftTail(line, idx, offset)
    !     -------------------------------------

    !     ! shift tail of line starting from idx by 'offset' characters
    !     ! offset > 0: insert 'offset' spaces at 'idx' position
    !     ! offset < 0: remove 'offset' characters starting from 'idx'

    character line*120
    integer idx, offset
    integer newlen, oldlen, i, ii

    oldlen = len_trim(line)
    newlen = oldlen + offset

    if (newlen > 120) then
        newlen = 120
        oldlen = newlen - offset
    endif

    if (offset > 0) then ! shift to right
        do i = oldlen, idx, -1
            ii = i + offset
            line(ii:ii) = line(i:i)
        enddo
        do i = idx, idx + offset - 1
            line(i:i) = ' '
        enddo
    else if (offset < 0) then ! shift to left
        do i = idx, newlen
            ii = i - offset
            line(i:i) = line(ii:ii)
        enddo
        do i = newlen + 1, oldlen
            line(i:i) = ' '
        enddo
    endif

    return
end

end module parsing_m