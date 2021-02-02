! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module kinds_m
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

    implicit none
    private

    public :: i4, i8, r4, r8

    integer(int32), parameter :: i4 = int32
    integer(int32), parameter :: i8 = int64
    integer(int32), parameter :: r4 = real32
    integer(int32), parameter :: r8 = real64

end module kinds_m