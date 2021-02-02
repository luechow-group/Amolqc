! Copyright (C) 2018 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module compilerStrings_m
   use, intrinsic :: iso_fortran_env, only: compiler_options, compiler_version
   use parsing_m, only: Py_split

   implicit none
   private

   public :: getCompilerOptions, getCompilerVersion

contains
   function getCompilerOptions() result(options)
      character(len=:), allocatable :: options
      character(len=:), allocatable :: words(:)
      integer :: i

      options = compiler_options()

      words = Py_split(options)

      options = ''

      i = 1
      do while (i < SIZE(words)/LEN(words))
         if (words(i) == '-I' &
        .or. words(i) == '-J' &
        .or. words(i) == '-D' &
        .or. words(i) == '-o' &
        .or. words(i) == '-module') then
            i = i + 1
         else if ( .not. (index(words(i), '.o') /= 0 &
                     .or. index(words(i), '-I') == 1 &
                     .or. index(words(i), '-D') == 1 &
                     .or. words(i) == '-c')) then
            if (options /= '') options = options // ' '
            options = options // trim(words(i))
         end if
         i = i + 1
      end do
   end function getCompilerOptions

   function getCompilerVersion() result(version)
      character(len=:), allocatable :: version

      version = compiler_version()

   end function getCompilerVersion

end module compilerStrings_m
