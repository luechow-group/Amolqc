! Copyright (C) 2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module plateau_m
   use kinds_m, only : r8, i4
   use genericFilter_m, only: Generic_filter, Std_deviation
   implicit none
   private

   public :: Find_plateau

contains
   function Find_plateau(data) result(plateau)
      real(r8), intent(in) :: data(:)
      logical :: plateau(SIZE(data))

      real(r8) :: filtered_data(SIZE(data))
      integer :: i, a, b, startA, startB

      ! standard deviation filter
      filtered_data = Generic_filter(data, Std_deviation, 3)

      ! normalization
      filtered_data = filtered_data / ( SUM(filtered_data)/SIZE(filtered_data) )

      ! get array of logicals
      plateau = filtered_data < 0.5_r8

      ! remove possible starting plateau
      i = 1
      do while (plateau(i))
         plateau(i) = .false.
         i = i + 1
      end do

      ! reduce to largest plateau, if several are found
      a = 0
      b = 0
      startA = 0
      startB = 0
      do i = 1, SIZE(plateau)
         if (plateau(i)) then
            if (a == 0) then
               startA = i
            end if
            a = a + 1
         else
            a = 0
         end if

         if (a > b) then
            b = a
            startB = startA
         end if
      end do

      if (b /= 0) then
         plateau = .false.
         do i = startB, startB + b - 1
            plateau(i) = .true.
         end do
      end if

   end function Find_plateau
end module plateau_m
