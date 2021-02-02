! Copyright (C) 2018-2019 Leonard Reuter
!
! SPDX-License-Identifier: GPL-3.0-or-later

module genericFilter_m
   use kinds_m, only : r8, i4
   use error_m, only: assert
   implicit none
   private

   public :: Std_deviation, Generic_filter
contains
   function Std_deviation(data, data_size) result(std_dev)
      integer(i4), intent(in) :: data_size
      real(r8), intent(in) :: data(data_size)
      real(r8) :: std_dev

      real(r8) :: mean
      integer :: i

      std_dev = 0._r8

      mean = SUM(data)/SIZE(data)
      do i=1, SIZE(data)
         std_dev = std_dev + (mean - data(i))**2
      end do

      std_dev = SQRT(std_dev / SIZE(data))
   end function Std_deviation

   ! implementation of pythons scipy.ndimage.filters.generic_filter
   function Generic_filter(data, Filter_procedure, filter_size) result(filtered_data)
      real(r8), intent(in) :: data(:)
      integer(i4), intent(in) :: filter_size
      real(r8), external :: Filter_procedure

      real(r8) :: filtered_data(SIZE(data))

      real(r8) :: temp_array(filter_size)
      integer(i4) :: i, l_range, r_range

      call assert(filter_size > 0, 'Generic_filter: filter_size is < 1.')
      call assert(filter_size <= SIZE(data), 'Generic_filter: filter_size is larger than data size.')

      temp_array = 0._r8

      filtered_data(SIZE(data)) = 0._r8

      ! always: filter_size = l_range + r_range
      ! if filter_size is odd: l_range = r_range
      ! else: l_range = r_range + 1
      l_range = filter_size / 2  ! 0, 1, 1, 2, 2 ...
      r_range = (filter_size - 1) / 2 ! 0, 0, 1, 1, 2, 2 ...

      do i = 1, SIZE(data)
         if (i < l_range + 1) then
            ! values for the filter are mirrored on the data border, with first value used twice
            temp_array(1:i + r_range) = data(1:i + r_range)
            temp_array(i + r_range + 1:filter_size) = data(1:filter_size - (i + r_range) )
            filtered_data(i) = Filter_procedure(temp_array, filter_size)
         else if (size(data) - i < r_range) then
            ! values for the filter are mirrored on the data border, with last value used twice
            temp_array(1:size(data) - i + 1 + l_range) = data(i - l_range:SIZE(data))
            temp_array(size(data) - i + 2 + l_range:filter_size) = data(2 * SIZE(data) - filter_size - i + 2 + l_range:SIZE(data))
            filtered_data(i) = Filter_procedure(temp_array, filter_size)
         else
            filtered_data(i) = Filter_procedure( data(i - l_range: i + r_range) , filter_size)
         end if
      end do
   end function Generic_filter

end module genericFilter_m