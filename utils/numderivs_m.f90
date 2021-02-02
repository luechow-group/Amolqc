! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module numderivs_m

   use kinds_m, only: r8
   use error_m

   implicit none

contains

   ! functions assume h spaced function values with increasing x. f(x0) is the mid array value!
   function numFirstDerivative3Point(f, h) result(res)
      real(r8), intent(in)  :: f(:)
      real(r8), intent(in)  :: h
      real(r8) :: res

      res =  (f(3) - f(1)) / (2.0_r8 * h)
   end function numFirstDerivative3Point

   function numFirstDerivative5Point(f, h) result(res)
      real(r8), intent(in)  :: f(:)
      real(r8), intent(in)  :: h
      real(r8) :: res

      res =   ( 8.0_r8 * (f(4) - f(2)) - (f(5) - f(1)) ) / (12.0_r8 * h)
   end function numFirstDerivative5Point

   function numFirstDerivative7Point(f, h) result(res)
      real(r8), intent(in)  :: f(:)
      real(r8), intent(in)  :: h
      real(r8) :: res

      res = ( 45.0_r8 * (f(5) - f(3)) - 9.0_r8 * (f(6) - f(2)) + (f(7) - f(1)) ) / (60.0_r8 * h)
   end function numFirstDerivative7Point

   function numFirstDerivative9Point(f, h) result(res)
      real(r8), intent(in)  :: f(:)
      real(r8), intent(in)  :: h
      real(r8) :: res

      res = (672.0_r8 * (f(4) - f(6)) - 168.0_r8 * (f(3) - f(7)) + &
              32.0_r8 * (f(2) - f(8)) - 3.0_r8 * (f(1) - f(9))) / (840.0_r8 * h)
   end function numFirstDerivative9Point


   function numSecondDerivative3Point(f, h) result(res)
      real(r8), intent(in)  :: f(:)
      real(r8), intent(in)  :: h
      real(r8) :: res

      res =  (f(1) - 2.0_r8*f(2) + f(3)) / h**2
   end function numSecondDerivative3Point

   function numSecondDerivative5Point(f, h) result(res)
      real(r8), intent(in)  :: f(:)
      real(r8), intent(in)  :: h
      real(r8) :: res

      res = ( -30.0_r8 * f(3) + 16.0_r8 * (f(2)+f(4)) - (f(1)+f(5)) ) / (12.0_r8 * h**2)
   end function numSecondDerivative5Point

   function numSecondDerivative7Point(f, h) result(res)
      real(r8), intent(in)  :: f(:)
      real(r8), intent(in)  :: h
      real(r8) :: res

      res = ( -490.0_r8 * f(4) + 270.0_r8 * (f(3) + f(5)) - 27.0_r8 * (f(2) + f(6)) + 2.0_r8 * (f(1) + f(7)) ) / (180.0_r8 * h**2)
   end function numSecondDerivative7Point

   function numSecondDerivative9Point(f, h) result(res)
      real(r8), intent(in)  :: f(:)
      real(r8), intent(in)  :: h
      real(r8) :: res

      res =  (-71750.0_r8 * f(5) + 40320.0_r8 * (f(4) + f(6))&
              - 5040.0_r8 * (f(3) + f(7))&
              + 640.0_r8 * (f(2) + f(9))&
              - 45.0_r8 * (f(1) + f(9))) / (25200.0_r8 * h**2)
   end function numSecondDerivative9Point

   subroutine nevilleExtrapolation(t0, extrapolatedValue, estimatedError, iterations)
      real(r8), intent(in)           :: t0(0:)
      real(r8), intent(out)          :: extrapolatedValue
      real(r8), optional, intent(out):: estimatedError
      integer, optional, intent(out) :: iterations
      real(r8) :: oldValue, oldError, eps, bestValue, bestError
      real(r8), allocatable :: T(:,:)
      integer  :: i, j, bestIter

      allocate(T(0 : size(t0) - 1, 0 : size(t0) - 1))
      oldValue = huge(1.0_r8)
      bestError = huge(1.0_r8)
      bestIter = size(t0) - 1
      if (present(iterations)) iterations = size(t0) - 1
      do i = 0, size(t0) - 1
         T(i, 0) = t0(i)
         do j = 1, i
            T(i, j) = ( 4**j * T(i, j-1) - T(i-1, j-1) ) / (4**j - 1)
         end do
         eps = abs(T(i,i) - oldValue)
         !!!print*, i, T(i,i), eps
         if (eps < bestError) then
            bestValue = T(i,i)
            bestError = eps
            bestIter = i
         end if
         oldValue = T(i,i)
      end do
      if (present(estimatedError)) estimatedError = bestError
      if (present(iterations)) iterations = bestIter
      extrapolatedValue = bestValue

   end subroutine nevilleExtrapolation

end module numderivs_m
