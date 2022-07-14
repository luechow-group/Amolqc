! Copyright (C) 2017-2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module fPsi2_m
   use kinds_m, only: r8
   use error_m
   use elocData_m, only: elPhi, elU, elxDrift, elyDrift, elzDrift
   use eloc_m, only: eloc, ELOC_X_INF_ERR
   use eConfigs_m
   use fctn_module
   implicit none

   type, extends(Function_t) :: fctn_psi2
   contains
      procedure :: eval_fg => fpsi2_fg
      procedure :: eval => fpsi2
   end type fctn_psi2

contains

   subroutine fpsi2_fg(this, x, f, g, mask)
      class(fctn_psi2), intent(in)  :: this
      real(r8), intent(in) :: x(:)      ! electron vector (without elecs at cores)
      real(r8), intent(out) :: f      ! function value
      real(r8), intent(out) :: g(:)   ! gradient
      logical, intent(in), optional :: mask(SIZE(x))
      type(eConfigArray)  :: eca
      real(r8), allocatable :: xx(:), yy(:), zz(:)
      integer n, i, error_code

      if (asserts) then 
         call assert(mod(size(x),3) == 0 .and. size(x) == size(g), "fpsi2_fg: illegal size")
         call assert(all(abs(x)<huge(1.d0)), "fpsi2_fg: illegal x coords")
      end if

      n = size(x)/3
      allocate(xx(n), yy(n), zz(n))
      do i = 1, n
         xx(i) = x(3*i-2)
         yy(i) = x(3*i-1)
         zz(i) = x(3*i)
      end do

      call eConfigArray_new(eca, n, 1)
      call eConfigArray_set(eca, 1, xx, yy, zz)
      call eloc(0, eca, 'none', error_code)

      if (debug) then
         if (error_code == ELOC_X_INF_ERR) then
            call abortp("fpsi2_fg: eloc returns ELOC_X_INF_ERR")
         end if 
      end if

      f = - 2.d0*(log(abs(elPhi(1))) + elU(1))

      do i = 1, n
         g(3*i-2) = - 2.d0 * elxDrift(i, 1)
         g(3*i-1) = - 2.d0 * elyDrift(i, 1)
         g(3*i)   = - 2.d0 * elzDrift(i, 1)
      end do
      call eConfigArray_destroy(eca)
   end subroutine fpsi2_fg


   function fpsi2(this, x) result(f)
      class(fctn_psi2), intent(in)  :: this
      real(r8), intent(in) :: x(:)      ! electron vector (without elecs at cores)
      real(r8)             :: f      ! function value
      type(eConfigArray)  :: eca
      real(r8), allocatable :: xx(:), yy(:), zz(:)
      integer n, i, error_code

      if (asserts) then 
         call assert(mod(size(x),3) == 0, "fpsi2: illegal size")
         call assert(all(abs(x)<huge(1.d0)), "fpsi2_fg: illegal x coords")
      end if         

      n = size(x)/3
      allocate(xx(n), yy(n), zz(n))
      do i = 1, n
         xx(i) = x(3*i-2)
         yy(i) = x(3*i-1)
         zz(i) = x(3*i)
      end do

      call eConfigArray_new(eca, n, 1)
      call eConfigArray_set(eca, 1, xx, yy, zz)
      call eloc(0, eca, 'none', error_code)
      if (debug) then
         if (error_code == ELOC_X_INF_ERR) then
            call abortp("fpsi2: eloc returns ELOC_X_INF_ERR")
         end if 
      end if

      f = - 2.d0*(log(abs(elPhi(1))) + elU(1))

      call eConfigArray_destroy(eca)
   end function fpsi2

end module fPsi2_m
