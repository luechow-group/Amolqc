! Copyright (C) 2019 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later!


! F90 module cspline_m. contains array of 1D cubic splines over [0:infinity[
!
! $Id: cspline_m.f,v 1.2 2008/02/07 17:00:53 luechow Exp $
!
! contains:
!    csplinit(...)              :  allocates arrays, set max # of spline points 
!    csplfinal()                :  deallocates arrays
!    csplxarray(np)             :  constructs array of spline knots
!                               :  (corresponds to spline evaluation on
!                               :   csplint)
!    cspline(ispl,...)          :  construct 'ispl'-th cubic spline, 
!                                  using 'n' (x,y) pairs with x points 
!                                  stored in array 'csplx(n)'
!    csplint(ispl,x)            :  evaluate 'ispl'-th spline at x
!
!
!     a cubic spline object is defined by number of points (knots) 
!     and the arrays for y and y''.
!
!     a static array of cubic spline objects is defined in the header file,
!     each objects is accessed by an index number (ispl).
!
!     equispaced values for x in [0,1[ mapped to [0,inf[ 
!     with x -> alpha*x/(1-x). 
!     More precisely:
!      h = 1.d0/(n-1)                         
!      do i=1,n-1
!         x(i) = alpha*(i-1)*h / (1.d0 - (i-1)*h))
!      enddo
!      x(n) = 1.d300                     ! better system parameter 'HUGE'
!
!     The alpha parameter (csalpha) scales the mapping.
!
!     The smaller alpha the closer are the spline knots for small x
!     (and the farther apart for large x) which is important for the 
!     exponential decay.
!
! this version saves the actual coefficients a_k,b_k,c_k,d_k 
! for the cubic polynomial to get faster spline interpolation.
!
!=========================================================

MODULE cspline_m

   use kinds_m, only : r8
   use error_m
   implicit none

   integer :: csnpmax = 0                      ! max. # of spline points (knots)
   integer :: csnsplmax = 0                    ! max. # of spline objects
   real(r8) :: csalpha = 1                      ! alpha param of mapping function
   integer :: csplnpnt = 0                     ! # of spline points

   real(r8), allocatable :: csplx(:)       ! array of x-values for all splines
   real(r8), allocatable :: cspla(:, :),&  ! array of A coefficients&
   csplb(:, :),&     ! array of B coefficients&
   csplc(:, :),&     ! array of C coefficients&
   cspld(:, :)       ! array of D coefficients


CONTAINS

   SUBROUTINE csplinit(nspl, np, alpha)
      ! csplinit allocates the spline arrays.

      integer nspl     ! (max) # spline functions
      integer np       ! (max) # spline points
      real(r8) alpha   ! alpha parameter for mapping [0:1] -> [0:infinity[
      integer ierr

      csnsplmax = nspl
      csnpmax = np
      csalpha = alpha

      allocate(cspla(csnsplmax, csnpmax), csplb(csnsplmax, csnpmax), &
             csplc(csnsplmax, csnpmax), cspld(csnsplmax, csnpmax), &
             csplx(csnpmax), stat = ierr)
      if (ierr /= 0) then
         if (allocated(cspla)) then
             call error("(csplinit): cspl arrays already allocated")
         else
             call error("(csplinit):allocation failed")
         endif
      endif

   END SUBROUTINE csplinit


   SUBROUTINE csplfinal()
      integer ierr

      deallocate(cspla, csplb, csplc, cspld, csplx, stat = ierr)
      if (ierr /= 0) then
         call error("csplfinal:deallocation failed")
      endif

   END SUBROUTINE csplfinal

   SUBROUTINE csplxarray(np)
      ! set the array 'csplx(j)' of np spline points
      ! for all functions (mapping of [0,1[ to [0,infty[ )
      integer np              ! # of spline knots
      integer i
      real(r8) h

      h = 1.d0 / (np - 1)
      do i = 1, np - 1
         csplx(i) = csalpha * (i - 1) * h / (1.d0 - (i - 1) * h)
      enddo
      csplx(np) = huge(1.d0)          ! "infinity"
      csplnpnt = np

   END SUBROUTINE csplxarray

   SUBROUTINE cspline(ispl, y, n, y1_x1, y11_x1)
      ! requires: # of spline points 'csplnpnt' and the array csplx(j),
      ! i.e. the points x_j, j=1,..,csplnpnt to be set.
      integer ispl        ! index of spline object
      integer n           ! size of array y, # of knots
      real(r8) y(n)       ! y(x_i), function values at spline points
      real(r8) y1_x1, y11_x1  ! y'(x_1) and y''(x_1)
      integer i, k
      real(r8) denom, qn, tmp, un
      real(r8) deltax, deltay, deltay2
      real(r8) y2(n), work(n)

      if (ispl < 1 .or. ispl > csnsplmax) call error("(cspline): illegal spline index")

      if (n /= csplnpnt) call error('(cspline): inconsistent use of cspline')

      y2(1) = y11_x1
      work(1) = (3.d0 / ( csplx(2) - csplx(1) ) ) * ( ( y(2) - y(1) )  &
              / ( csplx(2) - csplx(1) ) - y1_x1 )
      do i = 2, n - 1
         tmp = ( csplx(i) - csplx(i - 1) ) / ( csplx(i + 1) - csplx(i - 1) )
         denom = tmp * y2(i - 1) + 2.d0
         y2(i) = (tmp - 1.d0) / denom
         work(i) = ( 6.d0 * ( ( y(i + 1) - y(i) ) / ( csplx(i + 1) - csplx(i) ) -  &
                              ( y(i) - y(i - 1) ) / ( csplx(i) - csplx(i - 1) ) ) / &
                            ( csplx(i + 1) - csplx(i - 1) ) &
                     - tmp * work(i - 1) ) / denom
      end do
      tmp = -6.d0 * ( y(n) - y(n - 1) ) / ( csplx(n) - csplx(n - 1) )**2
      y2(n) = (tmp - work(n - 1)) / (y2(n - 1) + 2.d0)
      ! backsubstitution
      do k = n - 1, 1, -1
         y2(k) = y2(k) * y2(k + 1) + work(k)
      end do

      ! construct coefficients of cubic spline polynomial and save as
      ! spline object 'ispl'
      do k = 1, n - 1
         cspla(ispl, k) = y(k)

         deltax = csplx(k + 1) - csplx(k)
         deltay = y(k + 1) - y(k)
         deltay2 = y2(k + 1) - y2(k)
         csplb(ispl, k) = deltay / deltax - deltax * (deltay2 / 6.d0 + y2(k) / 2.d0)
         csplc(ispl, k) = y2(k) / 2.d0
         cspld(ispl, k) = deltay2 / (6.d0 * deltax)
      end do

   END SUBROUTINE cspline

   real(r8) FUNCTION csplint(ispl, x)
      ! requires: previous creation of spline points csplx()
      ! and # of spline points csplnpt
      ! requires: previous call to cspline to set
      ! coeff-arrays cspla,csplb,csplc,cspld

      ! cspline interpolation: assumes the knots at
      !     x_i = alpha*(i-1)*h / (1 - (i-1)*h), i=1,..,n-1, x_n = HUGE (1.d300)
      ! (calculated with 'csplxarray(np)'
      ! inversion gives x in [x_k,x_k+1] with k = int((n-1)*x/(x+alpha))
      !

      integer ispl
      real(r8) x
      integer j
      real(r8) dx

      j = (csplnpnt - 1) * x / (csalpha + x) + 1

      dx = x - csplx(j)
      csplint = cspla(ispl, j) + dx * (csplb(ispl, j) + dx * (csplc(ispl, j)&
             + dx * cspld(ispl, j)))
   END FUNCTION csplint

END MODULE cspline_m