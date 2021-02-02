! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module myfctn_tm
   use kinds_m, only: r8
   use error_m
   use myfctn_module
   use fctn_module
   implicit none
contains
   subroutine myfctn_test()

      real(r8), allocatable       :: x(:), g(:), g1(:), H(:,:)
      real(r8)                    :: f, gnorm, gerror, result
      integer                   :: n, i, j
      type(myFunction_t), target :: fg
      type(mysingFunction_t), target :: fgs
      class(Function_t), pointer :: fg_p => null()

      !!! test of myfctn

      n = 3
      allocate(x(n), g(n), H(n,n))
      x = [1, 2, 3]

      fg = myFunction_t(a=2.0d0)
      result = fg%eval(x)
      call assertEqualRelative(result, 196.d0, 1.d-14, "myfctn_mt: eval function value test failed")

      call fg%eval_fg(x, f, g)
      call assertEqualRelative(f, 196.d0, 1.d-14, "myfctn_mt: eval_fg function value test failed")

      g = g - [8.d0, 64.d0, 216.d0]
      gnorm = dot_product(g, g)
      call assertEqualRelative(gnorm, 0.d0, 1.d-13, "myfctn_mt: grad test failed")

      fg_p => fg
      result = fg_p%eval(x)
      call assertEqualRelative(result, 196.d0, 1.d-14, "myfctn_mt: eval function value test failed")

      call fg_p%eval_fg(x, f, g)
      call assertEqualRelative(f, 196.d0, 1.d-14, "myfctn_mt: base function test failed")

      g = g - [8.d0, 64.d0, 216.d0]
      gnorm = dot_product(g, g)
      call assertEqualAbsolute(gnorm, 0.d0, 1.d-13, "myfctn_mt: base pointer grad test failed")

      call fg%eval_fg_num(x, f, g)
      g = g - [8.d0, 64.d0, 216.d0]
      gnorm = dot_product(g, g)
      call assertEqualAbsolute(gnorm, 0.d0, 1.d-6, "myfctn_mt: numerical grad test failed")

      call fg%eval_fgh_num(x, f, g, H)
      do i=1,n
         do j=1,n
            if (i==j) then
               call assertEqualAbsolute(H(i,j), 12*2*x(i)**2, 1.d-6, "myfctn_mt: numerical hessian test failed")
            else
               call assertEqualAbsolute(H(i,j), 0.d0, 1.d-6, "myfctn_mt: numerical hessian test failed")
            end if
         end do
      end do

      g = g - [8.d0, 64.d0, 216.d0]
      gnorm = dot_product(g, g)

      call fg_p%eval_fg_num(x, f, g)
      g = g - [8.d0, 64.d0, 216.d0]
      gnorm = dot_product(g, g)
      call assertEqualAbsolute(gnorm, 0.d0, 1.d-5, "myfctn_mt: base pointer numerical grad test failed")


      !!! test of mysingfctn

      n = 6
      deallocate(x, g)
      allocate(x(n), g(n), g1(n))
      x = [3.d0, -2.d0, 0.1d0, 0.d0, 4.d0, 1.1d0]

      fgs = mysingFunction_t(a=2.0d0, b=2.0d0)
      fg_p => fgs

      call fg_p%eval_fg(x, f, g)
      call assertEqualRelative(f, 73.423725770523959d0, 1.d-14, "myfctn_mt: singfctn eval_fg function value test failed")
      call assertEqualRelative(sum(g), 4.5290299503941833d0, 1.d-14, "myfctn_mt: singfctn eval_fg grad sum test failed")

      call fg_p%set_numerical_denominator(h=1.e-7_r8)
      call fg_p%eval_fg_num(x, f, g1)

      gerror = sum(abs(g1 - g))
      call assertEqualAbsolute(gerror, 0.d0, 1.d-6, "myfctn_mt: singfctn numerical grad test failed")

   end subroutine myfctn_test
end module myfctn_tm