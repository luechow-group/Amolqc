! Copyright (C) 2012 Alexander Sturm
!
! SPDX-License-Identifier: GPL-3.0-or-later

module aoMoTask_m

  use kinds_m, only: r8
  use wfData_m
  use aosData_m
  use mos_m
  use error_m
  use OMP_LIB

  implicit none

  integer, parameter :: BF_S = 1
  integer, parameter :: BF_P = 2
  integer, parameter :: BF_D = 3
  integer, parameter :: BF_F = 4

  ! type definitions
  type bf_t
    logical :: allocated = .false.
    integer :: ngto = 0
    integer :: atom = 0
    integer :: type = 0
    real(r8) :: aocut = 0d0
    real(r8), allocatable :: alpha(:), coeff(:), moc(:)
  end type

  integer, parameter :: BFS_DIM = 5
  type bfs_t
    integer :: snum = 0, pnum = 0, dnum = 0, fnum = 0
    type(bf_t) :: s(BFS_DIM), p(BFS_DIM), d(BFS_DIM), f(BFS_DIM)
  end type

  type elec_bf_t
    integer, allocatable, dimension(:) :: bfs
    integer :: cnt = 0
    logical :: allocated = .false.
  end type

contains

  subroutine bf_init(bf, ngto)
    type(bf_t), intent(inout) :: bf
    integer, intent(in) :: ngto

    integer :: istat

    call assert(.not. bf%allocated, "bf already initialized")
    bf%allocated = .true.
    bf%ngto = ngto

    allocate(bf%alpha(ngto), bf%coeff(ngto), bf%moc(10*norb), stat=istat)
    call assert(istat == 0, "BF instance allocation in aomo_calc_task failed")
  end subroutine

  subroutine bf_destroy(bf)
    type(bf_t), intent(inout) :: bf
    if (.not. bf%allocated) return

    deallocate(bf%alpha, bf%coeff)
    deallocate(bf%moc)
    bf%allocated = .false.
  end subroutine

  pure logical function bf_cutoff(bf, rai)
    type(bf_t), intent(in) :: bf
    real(r8), intent(in) :: rai

    bf_cutoff = .false.

    if (cutao) then
      if (rai .gt. bf%aocut) bf_cutoff = .true.
    endif
  end function

  subroutine bf_build(bf, bfi, moc, nf, cutao)
    type(bf_t), intent(inout) :: bf
    integer, intent(in) :: bfi, moc, nf
    logical, intent(in) :: cutao

    integer :: i, idx

    call bf_init(bf, ngto(bfi))

    if(bl(bfi) == 'F' .and. .not. (evfmt == 'gau' .or. evfmt=='mol') ) then
      ! re-sort into Gaussian order
      ! using pointers instead of copying is slower
      idx = moc
      do i = 0, 10*(norb-1), 10
        bf%moc(i+1:i+3) = cmoa(idx:idx+2) ! xxx, yyy, zzz
        bf%moc(i+4) = cmoa(idx+5)  ! yyx
        bf%moc(i+5) = cmoa(idx+3)  ! xxy
        bf%moc(i+6) = cmoa(idx+4)  ! xxz
        bf%moc(i+7) = cmoa(idx+7)  ! zzx
        bf%moc(i+8) = cmoa(idx+8)  ! zzy
        bf%moc(i+9) = cmoa(idx+6)  ! yyz
        bf%moc(i+10) = cmoa(idx+9) ! xyz
        idx = idx + 10
      enddo
    else
      bf%moc(1:nf) = cmoa(moc:moc+nf-1)
    endif

    if (cutao) bf%aocut = aocuts(bfi)
    bf%alpha(1:bf%ngto) = cntrctn(1, 1:bf%ngto, bfi)
    bf%coeff(1:bf%ngto) = cntrctn(2, 1:bf%ngto, bfi)
  end subroutine bf_build

  subroutine bfs_init(bfs)
    type(bfs_t), intent(inout) :: bfs
    bfs%snum = 0; bfs%pnum = 0; bfs%dnum = 0; bfs%fnum = 0
  end subroutine bfs_init

  subroutine bfs_destroy(bfs)
    type(bfs_t), intent(inout) :: bfs

    integer :: i

    do i = 1, BFS_DIM
      call bf_destroy(bfs%s(i))
      call bf_destroy(bfs%p(i))
      call bf_destroy(bfs%d(i))
      call bf_destroy(bfs%f(i))
    enddo
  end subroutine bfs_destroy

  subroutine elecbf_init(elecBF, nbasf)
    type(elec_bf_t), intent(inout) :: elecBF
    integer, intent(in) :: nbasf

    integer :: istat
    if(elecBF%allocated) return

    allocate(elecBF%bfs(nbasf), stat=istat)
    call assert(istat == 0, "elecBF initialization failed in aomo_calc_task")

    elecBF%allocated = .true.
  end subroutine elecbf_init

  subroutine elecbf_destroy(elecBF)
    type(elec_bf_t), intent(inout) :: elecBF
    if(.not. elecBF%allocated) return

    deallocate(elecBF%bfs)
    elecBF%allocated = .false.
  end subroutine elecbf_destroy

  subroutine aomo_calc_task(ie, x, y, z, rai)
    ! electron index, 0 for all electrons
    integer, intent(in)                 :: ie
    ! electron positions
    real(r8), dimension(:), intent(inout) :: x, y, z

    ! electron-nucleus distances
    real(r8), dimension(:,:) :: rai

    ! local variables
    integer :: a, i, i1, i2
    integer :: istat

    type(bfs_t), dimension(:), allocatable, target :: bfs

    i1 = 1
    i2 = ne

    if (ie /= 0) then
      i1 = ie
      i2 = ie
    endif

    mat(1:norb,i1:i2,1)   = 0d0
    mat1x(1:norb,i1:i2,1) = 0d0
    mat1y(1:norb,i1:i2,1) = 0d0
    mat1z(1:norb,i1:i2,1) = 0d0
    mat2(1:norb,i1:i2,1)  = 0d0


    ! necessary every call? (ie, do they change?)
    allocate(bfs(ncenter), stat=istat)
    call assert(istat == 0, "BF structure allocation in aomo_calc_task failed")
    call sortBFs(bfs)

    !$omp parallel default(none) &
    !$omp& shared(i1,i2,ncenter) &
    !$omp& shared(bfs,rai,x,y,z) &
    !$omp& private(a,i)
    !$omp do schedule(auto)
    do i = i1, i2
      do a = 1, ncenter
        call calcAOMO(a, bfs, i, rai(a, i), x(i), y(i), z(i))
      enddo
    enddo
    !$omp end do
    !$omp end parallel

    do a = 1, ncenter
      call bfs_destroy(bfs(a))
    enddo
    deallocate(bfs)
  end subroutine aomo_calc_task

  subroutine aomo_calc_task_pair(ie, x, y, z, rai)
    ! electron index, 0 for all electrons
    integer, intent(in)                 :: ie
    ! electron positions
    real(r8), dimension(:), intent(inout) :: x, y, z

    ! electron-nucleus distances
    real(r8), dimension(:,:) :: rai

    ! local variables
    integer :: i, i1, i2
    integer :: istat

    type(bf_t), dimension(:), allocatable, target :: bfs
    type(elec_bf_t), dimension(:), allocatable, target :: elecBFs

    i1 = 1
    i2 = ne

    if (ie /= 0) then
      i1 = ie
      i2 = ie
    endif

    mat(1:norb,i1:i2,1)   = 0d0
    mat1x(1:norb,i1:i2,1) = 0d0
    mat1y(1:norb,i1:i2,1) = 0d0
    mat1z(1:norb,i1:i2,1) = 0d0
    mat2(1:norb,i1:i2,1)  = 0d0

    ! necessary every call? (ie, do they change?)
    allocate(elecBFs(i1:i2), bfs(nbasf), stat=istat)
    call assert(istat == 0, "elecBF allocation failed in aomo_calc_task")
    call sortBFsElec(elecBFs, bfs, i1, i2, rai)

    !$omp parallel default(none) &
    !$omp& shared(i1,i2,ncenter) &
    !$omp& shared(bfs,rai,x,y,z) &
    !$omp& shared(elecBFs) &
    !$omp& private(i)
    !$omp do schedule(auto)
    do i = i1, i2
      call calcAOMOelec(bfs, elecBFs(i), i, rai, x(i), y(i), z(i))
    enddo
    !$omp end do
    !$omp end parallel

    do i = i1, i2
      call elecbf_destroy(elecBFs(i))
    enddo
    deallocate(elecBFs)
    do i = 1, nbasf
      call bf_destroy(bfs(i))
    enddo
    deallocate(bfs)
  end subroutine aomo_calc_task_pair

  subroutine sortBFs(bfs)
    type(bfs_t), intent(inout), target :: bfs(:)
    integer :: nf,moc,bfi,a
    type(bfs_t), pointer :: f
    type(bf_t), pointer :: bf

    ! sort by function type (SPDF) and atom center
    moc = 1
    do bfi = 1, nbasf
      a = bc(bfi)
      f => bfs(a)

      select case(bl(bfi))
      case('S')
        f%snum = f%snum + 1
        call assert(f%snum <= BFS_DIM, "aomo_calc_task: too many s functions")

        bf => f%s(f%snum)
        nf = norb
      case('P')
        f%pnum = f%pnum + 1
        call assert(f%pnum <= BFS_DIM, "aomo_calc_task: too many p functions")

        bf => f%p(f%pnum)
        nf = 3*norb
      case('D')
        f%dnum = f%dnum + 1
        call assert(f%dnum <= BFS_DIM, "aomo_calc_task: too many d functions")

        bf => f%d(f%dnum)
        nf = 6*norb
      case('F')
        f%fnum = f%fnum + 1
        call assert(f%fnum <= BFS_DIM, "aomo_calc_task: too many f functions")

        bf => f%f(f%fnum)
        nf = 10*norb
      case default
        call abortp('(getaos): wrong GTO')
        ! just to silence compiler warnings about uninitialized variables
        return
      end select

      call bf_build(bf, bfi, moc, nf, cutao)
      moc = moc + nf
    enddo
  end subroutine sortBFs

  subroutine calcAOMO(a, bfs, i, rr, xx, yy, zz)
    type(bfs_t), intent(in), target :: bfs(:)
    integer, intent(in) :: a, i
    real(r8), intent(in) :: rr, xx, yy, zz

    real(r8) :: r2, dx, dy, dz

    type(bfs_t), pointer :: f
    type(bf_t), pointer :: bf

    integer :: si, pi, di, fi

    f => bfs(a)
    r2 = rr*rr
    dx = xx-atoms(a)%cx
    dy = yy-atoms(a)%cy
    dz = zz-atoms(a)%cz

    do si = 1, f%snum
      bf => f%s(si)

      if (bf_cutoff(bf, rr)) cycle

      call s_calc(bf, i, r2, dx, dy, dz)
    enddo

    do pi = 1, f%pnum
      bf => f%p(pi)

      if (bf_cutoff(bf, rr)) cycle

      call p_calc(bf, i, r2, dx, dy, dz)
    enddo

    do di = 1, f%dnum
      bf => f%d(di)

      if (bf_cutoff(bf, rr)) cycle

      call d_calc(bf, i, r2, dx, dy, dz)
    enddo

    do fi = 1, f%fnum
      bf => f%f(fi)

      if (bf_cutoff(bf, rr)) cycle

      call f_calc(bf, i, r2, dx, dy, dz)
    enddo
  end subroutine calcAOMO

  subroutine sortBFsElec(elecBFs, bfs, i1, i2, rai)
    type(elec_bf_t), dimension(:), allocatable, intent(inout) :: elecBFs
    type(bf_t), dimension(:), allocatable, target, intent(inout) :: bfs
    integer, intent(in) :: i1, i2
    real(r8), intent(in) :: rai(:, :)

    type(bf_t), pointer :: bf

    integer :: nf,moc,bfi,i,j,a
    real(r8) :: rr  !, ers
    real(r8) :: moc_max, alp_min, crit

    do i = i1, i2
      call elecBF_init(elecBFs(i), nbasf)
    enddo

    moc = 1
    do bfi = 1, nbasf
      a = bc(bfi)
      bf => bfs(bfi)

      select case(bl(bfi))
      case('S')
        nf = norb
        bf%type = BF_S
      case('P')
        nf = 3*norb
        bf%type = BF_P
      case('D')
        nf = 6*norb
        bf%type = BF_D
      case('F')
        nf = 10*norb
        bf%type = BF_F
      case default
        call abortp('(aomo_task): wrong GTO')
        ! just to silence compiler warnings about uninitialized variables
        return
      end select

      call bf_build(bf, bfi, moc, nf, cutao)
      bf%atom = a
      moc = moc + nf

      ! maxval is surprisingly slow, so instead of using it, we roll our own
      ! loop to get the maximal MO coefficient to use it for the product cutoff
      moc_max = 0
      do j = 1, norb
       if(bf%moc(j) > moc_max) moc_max = abs(bf%moc(j))
      enddo

      alp_min = 1d100
      do j = 1, bf%ngto
        if(bf%alpha(j) < alp_min) alp_min = bf%alpha(j)
      enddo
      ! moc*exp(-a*r^2) < cutoff
      ! <=> -a*r^2 < ln(cutoff/moc)
      ! <=> r^2 > -ln(cutoff/moc) / a
      crit = sqrt(-log(prodcutoff / moc_max) / alp_min)

      elec: do i = i1, i2
        if(cutao) then
          rr = rai(a, i)
          if (rr .gt. bf%aocut) then
            cycle elec
          endif
        endif

        if(prodcutoff > 0) then
          rr = rai(a, i)
          ! the check below is equivalent to
          ! ers = exp(-alp_min * rr*rr)
          ! if(moc_max * ers < prodcutoff) then
          if(rr .gt. crit) then
            cycle elec
          endif
        endif

        elecBFs(i)%cnt = elecBFs(i)%cnt + 1
        elecBFs(i)%bfs(elecBFs(i)%cnt) = bfi
      enddo elec
    enddo
  end subroutine sortBFsElec

  subroutine calcAOMOelec(bfs, elecBF, i, rai, xx, yy, zz)
    type(bf_t), dimension(:), target :: bfs
    type(elec_bf_t), intent(in) :: elecBF
    integer, intent(in) :: i
    real(r8), intent(in) :: rai(:, :)
    real(r8), intent(in) :: xx, yy, zz

    type(bf_t), pointer :: bf
    real(r8) :: rr, r2, dx, dy, dz
    integer :: a, f

    do f = 1, elecBF%cnt
      bf => bfs(elecBF%bfs(f))

      a = bf%atom
      rr = rai(a, i)
      r2 = rr*rr
      dx = xx-atoms(a)%cx
      dy = yy-atoms(a)%cy
      dz = zz-atoms(a)%cz

      select case(bf%type)
      case(BF_S)
        call s_calc(bf, i, r2, dx, dy, dz)
      case(BF_P)
        call p_calc(bf, i, r2, dx, dy, dz)
      case(BF_D)
        call d_calc(bf, i, r2, dx, dy, dz)
      case(BF_F)
        call f_calc(bf, i, r2, dx, dy, dz)
      end select
    enddo
  end subroutine calcAOMOelec

  subroutine s_calc(bf, i, r2, dx, dy, dz)
    type(bf_t), intent(inout) :: bf
    real(r8), intent(in) :: r2, dx, dy, dz
    integer, intent(in) :: i

    integer :: j
    real(r8) :: alp, u, ux
    real(r8) :: tmp, tmps(5)

    tmps = 0d0

    ! loop over contraction
    do j = 1, bf%ngto
      alp = bf%alpha(j)
      u = bf%coeff(j) * exp(-alp*r2)
      ux = -2d0*alp*u

      tmps(1) = tmps(1) + u
      tmps(2) = tmps(2) + ux*dx
      tmps(3) = tmps(3) + ux*dy
      tmps(4) = tmps(4) + ux*dz
      tmps(5) = tmps(5) + ux*(3d0-2d0*alp*r2)
    enddo

    ! MO calculation
    do j = 1, norb
      tmp = bf%moc(j)
      mat(j,i,1)   = mat(j,i,1)   + tmp*tmps(1)
      mat1x(j,i,1) = mat1x(j,i,1) + tmp*tmps(2)
      mat1y(j,i,1) = mat1y(j,i,1) + tmp*tmps(3)
      mat1z(j,i,1) = mat1z(j,i,1) + tmp*tmps(4)
      mat2(j,i,1)  = mat2(j,i,1)  + tmp*tmps(5)
    enddo
  end subroutine s_calc

  subroutine p_calc(bf, i, r2, dx, dy, dz)
    type(bf_t), intent(inout) :: bf
    real(r8), intent(in) :: r2, dx, dy, dz
    integer, intent(in) :: i

    integer :: j, d
    real(r8) :: alp,u,ux
    real(r8) :: dx2, dy2, dz2, dxdy, dxdz, dydz
    real(r8) :: tmp,tmpp(3,5)

    ! do all 3 P simultaneously (same exponent is required)
    ! order p_x,p_y,p_z

    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz

    dxdy = dx*dy
    dxdz = dx*dz
    dydz = dy*dz

    tmpp = 0d0

    do j=1, bf%ngto
      alp = bf%alpha(j)
      u = bf%coeff(j) * exp(-alp*r2)
      ux = -2d0*alp*u

      tmpp(1,1) = tmpp(1,1) + dx*u
      tmpp(2,1) = tmpp(2,1) + dy*u
      tmpp(3,1) = tmpp(3,1) + dz*u

      tmpp(1,2) = tmpp(1,2) + u + ux*dx2
      tmpp(2,2) = tmpp(2,2) + ux*dxdy
      tmpp(3,2) = tmpp(3,2) + ux*dxdz
      tmpp(1,3) = tmpp(1,3) + ux*dxdy
      tmpp(2,3) = tmpp(2,3) + u + ux*dy2
      tmpp(3,3) = tmpp(3,3) + ux*dydz
      tmpp(1,4) = tmpp(1,4) + ux*dxdz
      tmpp(2,4) = tmpp(2,4) + ux*dydz
      tmpp(3,4) = tmpp(3,4) + u + ux*dz2

      tmp = (5d0-2d0*alp*r2)*ux
      tmpp(1,5) = tmpp(1,5) + tmp*dx
      tmpp(2,5) = tmpp(2,5) + tmp*dy
      tmpp(3,5) = tmpp(3,5) + tmp*dz
    enddo

    ! MO calculation
    do j=1, norb
      do d=1, 3
        tmp = bf%moc(3*j+d-3)
        mat(j,i,1)   = mat(j,i,1)   + tmp*tmpp(d,1)
        mat1x(j,i,1) = mat1x(j,i,1) + tmp*tmpp(d,2)
        mat1y(j,i,1) = mat1y(j,i,1) + tmp*tmpp(d,3)
        mat1z(j,i,1) = mat1z(j,i,1) + tmp*tmpp(d,4)
        mat2(j,i,1)  = mat2(j,i,1)  + tmp*tmpp(d,5)
      enddo
    enddo
  end subroutine p_calc

  subroutine d_calc(bf, i, r2, dx, dy, dz)
    type(bf_t), intent(inout) :: bf
    real(r8), intent(in) :: r2, dx, dy, dz
    integer, intent(in) :: i

    integer :: ic, j, d
    real(r8) :: alp,u,ux
    real(r8) :: dx2,dy2,dz2,dxdy,dxdz,dydz,dxdydz,dy2dx,dx2dy,dx2dz,dz2dx,dy2dz,dz2dy
    real(r8) :: tmp,tmpd(6,5)

    real(r8), parameter :: sqr3 = 1.73205080756887729d0

    ! do all 6 D simultaneously (same exponent is required)
    ! order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)

    dx2   = dx*dx
    dy2   = dy*dy
    dz2   = dz*dz

    dxdy  = dx*dy
    dxdydz =dxdy*dz
    dy2dx = dxdy*dy
    dx2dy = dxdy*dx
    dxdz  = dx*dz
    dx2dz = dxdz*dx
    dz2dx = dxdz*dz
    dydz  = dy*dz
    dy2dz = dydz*dy
    dz2dy = dydz*dz

    tmpd = 0d0

    do ic=1, bf%ngto
      alp = bf%alpha(ic)
      u = bf%coeff(ic) * exp(-alp*r2)
      ux = -2d0*alp*u

      tmpd(1,1) = tmpd(1,1) + dx2*u
      tmpd(2,1) = tmpd(2,1) + dy2*u
      tmpd(3,1) = tmpd(3,1) + dz2*u

      tmpd(1,2) = tmpd(1,2) + (2d0*u + ux*dx2)*dx
      tmpd(2,2) = tmpd(2,2) + dy2dx*ux
      tmpd(3,2) = tmpd(3,2) + dz2dx*ux
      tmpd(1,3) = tmpd(1,3) + dx2dy*ux
      tmpd(2,3) = tmpd(2,3) + (2d0*u + ux*dy2)*dy
      tmpd(3,3) = tmpd(3,3) + dz2dy*ux
      tmpd(1,4) = tmpd(1,4) + dx2dz*ux
      tmpd(2,4) = tmpd(2,4) + dy2dz*ux
      tmpd(3,4) = tmpd(3,4) + (2d0*u + ux*dz2)*dz
      tmp       = (7d0 - 2d0*alp*r2)*ux
      tmpd(1,5) = tmpd(1,5) + 2d0*u + dx2*tmp
      tmpd(2,5) = tmpd(2,5) + 2d0*u + dy2*tmp
      tmpd(3,5) = tmpd(3,5) + 2d0*u + dz2*tmp

      ! correction of norm, N(dxx)*sqr3 = N(dxy)
      u = sqr3*u
      ux = sqr3*ux

      tmpd(4,1) = tmpd(4,1) + dxdy*u
      tmpd(5,1) = tmpd(5,1) + dxdz*u
      tmpd(6,1) = tmpd(6,1) + dydz*u
      tmp = ux*dxdydz
      tmpd(4,2) = tmpd(4,2) + (u + ux*dx2)*dy
      tmpd(5,2) = tmpd(5,2) + (u + ux*dx2)*dz
      tmpd(6,2) = tmpd(6,2) + tmp
      tmpd(4,3) = tmpd(4,3) + (u + ux*dy2)*dx
      tmpd(5,3) = tmpd(5,3) + tmp
      tmpd(6,3) = tmpd(6,3) + (u + ux*dy2)*dz
      tmpd(4,4) = tmpd(4,4) + tmp
      tmpd(5,4) = tmpd(5,4) + (u + ux*dz2)*dx
      tmpd(6,4) = tmpd(6,4) + (u + ux*dz2)*dy
      tmp = (7d0 - 2d0*alp*r2)*ux
      tmpd(4,5) = tmpd(4,5) + tmp*dxdy
      tmpd(5,5) = tmpd(5,5) + tmp*dxdz
      tmpd(6,5) = tmpd(6,5) + tmp*dydz
    enddo

    ! MO calculation
    do j=1, norb
      do d=1, 6
        tmp = bf%moc(6*j+d-6)
        mat(j,i,1)   = mat(j,i,1)   + tmp*tmpd(d,1)
        mat1x(j,i,1) = mat1x(j,i,1) + tmp*tmpd(d,2)
        mat1y(j,i,1) = mat1y(j,i,1) + tmp*tmpd(d,3)
        mat1z(j,i,1) = mat1z(j,i,1) + tmp*tmpd(d,4)
        mat2(j,i,1)  = mat2(j,i,1)  + tmp*tmpd(d,5)
      enddo
    enddo
  end subroutine d_calc

  subroutine f_calc(bf, i, r2, dx, dy, dz)
    type(bf_t), intent(inout) :: bf
    real(r8), intent(in) :: r2, dx, dy, dz
    integer, intent(in) :: i

    integer :: ic, j, d
    real(r8) :: alp,u,ux
    real(r8) :: dx2,dy2,dz2,dxyz
    real(r8) :: tmp,tmpf(10,5)

    real(r8), parameter :: sqr3 = 1.73205080756887729d0
    real(r8), parameter :: sqr5 = 2.236067977499789696d0

    ! do all 10 F simultaneously (same exponent is required)

    tmpf = 0d0

    do ic=1, bf%ngto
      alp = bf%alpha(ic)
      u = bf%coeff(ic) * exp(-alp*r2)
      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz
      dxyz = dx*dy*dz
      ux = -2d0*alp*u

      ! f_xxx, f_yyy, f_zzz
      tmpf(1,1) = tmpf(1,1) + dx2*dx*u
      tmpf(2,1) = tmpf(2,1) + dy2*dy*u
      tmpf(3,1) = tmpf(3,1) + dz2*dz*u

      tmpf(1,2) = tmpf(1,2) + (3d0*u + ux*dx2)*dx2
      tmpf(2,2) = tmpf(2,2) + dy2*dy*ux*dx
      tmpf(3,2) = tmpf(3,2) + dz2*dz*ux*dx
      tmpf(1,3) = tmpf(1,3) + dx2*dx*ux*dy
      tmpf(2,3) = tmpf(2,3) + (3d0*u + ux*dy2)*dy2
      tmpf(3,3) = tmpf(3,3) + dz2*dz*ux*dy
      tmpf(1,4) = tmpf(1,4) + dx2*dx*ux*dz
      tmpf(2,4) = tmpf(2,4) + dy2*dy*ux*dz
      tmpf(3,4) = tmpf(3,4) + (3d0*u + ux*dz2)*dz2
      tmp = (9d0 - 2d0*alp*r2)*ux
      tmpf(1,5) = tmpf(1,5) + (6d0*u + dx2*tmp)*dx
      tmpf(2,5) = tmpf(2,5) + (6d0*u + dy2*tmp)*dy
      tmpf(3,5) = tmpf(3,5) + (6d0*u + dz2*tmp)*dz

      ! f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
      ! correction of norm, N(fxxx)*sqrt(5) = N(fxxy)
      u = sqr5*u
      ux = sqr5*ux

      tmpf(5,1) = tmpf(5,1) + dx2*dy*u
      tmpf(6,1) = tmpf(6,1) + dx2*dz*u
      tmpf(4,1) = tmpf(4,1) + dy2*dx*u
      tmpf(9,1) = tmpf(9,1) + dy2*dz*u
      tmpf(7,1) = tmpf(7,1) + dz2*dx*u
      tmpf(8,1) = tmpf(8,1) + dz2*dy*u

      ! derivatives
      tmp = ux*dxyz
      tmpf(5,2) = tmpf(5,2) + (2d0*u + ux*dx2)*dx*dy
      tmpf(6,2) = tmpf(6,2) + (2d0*u + ux*dx2)*dx*dz
      tmpf(4,2) = tmpf(4,2) + (u + ux*dx2)*dy2
      tmpf(9,2) = tmpf(9,2) + tmp*dy
      tmpf(7,2) = tmpf(7,2) + (u + ux*dx2)*dz2
      tmpf(8,2) = tmpf(8,2) + tmp*dz
      tmpf(5,3) = tmpf(5,3) + (u + ux*dy2)*dx2
      tmpf(6,3) = tmpf(6,3) + tmp*dx
      tmpf(4,3) = tmpf(4,3) + (2d0*u + ux*dy2)*dx*dy
      tmpf(9,3) = tmpf(9,3) + (2d0*u + ux*dy2)*dy*dz
      tmpf(7,3) = tmpf(7,3) + tmp*dz
      tmpf(8,3) = tmpf(8,3) + (u + ux*dy2)*dz2
      tmpf(5,4) = tmpf(5,4) + tmp*dx
      tmpf(6,4) = tmpf(6,4) + (u + ux*dz2)*dx2
      tmpf(4,4) = tmpf(4,4) + tmp*dy
      tmpf(9,4) = tmpf(9,4) + (u + ux*dz2)*dy2
      tmpf(7,4) = tmpf(7,4) + (2d0*u + ux*dz2)*dx*dz
      tmpf(8,4) = tmpf(8,4) + (2d0*u + ux*dz2)*dy*dz
      ! laplacians
      tmp = (9d0 - 2d0*alp*r2)*ux
      tmpf(5,5) = tmpf(5,5) + (2d0*u + dx2*tmp)*dy
      tmpf(6,5) = tmpf(6,5) + (2d0*u + dx2*tmp)*dz
      tmpf(4,5) = tmpf(4,5) + (2d0*u + dy2*tmp)*dx
      tmpf(9,5) = tmpf(9,5) + (2d0*u + dy2*tmp)*dz
      tmpf(7,5) = tmpf(7,5) + (2d0*u + dz2*tmp)*dx
      tmpf(8,5) = tmpf(8,5) + (2d0*u + dz2*tmp)*dy

      ! f_xyz
      ! correction of norm, N(fxxx)*sqrt(15) = N(fxxy)*sqrt(3) = N(fxyz)
      u  = sqr3*u
      ux = sqr3*ux

      tmpf(10,1) = tmpf(10,1) + dxyz*u

      tmpf(10,2) = tmpf(10,2) + (u + ux*dx2)*dy*dz
      tmpf(10,3) = tmpf(10,3) + (u + ux*dy2)*dx*dz
      tmpf(10,4) = tmpf(10,4) + (u + ux*dz2)*dx*dy
      tmp = (9d0 - 2d0*alp*r2)*ux
      tmpf(10,5) = tmpf(10,5) + dxyz*tmp

    enddo

    ! MO calculation
    do j=1, norb
      do d=1, 10
        tmp = bf%moc(10*j+d-10)
        mat(j,i,1)   = mat(j,i,1)   + tmp*tmpf(d,1)
        mat1x(j,i,1) = mat1x(j,i,1) + tmp*tmpf(d,2)
        mat1y(j,i,1) = mat1y(j,i,1) + tmp*tmpf(d,3)
        mat1z(j,i,1) = mat1z(j,i,1) + tmp*tmpf(d,4)
        mat2(j,i,1)  = mat2(j,i,1)  + tmp*tmpf(d,5)
      enddo
    enddo
  end subroutine f_calc
end module aoMoTask_m
