! Copyright (C) 2011-2013 Alexander Sturm
! Copyright (C) 2012-2013, 2015 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module jastrowDTN_m

  use kinds_m, only: r8
  use error_m
  use wfData_m
  use jastrowParamData_m
  use linkedList_m
  implicit none

  private
  public :: jasinput_dtn, jasoutput_dtn, jasoutput_dtn_new, jas_shortoutput_dtn, &
            jasChangeType_dtn, jasdtnall, jasdtnp, &
            getVectorLenDTN, getVectorDTN, putVectorDTN, jasdtnall_naiv, &
            jasgeneric_naiv, getNumConstraintsDTN, getConstraintsDTN

  ! max index of parameter
  integer :: umax, xmax, fenmax, feemax
  ! number of parameters per term
  integer :: unum, xnum, fnum
  ! cutoff behavior
  integer :: C
  ! cutoff lengths for e-e, e-n
  real(r8) :: cutOffEE, cutOffEN, cutOffEEN
  ! coefficients for power expansions
  logical :: allocatedParams = .false.
  real(r8), allocatable :: alpha(:, :) ! 2, kmax
  real(r8), allocatable :: beta(:, :, :) ! 2, amax, kmax
  real(r8), allocatable :: gamma(:,:,:,:,:) ! 2, amax, 0:kmax, 0:kmax, 0:kmax
  ! last nucleus with correlation
  integer :: nclast
  ! whether electron nucleus should be enforced
  logical :: nucCusp
  ! XXX current optmode
  integer :: curOptMode = 1
  ! spin dependencies in u, x and f
  integer :: spinU = 1, spinX = 1, spinF = 1

contains

!=======================================================================

subroutine allocateParams()
  integer :: istat
  if(allocatedParams) return

  allocate(alpha(2, umax), beta(2, nclast, xmax), &
                gamma(2, nclast, 0:fenmax, 0:fenmax, 0:feemax), stat=istat)
  call assert(istat == 0, "jastrowDTN_m: parameter array allocation failed")
  alpha = 0
  beta = 0
  gamma = 0
  allocatedParams = .true.
end subroutine allocateParams

subroutine deallocateParams()
  if(.not. allocatedParams) return
  deallocate(alpha, beta, gamma)
  allocatedParams = .false.
end subroutine deallocateParams

! read Jastrow related input from 'lines'
subroutine jasinput_dtn(lines,nl)
  character(len=*), intent(in) :: lines(:)! lines array
  integer, intent(in)          :: nl      ! actual # of lines
  integer :: i, j, k, m, s, idx
  character :: t

  read(lines(2), *) nclast
  if(nclast == 0) nclast = atoms_getNCLast(atoms)

  read(lines(3), *) C
  read(lines(4), *) nucCusp
  read(lines(5), *) spinU, spinX, spinF
  spinU = spinU + 1
  spinX = spinX + 1
  spinF = spinF + 1

  read(lines(6), *) umax, xmax, fenmax, feemax

  ! deallocate the params arrays first, in case they are already allocated
  ! this can happen e.g. during changejastrow
  call deallocateParams()
  call allocateParams()
  ! fenmax and feemax are 0-based
  fenmax = fenmax - 1
  feemax = feemax - 1

  idx = 7
  read(lines(idx), "(A)") t
  idx = idx + 1
  if(t /= 'x') then
    idx = idx - 1
    if(umax > 0 .and. kmax >= umax) then
      do s = 1, spinU
        read(lines(idx+s-1), *) (alpha(s, m), m = 1, umax)
      enddo
      idx = idx + spinU
    elseif(umax < 0) then
      call abortp("(inputwf): umax out of range")
    endif

    if(xmax > 0 .and. kmax >= xmax) then
      do s = 1, spinX
        do m = 1, xmax
          read(lines(idx+m-1), *) (beta(s, i, m), i = 1, nclast)
        enddo
        idx = idx + xmax
      enddo
    elseif(xmax < 0) then
      call abortp("(inputwf): xmax out of range")
    endif

    if(fenmax >= -1 .and. kmax >= fenmax) then
      if(feemax >= -1 .and. kmax >= feemax) then
        do s = 1, spinF
          do i = 0, fenmax
            do j = 0, fenmax
              do k = 0, feemax
                read(lines(idx+k), *) (gamma(s, m, i, j, k), m = 1, nclast)
              enddo
              idx = idx + feemax + 1
            enddo
          enddo
        enddo
      else
        call abortp("(inputwf): feemax out of range")
      endif
    else
      call abortp("(inputwf): fenmax out of range")
    endif
  else
    alpha = 0
    beta = 0
    gamma = 0
  endif

  read(lines(idx), *) cutOffEE, cutOffEN, cutOffEEN

  fnum = nclast * (fenmax + 1) * (fenmax + 1) * (feemax + 1)
  xnum = nclast * xmax
  unum = umax

  !call checkCusp()
end subroutine jasinput_dtn

!=======================================================================

! write Jastrow related terms to log file unit 'iul'
subroutine jasoutput_dtn
  integer :: a, i, j, k, s
  if(mytid /= 0) return

  write(iul, "(/A)") "Jastrow part:"
  if(umax > 0) then
    write(iul, "(/A,I4)") "No. of el-el Jastrow terms = ", umax * spinU
    do s = 1, spinU
      write(iul, "(A,I4)") "Spin function: ", s
      do i = 1, umax
        write(iul, "(F22.15)", advance="no") alpha(s, i)
      enddo
      write(iul, *)
    enddo
  endif

  write(iul, *)

  if(xmax > 0) then
    write(iul, "(/A,I4)") "No. of el-nuc Jastrow terms = ", xmax * spinX
    do s = 1, spinX
      write(iul, "(A,I4)") "Spin function: ", s
      do i = 1, xmax
        do a = 1, nclast
          write(iul, "(F22.15)", advance="no") beta(s, a, i)
        enddo
        write(iul, *)
      enddo
    enddo
  endif

  write(iul, *)

  if(fenmax >= 0 .and. feemax >= 0) then
    write(iul, "(/A,I4)") "No. of el-nuc Jastrow f-terms = ", (fenmax+1) * spinF
    write(iul, "(A,I4)") "No. of el-el Jastrow f-terms = ", (feemax+1) * spinF
    do s = 1, spinF
      write(iul, "(A,I4)") "Spin function: ", s
      do i = 0, fenmax
        do j = 0, fenmax
          do k = 0, feemax
            do a = 1, nclast
              write(iul, "(F22.15)", advance="no") gamma(s, a, i, j, k)
            enddo
            write(iul, *)
          enddo
        enddo
      enddo
    enddo
  endif

  write(iul, *)

end subroutine jasoutput_dtn

subroutine jasoutput_dtn_new(iu)
  integer, intent(in) :: iu
  integer :: i, j, k, m, s

  write(iu, "(I4)") nclast
  write(iu, "(I4)") C
  write(iu, "(L1)") nucCusp
  write(iu, "(3I4)") spinU - 1, spinX - 1, spinF -1

  write(iu, "(4I4)") umax, xmax, fenmax + 1, feemax + 1
  if(umax > 0) then
    do s = 1, spinU
      do m = 1, umax
        write(iu, "(F12.6)", advance="no") alpha(s, m)
      enddo
      write(iu, *)
    enddo
  endif

  if(xmax > 0) then
    do s = 1, spinX
      do m = 1, xmax
        do i = 1, nclast
          write(iu, "(F12.6)", advance="no") beta(s, i, m)
        enddo
        write(iu, *)
      enddo
    enddo
  endif

  if(fenmax >= -1 .and. feemax >= -1) then
    do s = 1, spinF
      do i = 0, fenmax
        do j = 0, fenmax
          do k = 0, feemax
            do m = 1, nclast
              write(iu, "(F12.6)", advance="no") gamma(s, m, i, j, k)
            enddo
            write(iu, *)
          enddo
        enddo
      enddo
    enddo
  endif

  write(iu, "(3F12.6)") cutOffEE, cutOffEN, cutOffEEN
end subroutine jasoutput_dtn_new


subroutine jas_shortoutput_dtn(iu)
  integer, intent(in) :: iu
  integer :: i, j, k, m, s

  write(iu, "(4(a,I4))") ' umax=',umax, ' xmax=',xmax,' fenmax=',fenmax,' feemax=',feemax
end subroutine jas_shortoutput_dtn

!=======================================================================

! reset Jastrow parameters (to start e.g. the optimization)
subroutine jasChangeType_dtn()
  alpha = 0d0
  beta = 0d0
  gamma = 0d0
  cutOffEE = 1d0
  cutOffEN = 1d0
  cutOffEEN = 1d0
  curOptMode = 1
  nucCusp = .false.
  spinU = 1
  spinX = 1
  spinF = 1
  nclast = atoms_getNCLast(atoms)
end subroutine jasChangeType_dtn

!=======================================================================


subroutine jasdtnall(x, y, z, rai, rij, optType, ju, jud, julapl, julapli)

  real(r8), intent(in)  ::  x(:), y(:), z(:)     ! cartesian coords
  real(r8), intent(in)  ::  rai(:, :), rij(:, :) ! distances r_ai, r_ij
  character(len=*)    ::  optType              ! parameter optimization
  real(r8), intent(out) ::  ju                   ! U
  real(r8), intent(out) ::  jud(:)               ! nabla U
  real(r8), intent(out) ::  julapl               ! sum of laplacians
  real(r8), intent(out) ::  julapli(:)           ! laplacian U

  ! values of correlation terms
  real(r8) :: uTerm, xTerm, fTerm
  ! nabla of correlation terms
  real(r8) :: uDeriv(3*ne), xDeriv(3*ne), fDeriv(3*ne)
  ! laplacian
  real(r8) :: uLapli(ne), xLapli(ne), fLapli(ne)

  ! parameter derivatives
  logical :: deriveParams

  ! grids for electrons
  type(List), target, allocatable, dimension(:,:,:) :: electronGrid, nucleusGrid, eenGrid
  type(List), pointer :: cell, cellb
  type(ListNode), pointer :: node, nodeb
  integer :: coordx, coordy, coordz
  real(r8) :: minX, minY, minZ, maxX, maxY, maxZ, diffX, diffY, diffZ

  ! spin coefficient for e-e-term
  real(r8) :: spinCoeff(ne, ne)
  ! numerator in e-n term, used to enforce cusp if necessary
  real(r8) :: nucNum(nclast)

  ! x/y/z distances e-n and e-e
  real(r8) :: xai(nclast, ne), yai(nclast, ne), zai(nclast, ne)
  real(r8) :: xij(ne, ne), yij(ne, ne), zij(ne, ne)

  ! maximum power of cutoff distances needed

  ! powers of distances
  real(r8) :: eePowers(ne, ne, -2:max(umax, feemax, 1))
  real(r8) :: enPowers(nclast, ne, -2:max(xmax, fenmax, 1))

  ! derivatives of distances - 1st index: type (ai, aj, ij)
  real(r8) :: rDeriv(3, 3), rLapl(3), rComb(3) !XXX maybe reuse in e-e-n?

  ! (distance - cut-off length) and powers
  real(r8) :: dist, distj, distP(-1:C), distjP(-1:C)
  ! sum in x, u and f-term, derivative, nabla
  real(r8) :: sumr, sumrd, sumrn, sumpr
  ! commonly used expressions in sum/derivative/nabla/laplacian of u/x/f
  real(r8) :: uterms(5), xterms(5), fterms(8)
  real(r8) :: fsum(8), fsuml(2)
  ! for parameter derivatives
  real(r8) :: fsump(8), fsumpl(2)
  ! position of parameter for param derivatives
  integer :: pos, offset, coOffset
  ! array index for parameters of spin dependent functions
  integer :: spin

  ! commonly used terms
  real(r8) :: g, p, q, r, s, u, v, w, tmp, tmpa(3)

  ! loop variables
  integer :: a, i, j, l, m, n, t, ll, mm, nn, xx, yy, zz, np1, np2, npnl

  ! alloc stat
  integer :: error

  !write(iul, *) "start"


  !---------Initialization of arrays-------------------------------
  ju = 0
  jud = 0
  julapl = 0
  julapli = 0

  uTerm = 0
  xTerm = 0
  fTerm = 0
  uDeriv = 0
  xDeriv = 0
  fDeriv = 0
  uLapli = 0
  xLapli = 0
  fLapli = 0

  minX = min(minval(x), minval(atoms(:)%cx))
  minY = min(minval(y), minval(atoms(:)%cy))
  minZ = min(minval(z), minval(atoms(:)%cz))
  maxX = max(maxval(x), maxval(atoms(:)%cx))
  maxY = max(maxval(y), maxval(atoms(:)%cy))
  maxZ = max(maxval(z), maxval(atoms(:)%cz))
  !write(iul, *) "max: ", maxX, maxY, maxZ
  !write(iul, *) "min: ", minX, minY, minZ

  diffX = maxX - minX
  diffY = maxY - minY
  diffZ = maxZ - minZ

  !write(iul, *) "prealloc"
  ! increase size by 1 in each direction so we can check all surrounding cells
  ! for every electron - even if it's in the last cell
  call assert(cutOffEE > 0, "electron-electron cutoff <= 0")
  call assert(cutOffEN > 0, "electron-nucleus cutoff <= 0")
  call assert(cutOffEEN > 0, "electron-electron-nucleus cutoff <= 0")

  allocate(electronGrid(-1:int(diffX/cutOffEE)+1, &
                        -1:int(diffY/cutOffEE)+1, &
                        -1:int(diffZ/cutOffEE)+1), stat=error)
  call assert(error == 0, "jasdtnall failed to allocate electron grid")

  allocate(nucleusGrid(-1:int(diffX/cutOffEN)+1, &
                       -1:int(diffY/cutOffEN)+1, &
                       -1:int(diffZ/cutOffEN)+1), stat=error)
  call assert(error == 0, "jasdtnall failed to allocate electron-nucleus grid")

  allocate(eenGrid(-1:int(diffX/cutOffEEN)+1, &
                   -1:int(diffY/cutOffEEN)+1, &
                   -1:int(diffZ/cutOffEEN)+1), stat=error)
  call assert(error == 0, "jasdtnall failed to allocate e-e-n grid")
  !write(iul, *) "allocated"

  deriveParams = .false.
  offset = 0
  if(optType == "jastrow" .or. optType == "jas+ci" .or. optType == "jas+mo" &
                                           .or. optType == "jas+mo+ci"  ) then
    deriveParams = .true.
    call assert(allocated(uk), "jastrowParamData not allocated")
    call getVectorLenDTN(curOptMode,np1,np2,npnl)
    call assert(uk_params == np1+np2+npnl,"jastdtnall: mismatch of uk_params")
    if(curOptMode == 1) offset = getOffset()
    uk = 0
    ukgrad = 0
    uklapl = 0
    uklapli = 0
  endif

  !---------Precalculation of Terms---------------------------------------

  if(nucCusp) then
    nucNum = atoms(1:nclast)%za
  else
    nucNum = 0
  endif

  rDeriv = 0
  distP(-1) = 0
  distP(0) = 1
  distjP(-1) = 0
  distjP(0) = 1
  do i = 1, ne
    do a = 1, nclast
      ! x/y/z distances for e-n
      xai(a, i) = x(i) - atoms(a)%cx
      yai(a, i) = y(i) - atoms(a)%cy
      zai(a, i) = z(i) - atoms(a)%cz

      ! powers of e-n distances
      enPowers(a, i, -2) = 0d0
      enPowers(a, i, -1) = 0d0
      enPowers(a, i, 0) = 1d0
      enPowers(a, i, 1) = rai(a, i)
      do t = 2, max(xmax, fenmax)
        enPowers(a, i, t) = enPowers(a, i, t-1) * rai(a, i)
      enddo
    enddo
    do j = i + 1, ne
      ! x/y/z distance for e-e
      xij(i, j) = x(i) - x(j)
      yij(i, j) = y(i) - y(j)
      zij(i, j) = z(i) - z(j)

      ! powers of e-e distances
      eePowers(i, j, -2) = 0d0
      eePowers(i, j, -1) = 0d0
      eePowers(i, j, 0) = 1d0
      eePowers(i, j, 1) = rij(i, j)
      do t = 2, max(umax, feemax)
        eePowers(i, j, t) = eePowers(i, j, t-1) * rij(i, j)
      enddo
    enddo
  enddo

  ! spin coefficients for e-e correlation
  spinCoeff(1:nalpha, 1:nalpha) = 0.25d0
  spinCoeff(1:nalpha, 1 + nalpha:ne) = 0.5d0
  spinCoeff(1 + nalpha:ne, 1 + nalpha:ne) = 0.25d0


  ! To simplify cut-off length check sort electrons/nuclei into grid
  ! where each cell has x/y/z-length = cut-off length
  do i = 1, ne
    ! sort electrons into grid for electron-electron cutoff
    coordx = int((x(i) - minX)/cutOffEE)
    coordy = int((y(i) - minY)/cutOffEE)
    coordz = int((z(i) - minZ)/cutOffEE)
    cell => electronGrid(coordx, coordy, coordz)
    call insertNode(cell, i)

    ! e-n cutoff
    coordx = int((x(i) - minX)/cutOffEN)
    coordy = int((y(i) - minY)/cutOffEN)
    coordz = int((z(i) - minZ)/cutOffEN)
    cell => nucleusGrid(coordx, coordy, coordz)
    call insertNode(cell, i)

    ! e-e-n cutoff
    coordx = int((x(i) - minX)/cutOffEEN)
    coordy = int((y(i) - minY)/cutOffEEN)
    coordz = int((z(i) - minZ)/cutOffEEN)
    cell => eenGrid(coordx, coordy, coordz)
    call insertNode(cell, i)
  enddo

  !write(iul, *) "postinit"

  !---------Electron-Electron-Correlation Terms---------------------------

  r = (-cutOffEE)**C
  v = C / cutOffEE

  do i = 1, ne - 1
    coordx = int((x(i) - minX)/cutOffEE)
    coordy = int((y(i) - minY)/cutOffEE)
    coordz = int((z(i) - minZ)/cutOffEE)

    ! check adjacent grid cells
    do l = -1, 1
      do m = -1, 1
        do n = -1, 1

          cell => electronGrid(coordx + l, coordy + m, coordz + n)
          node => cell%node

          ! iterate over electrons in current cell
          do while(associated(node))
            j = node%idx

            dist = rij(i, j) - cutOffEE ! XXX actually the negative distance
            ! only calculate correlation if electron-electron distance is <= cutoff
            if(i < j .and. dist < 0) then
              if(spinU == 1 .or. ((i <= nalpha) .eqv. (j <= nalpha))) then
                spin = 1
              else
                spin = 2
              endif

              s = alpha(spin, 1) * v

              sumrd = spinCoeff(i, j)/r + s
              sumr = alpha(spin, 1) + sumrd * rij(i, j)
              sumrn = 0

              ! calculation of power expansion sum/derivatives
              do t = 2, umax
                sumr = sumr + alpha(spin, t) * eePowers(i, j, t)
                sumrd = sumrd + t * alpha(spin, t) * eePowers(i, j, t-1)
                sumrn = sumrn + t * (t-1) * alpha(spin, t) * eePowers(i, j, t-2)
              enddo

              ! derivatives dr/d{x,y,z}
              rDeriv(3, 1) = xij(i, j) / rij(i, j)
              rDeriv(3, 2) = yij(i, j) / rij(i, j)
              rDeriv(3, 3) = zij(i, j) / rij(i, j)
              ! laplacian of r
              rLapl(3) = 2 / rij(i, j)

              distP(1) = dist
              do t = 2, C
                distP(t) = dist * distP(t-1)
              enddo

              uterms(1) = distP(C)
              uterms(2) = C * distP(C-1)
              uterms(3) = C * (C-1) * distP(C-2)
              uterms(4) = 2 * C * distP(C-1)
              uterms(5) = C * (C-1) * (C-2) * distP(C-3)
              u = uterms(1) * sumrd + uterms(2) * sumr

              uTerm = uTerm + uterms(1) * sumr

              uDeriv(3*i-2:3*i) = uDeriv(3*i-2:3*i) + u * rDeriv(3, :)
              uDeriv(3*j-2:3*j) = uDeriv(3*j-2:3*j) - u * rDeriv(3, :)

              tmp = u * rLapl(3) + &
                    uterms(3) * sumr + &
                    uterms(4) * sumrd + &
                    uterms(1) * sumrn
              uLapli(i) = uLapli(i) + tmp
              uLapli(j) = uLapli(j) + tmp

              if(deriveParams .and. unum > 0) then
                if(curOptMode == 1) then
                  ! cutoff derivatives
                  p = - (C * spinCoeff(i, j)/r + s)/cutOffEE
                  q = p * rij(i, j)
                  u = uterms(1) * p + &
                      uterms(2) * q - &
                      uterms(2) * sumrd - &
                      uterms(3) * sumr
                  uk(1) = uk(1) - uterms(2) * sumr + uterms(1) * q

                  tmpa = u * rDeriv(3, :)
                  ukgrad(3*i-2:3*i, 1) = ukgrad(3*i-2:3*i, 1) + tmpa
                  ukgrad(3*j-2:3*j, 1) = ukgrad(3*j-2:3*j, 1) - tmpa

                  tmp = u * rLapl(3) + &
                        uterms(3) * q + &
                        uterms(4) * p - &
                        uterms(5) * sumr - &
                        uterms(2) * sumrn - &
                        2 * uterms(3) * sumrd

                  uklapli(i, 1) = uklapli(i, 1) + tmp
                  uklapli(j, 1) = uklapli(j, 1) + tmp
                  uklapl(1) = uklapl(1) + 2 * tmp
                endif
                ! alpha(1) derivative
                pos = 1 + offset + unum * (spin - 1)

                w = 1 + v * rij(i, j)
                u = uterms(1) * v + uterms(2) * w
                uk(pos) = uk(pos) + uterms(1) * w

                tmpa = u * rDeriv(3, :)
                ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + tmpa
                ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) - tmpa

                tmp = u * rLapl(3) + &
                      uterms(3) * w + &
                      uterms(4) * v
                uklapli(i, pos) = uklapli(i, pos) + tmp
                uklapli(j, pos) = uklapli(j, pos) + tmp
                uklapl(pos) = uklapl(pos) + 2 * tmp

                do t = 2, umax
                  pos = offset + t + umax * (spin - 1)
                  uk(pos) = uk(pos) + uterms(1) * eePowers(i, j, t)

                  tmpa(1:3) = (uterms(2) * eePowers(i, j, t) + &
                               uterms(1) * t * eePowers(i, j, t-1)) * rDeriv(3, :)
                  ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + tmpa
                  ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) - tmpa


                  tmp = (uterms(1) * t * eePowers(i, j, t-1) + &
                         uterms(2) * eePowers(i, j, t)) * rLapl(3) + &
                        uterms(3) * eePowers(i, j, t) + &
                        uterms(4) * t * eePowers(i, j, t-1) + &
                        uterms(1) * t * (t-1) * eePowers(i, j, t-2)
                  uklapli(i, pos) = uklapli(i, pos) + tmp
                  uklapli(j, pos) = uklapli(j, pos) + tmp
                  uklapl(pos) = uklapl(pos) + 2 * tmp

                enddo
              endif

            endif

            node => node%next
          enddo
        enddo
      enddo
    enddo
  enddo

  !write(iul, *) "e-e corr complete"
  !---------Electron-Nucleus-Correlation-Terms----------------------------

  r = (-cutOffEN)**C
  v = C / cutOffEN
  do a = 1, nclast
    coordx = int((atoms(a)%cx - minX)/cutOffEN)
    coordy = int((atoms(a)%cy - minY)/cutOffEN)
    coordz = int((atoms(a)%cz - minZ)/cutOffEN)

    do l = -1, 1
      do m = -1, 1
        do n = -1, 1

          cell => nucleusGrid(coordx + l, coordy + m, coordz + n)
          node => cell%node

          do while(associated(node))
            i = node%idx

            dist = rai(a, i) - cutOffEN
            if(dist < 0) then
              if(spinX == 1 .or. i <= nalpha) then
                spin = 1
              else
                spin = 2
              endif

              s = beta(spin, a, 1) * v
              p = -nucNum(a)/r + s
              w = (C * nucNum(a)/r - s)/cutOffEN

              sumr = beta(spin, a, 1) + p * rai(a, i)
              sumrd = p
              sumrn = 0

              do t = 2, xmax
                sumr = sumr + beta(spin, a, t) * enPowers(a, i, t)
                sumrd = sumrd + t * beta(spin, a, t) * enPowers(a, i, t-1)
                sumrn = sumrn + t * (t-1) * beta(spin, a, t) * enPowers(a, i, t-2)
              enddo

              rDeriv(1, 1) = xai(a, i) / rai(a, i)
              rDeriv(1, 2) = yai(a, i) / rai(a, i)
              rDeriv(1, 3) = zai(a, i) / rai(a, i)
              rLapl(1) = 2 / rai(a, i)

              distP(1) = dist
              do t = 2, C
                distP(t) = dist * distP(t-1)
              enddo

              xterms(1) = distP(C)
              xterms(2) = C * distP(C-1)
              xterms(3) = C * (C-1) * distP(C-2)
              xterms(4) = 2 * C * distP(C-1)
              xterms(5) = C * (C-1) * (C-2) * distP(C-3)
              u = xterms(1) * sumrd + xterms(2) * sumr

              xTerm = xTerm + xterms(1) * sumr
              xDeriv(3*i-2:3*i) = xDeriv(3*i-2:3*i) + u * rDeriv(1, :)
              xLapli(i) = xLapli(i) + &
                          u * rLapl(1) + &
                          xterms(3) * sumr + &
                          xterms(4) * sumrd + &
                          xterms(1) * sumrn

              if(deriveParams .and. xnum > 0) then
                if(curOptMode == 1) then
                  pos = 1
                  if(unum > 0) pos = 2
                  ! cutoff derivatives
                  q = w * rai(a, i)
                  u = xterms(1) * w + &
                      xterms(2) * q - &
                      xterms(2) * sumrd - &
                      xterms(3) * sumr

                  uk(pos) = uk(pos) - xterms(2) * sumr + xterms(1) * q

                  ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + u * rDeriv(1, :)

                  tmp = u * rLapl(1) + &
                        xterms(3) * q + &
                        xterms(4) * w - &
                        xterms(5) * sumr - &
                        xterms(2) * sumrn - &
                        2 * xterms(3) * sumrd

                  uklapli(i, pos) = uklapli(i, pos) + tmp
                  uklapl(pos) = uklapl(pos) + tmp
                endif

                ! beta(1) derivatives
                q = 1 + v * rai(a, i)
                u = xterms(1) * v + xterms(2) * q
                pos = offset + umax * spinU + a + xnum * (spin - 1)

                uk(pos) = uk(pos) + xterms(1) * q

                ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + u * rDeriv(1, :)

                tmp = u * rLapl(1) + &
                      xterms(3) * q + &
                      xterms(4) * v
                uklapli(i, pos) = uklapli(i, pos) + tmp
                uklapl(pos) = uklapl(pos) + tmp

                do t = 2, xmax
                  pos = offset + umax * spinU + (t-1) * nclast + a + xnum * (spin - 1)

                  uk(pos) = uk(pos) + xterms(1) * enPowers(a, i, t)

                  ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + &
                                           (xterms(2) * enPowers(a, i, t) + &
                                            xterms(1) * t * enPowers(a, i, t-1)) * rDeriv(1, :)

                  tmp = (xterms(1) * t * enPowers(a, i, t-1) + &
                         xterms(2) * enPowers(a, i, t)) * rLapl(1) + &
                        xterms(3) * enPowers(a, i, t) + &
                        xterms(4) * t * enPowers(a, i, t-1) + &
                        xterms(1) * t * (t-1) * enPowers(a, i, t-2)
                  uklapli(i, pos) = uklapli(i, pos) + tmp
                  uklapl(pos) = uklapl(pos) + tmp
                enddo
              endif
            endif

            node => node%next
          enddo
        enddo
      enddo
    enddo
  enddo

  !write(iul, *) "e-n corr complete"
  !---------Electron-Electron-Nucleus-Correlation-Terms-------------------
  if(fnum > 0) then
    do a = 1, nclast
      coordx = int((atoms(a)%cx - minX)/cutOffEEN)
      coordy = int((atoms(a)%cy - minY)/cutOffEEN)
      coordz = int((atoms(a)%cz - minZ)/cutOffEEN)

      do xx = -1, 1
        do yy = -1, 1
          do zz = -1, 1

            cell => eenGrid(coordx + xx, coordy + yy, coordz + zz)
            node => cell%node

            do while(associated(node))
              i = node%idx
              dist = rai(a, i) - cutOffEEN

              if(dist < 0) then

                do ll = -1, 1
                  do mm = -1, 1
                    do nn = -1, 1

                      cellb => eenGrid(coordx + ll, coordy + mm, coordz + nn)
                      nodeb => cellb%node

                      ! XXX maybe generate only if(associated(nodeb))
                      distP(1) = dist
                      do t = 2, C
                        distP(t) = dist * distP(t-1)
                      enddo

                      do while(associated(nodeb))
                        j = nodeb%idx
                        distj = rai(a, j) - cutOffEEN

                        if(i < j .and. distj < 0) then
                          if(spinF == 1 .or. ((i <= nalpha) .eqv. (j <= nalpha))) then
                            spin = 1
                          else
                            spin = 2
                          endif

                          sumr = 0

                          distjP(1) = distj
                          do t = 2, C
                            distjP(t) = distj * distjP(t-1)
                          enddo

                          fsum = 0

                          do l = 0, fenmax
                            do m = 0, fenmax
                              do n = 0, feemax
                                g = gamma(spin, a, l, m, n)

                                sumr = sumr + g * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n)

                                fsum(1) = fsum(1) + g * l * enPowers(a, i, l-1) * enPowers(a, j, m)   * eePowers(i, j, n)
                                fsum(2) = fsum(2) + g * m * enPowers(a, i, l)   * enPowers(a, j, m-1) * eePowers(i, j, n)
                                fsum(3) = fsum(3) + g * n * enPowers(a, i, l)   * enPowers(a, j, m)   * eePowers(i, j, n-1)
                                fsum(4) = fsum(4) + g * l * (l-1) * enPowers(a, i, l-2) * enPowers(a, j, m)   * eePowers(i, j, n)
                                fsum(5) = fsum(5) + g * m * (m-1) * enPowers(a, i, l)   * enPowers(a, j, m-2) * eePowers(i, j, n)
                                fsum(6) = fsum(6) + g * n * (n-1) * enPowers(a, i, l)   * enPowers(a, j, m)   * eePowers(i, j, n-2)
                                fsum(7) = fsum(7) + g * 2 * l * n * enPowers(a, i, l-1) * enPowers(a, j, m)   * eePowers(i, j, n-1)
                                fsum(8) = fsum(8) + g * 2 * m * n * enPowers(a, i, l) * enPowers(a, j, m-1)   * eePowers(i, j, n-1)
                              enddo
                            enddo
                          enddo

                          rDeriv(1, 1) = xai(a, i) / rai(a, i)
                          rDeriv(1, 2) = yai(a, i) / rai(a, i)
                          rDeriv(1, 3) = zai(a, i) / rai(a, i)
                          rLapl(1) = 2 / rai(a, i)

                          rDeriv(2, 1) = xai(a, j) / rai(a, j)
                          rDeriv(2, 2) = yai(a, j) / rai(a, j)
                          rDeriv(2, 3) = zai(a, j) / rai(a, j)
                          rLapl(2) = 2 / rai(a, j)

                          rDeriv(3, 1) = xij(i, j) / rij(i, j)
                          rDeriv(3, 2) = yij(i, j) / rij(i, j)
                          rDeriv(3, 3) = zij(i, j) / rij(i, j)
                          rLapl(3) = 2 / rij(i, j)

                          rComb(1) = sum(rDeriv(1, :) * rDeriv(3, :))
                          rComb(2) = sum(rDeriv(2, :) * rDeriv(3, :))

                          fterms(1) = distP(C) * distjP(C)
                          fterms(2) = C * distP(C-1) * distjP(C)
                          fterms(3) = C * distP(C) * distjP(C-1)
                          fterms(4) = C * (C-1) * distP(C-2) * distjP(C)
                          fterms(5) = C * (C-1) * distP(C) * distjP(C-2)
                          fterms(6) = 2 * C * distP(C-1) * distjP(C)
                          fterms(7) = 2 * C * distP(C) * distjP(C-1)
                          fterms(8) = C * C * distP(C-1) * distjP(C-1)

                          fTerm = fTerm + fterms(1) * sumr

                          fDeriv(3*i-2:3*i) = fDeriv(3*i-2:3*i) + &
                                              fterms(1) * fsum(1) * rDeriv(1, 1:3) + &
                                              fterms(1) * fsum(3) * rDeriv(3, 1:3) + &
                                              fterms(2) * sumr * rDeriv(1, 1:3)
                          fDeriv(3*j-2:3*j) = fDeriv(3*j-2:3*j) + &
                                              fterms(1) * fsum(2) * rDeriv(2, 1:3) - &
                                              fterms(1) * fsum(3) * rDeriv(3, 1:3) + &
                                              fterms(3) * sumr * rDeriv(2, 1:3)

                          fsuml(1) = fsum(1) * rLapl(1) + &
                                     fsum(3) * rLapl(3) + &
                                     fsum(4) + &
                                     fsum(6) + &
                                     fsum(7) * rComb(1) ! XXX 2*?
                          fsuml(2) = fsum(2) * rLapl(2) + &
                                     fsum(3) * rLapl(3) + &
                                     fsum(5) + &
                                     fsum(6) - &
                                     fsum(8) * rComb(2) ! XXX 2*?

                          fLapli(i) = fLapli(i) + &
                                      fterms(1) * fsuml(1) + &
                                      fterms(2) * sumr * rLapl(1) + &
                                      fterms(4) * sumr + &
                                      fterms(6) * (fsum(1) + fsum(3) * rComb(1))
                          fLapli(j) = fLapli(j) + &
                                      fterms(1) * fsuml(2) + &
                                      fterms(3) * sumr * rLapl(2) + &
                                      fterms(5) * sumr + &
                                      fterms(7) * (fsum(2) - fsum(3) * rComb(2))

                          if(deriveParams) then
                            if(curOptMode == 1) then
                              pos = 1
                              if(unum > 0) pos = pos + 1
                              if(xnum > 0) pos = pos + 1
                              ! cutoff derivatives
                              uk(pos) = uk(pos) - (fterms(2) + fterms(3)) * sumr

                              ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) - &
                                                       (fterms(4) + fterms(8)) * sumr * rDeriv(1, 1:3) - &
                                                       (fterms(2) + fterms(3)) * fsum(1) * rDeriv(1, 1:3) - &
                                                       (fterms(2) + fterms(3)) * fsum(3) * rDeriv(3, 1:3)
                              ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) - &
                                                       (fterms(5) + fterms(8)) * sumr * rDeriv(2, 1:3) - &
                                                       (fterms(2) + fterms(3)) * fsum(2) * rDeriv(2, 1:3) + &
                                                       (fterms(2) + fterms(3)) * fsum(3) * rDeriv(3, 1:3)

                              tmp = (fterms(2) + fterms(3)) * fsuml(1) + &
                                    (fterms(4) + fterms(8)) * sumr * rLapl(1) + &
                                    (C * (C-1) * (C-2) * distP(C-3) * distjP(C) + &
                                     C * C * (C-1) * distP(C-2) * distjP(C-1)) * sumr + &
                                    2 * (fterms(4) + fterms(8)) * (fsum(1) + fsum(3) * rComb(1))
                              uklapli(i, pos) = uklapli(i, pos) - tmp
                              uklapl(pos) = uklapl(pos) - tmp

                              tmp = (fterms(2) + fterms(3)) * fsuml(2) + &
                                    (fterms(5) + fterms(8)) * sumr * rLapl(2) + &
                                    (C * (C-1) * (C-2) * distP(C) * distjP(C-3) + &
                                     C * C * (C-1) * distP(C-1) * distjP(C-2)) * sumr + &
                                    2 * (fterms(5) + fterms(8)) * (fsum(2) - fsum(3) * rComb(2))
                              uklapli(j, pos) = uklapli(j, pos) - tmp
                              uklapl(pos) = uklapl(pos) - tmp
                            endif

                            do l = 0, fenmax
                              do m = 0, fenmax
                                do n = 0, feemax
                                  pos = offset + umax * spinU + xnum * spinX + fnum * (spin - 1) + oneDim(l, m, n, a)

                                  fsump(1) = l * enPowers(a, i, l-1) * enPowers(a, j, m) * eePowers(i, j, n)
                                  fsump(2) = m * enPowers(a, i, l) * enPowers(a, j, m-1) * eePowers(i, j, n)
                                  fsump(3) = n * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n-1)
                                  fsump(4) = l * (l-1) * enPowers(a, i, l-2) * enPowers(a, j, m) * eePowers(i, j, n)
                                  fsump(5) = m * (m-1) * enPowers(a, i, l) * enPowers(a, j, m-2) * eePowers(i, j, n)
                                  fsump(6) = n * (n-1) * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n-2)
                                  fsump(7) = 2 * l * n * enPowers(a, i, l-1) * enPowers(a, j, m) * eePowers(i, j, n-1)
                                  fsump(8) = 2 * m * n * enPowers(a, i, l) * enPowers(a, j, m-1) * eePowers(i, j, n-1)

                                  sumpr = enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n)

                                  uk(pos) = uk(pos) + fterms(1) * sumpr

                                  ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + &
                                                           fterms(1) * fsump(1) * rDeriv(1, 1:3) + &
                                                           fterms(1) * fsump(3) * rDeriv(3, 1:3) + &
                                                           fterms(2) * sumpr * rDeriv(1, 1:3)
                                  ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) + &
                                                           fterms(1) * fsump(2) * rDeriv(2, 1:3) - &
                                                           fterms(1) * fsump(3) * rDeriv(3, 1:3) + &
                                                           fterms(3) * sumpr * rDeriv(2, 1:3)

                                  fsumpl(1) = fsump(1) * rLapl(1) + &
                                              fsump(3) * rLapl(3) + &
                                              fsump(4) + &
                                              fsump(6) + &
                                              fsump(7) * rComb(1)
                                  fsumpl(2) = fsump(2) * rLapl(2) + &
                                              fsump(3) * rLapl(3) + &
                                              fsump(5) + &
                                              fsump(6) - &
                                              fsump(8) * rComb(2)

                                  tmp = fterms(1) * fsumpl(1) + &
                                        fterms(2) * sumpr * rLapl(1) + &
                                        fterms(4) * sumpr + &
                                        fterms(6) * (fsump(1) + fsump(3) * rComb(1))
                                  uklapli(i, pos) = uklapli(i, pos) + tmp
                                  uklapl(pos) = uklapl(pos) + tmp

                                  tmp = fterms(1) * fsumpl(2) + &
                                        fterms(3) * sumpr * rLapl(2) + &
                                        fterms(5) * sumpr + &
                                        fterms(7) * (fsump(2) - fsump(3) * rComb(2))
                                  uklapli(j, pos) = uklapli(j, pos) + tmp
                                  uklapl(pos) = uklapl(pos) + tmp

                                enddo
                              enddo
                            enddo
                          endif
                        endif

                        nodeb => nodeb%next
                      enddo
                    enddo
                  enddo
                enddo
              endif

              node => node%next
            enddo
          enddo
        enddo
      enddo
    enddo
  endif

  !write(iul, *) "e-e-n complete"
  !---------Deallocate/clean up-------------------------------------------

  call deallocateList(electronGrid)
  call deallocateList(nucleusGrid)
  call deallocateList(eenGrid)

  deallocate(electronGrid)
  deallocate(nucleusGrid)
  deallocate(eenGrid)

  !---------Calculation of U and its Deriviatives-------------------------

  ju = uTerm + xTerm + fTerm
  jud = uDeriv + xDeriv + fDeriv
  julapli = uLapli + xLapli + fLapli
  julapl = sum(julapli(1:ne))

  !write(iul, *) uTerm, xTerm, fTerm
  !write(iul, *) ju, julapl

  !---------Parameter derivatives-----------------------------------------

end subroutine jasdtnall

!=========================================================================


subroutine jasdtnall_naiv(x, y, z, rai, rij, optType, ju, jud, julapl, julapli)
  ! main calculations for all electrons moved

  real(r8), intent(in)  ::  x(:), y(:), z(:)     ! cartesian coords
  real(r8), intent(in)  ::  rai(:, :), rij(:, :) ! distances r_ai, r_ij
  character(len=*)    ::  optType              ! parameter optimization
  real(r8), intent(out) ::  ju                   ! U
  real(r8), intent(out) ::  jud(:)               ! nabla U
  real(r8), intent(out) ::  julapl               ! sum of laplacians
  real(r8), intent(out) ::  julapli(:)           ! laplacian U

  ! values of correlation terms
  real(r8) :: uTerm, xTerm, fTerm
  ! nabla of correlation terms
  real(r8) :: uDeriv(3*ne), xDeriv(3*ne), fDeriv(3*ne)
  ! laplacian
  real(r8) :: uLapli(ne), xLapli(ne), fLapli(ne)

  ! parameter derivatives
  logical :: deriveParams

  ! spin coefficient for e-e-term
  real(r8) :: spinCoeff(ne, ne)
  ! numerator in e-n term, used to enforce cusp if necessary
  real(r8) :: nucNum(nclast)

  ! x/y/z distances e-n and e-e
  real(r8) :: xai(nclast, ne), yai(nclast, ne), zai(nclast, ne)
  real(r8) :: xij(ne, ne), yij(ne, ne), zij(ne, ne)

  ! maximum power of cutoff distances needed

  ! powers of distances
  real(r8) :: eePowers(ne, ne, -2:max(umax, feemax, 1))
  real(r8) :: enPowers(nclast, ne, -2:max(xmax, fenmax, 1))

  ! derivatives of distances - 1st index: type (ai, aj, ij)
  real(r8) :: rDeriv(3, 3), rLapl(3), rComb(3) !XXX maybe reuse in e-e-n?

  ! (distance - cut-off length) and powers
  real(r8) :: dist, distj, distP(-1:C), distjP(-1:C)
  ! sum in x, u and f-term, derivative, nabla
  real(r8) :: sumr, sumrd, sumrn, sumpr
  ! commonly used expressions in sum/derivative/nabla/laplacian of u/x/f
  real(r8) :: uterms(5), xterms(5), fterms(8)
  real(r8) :: fsum(8), fsuml(2)
  ! for parameter derivatives
  real(r8) :: fsump(8), fsumpl(2)
  ! position of parameter for param derivatives
  integer :: pos, offset, coOffset
  integer :: spin

  ! commonly used terms
  real(r8) :: g, p, q, r, s, u, v, w, tmp, tmpa(3)

  ! loop variables
  integer :: a, i, j, l, m, n, t, ll, mm, nn, xx, yy, zz, np1, np2, npnl

  ! alloc stat
  integer :: error

  !write(iul, *) "start"

  !---------Initialization of arrays-------------------------------
  ju = 0
  jud = 0
  julapl = 0
  julapli = 0

  uTerm = 0
  xTerm = 0
  fTerm = 0
  uDeriv = 0
  xDeriv = 0
  fDeriv = 0
  uLapli = 0
  xLapli = 0
  fLapli = 0

  offset = 0
  deriveParams = .false.
  if(optType == "jastrow") then
    deriveParams = .true.
    call assert(allocated(uk), "jastrowParamData not allocated")
    call getVectorLenDTN(curOptMode,np1,np2,npnl)
    call assert(uk_params == np1+np2+npnl,"jastdtnall_naiv: mismatch of uk_params")
    if(curOptMode == 1) offset = getOffset()
    uk = 0
    ukgrad = 0
    uklapl = 0
    uklapli = 0
  endif

  !---------Precalculation of Terms---------------------------------------

  if(nucCusp) then
    nucNum = atoms(1:nclast)%za
  else
    nucNum = 0
  endif

  rDeriv = 0
  distP(-1) = 0
  distP(0) = 1
  distjP(-1) = 0
  distjP(0) = 1
  do i = 1, ne
    do a = 1, nclast
      ! x/y/z distances for e-n
      xai(a, i) = x(i) - atoms(a)%cx
      yai(a, i) = y(i) - atoms(a)%cy
      zai(a, i) = z(i) - atoms(a)%cz

      ! powers of e-n distances
      enPowers(a, i, -2) = 0d0
      enPowers(a, i, -1) = 0d0
      enPowers(a, i, 0) = 1d0
      enPowers(a, i, 1) = rai(a, i)
      do t = 2, max(xmax, fenmax)
        enPowers(a, i, t) = enPowers(a, i, t-1) * rai(a, i)
      enddo
    enddo
    do j = i + 1, ne
      ! x/y/z distance for e-e
      xij(i, j) = x(i) - x(j)
      yij(i, j) = y(i) - y(j)
      zij(i, j) = z(i) - z(j)

      ! powers of e-e distances
      eePowers(i, j, -2) = 0d0
      eePowers(i, j, -1) = 0d0
      eePowers(i, j, 0) = 1d0
      eePowers(i, j, 1) = rij(i, j)
      do t = 2, max(umax, feemax)
        eePowers(i, j, t) = eePowers(i, j, t-1) * rij(i, j)
      enddo
    enddo
  enddo

  ! spin coefficients for e-e correlation
  spinCoeff(1:nalpha, 1:nalpha) = 0.25d0
  spinCoeff(1:nalpha, 1 + nalpha:ne) = 0.5d0
  spinCoeff(1 + nalpha:ne, 1 + nalpha:ne) = 0.25d0


  !---------Electron-Electron-Correlation Terms---------------------------

  r = (-cutOffEE)**C
  v = C / cutOffEE

  do i = 1, ne
    do j = i + 1, ne

      dist = rij(i, j) - cutOffEE ! XXX actually the negative distance
      ! only calculate correlation if electron-electron distance is <= cutoff
      if(i < j .and. dist < 0) then
        if(spinU == 1 .or. ((i <= nalpha) .eqv. (j <= nalpha))) then
          spin = 1
        else
          spin = 2
        endif

        s = alpha(spin, 1) * v

        sumrd = spinCoeff(i, j)/r + s
        sumr = alpha(spin, 1) + sumrd * rij(i, j)
        sumrn = 0

        ! calculation of power expansion sum/derivatives
        do t = 2, umax
          sumr = sumr + alpha(spin, t) * eePowers(i, j, t)
          sumrd = sumrd + t * alpha(spin, t) * eePowers(i, j, t-1)
          sumrn = sumrn + t * (t-1) * alpha(spin, t) * eePowers(i, j, t-2)
        enddo

        ! derivatives dr/d{x,y,z}
        rDeriv(3, 1) = xij(i, j) / rij(i, j)
        rDeriv(3, 2) = yij(i, j) / rij(i, j)
        rDeriv(3, 3) = zij(i, j) / rij(i, j)
        ! laplacian of r
        rLapl(3) = 2 / rij(i, j)

        distP(1) = dist
        do t = 2, C
          distP(t) = dist * distP(t-1)
        enddo

        uterms(1) = distP(C)
        uterms(2) = C * distP(C-1)
        uterms(3) = C * (C-1) * distP(C-2)
        uterms(4) = 2 * C * distP(C-1)
        uterms(5) = C * (C-1) * (C-2) * distP(C-3)
        u = uterms(1) * sumrd + uterms(2) * sumr

        uTerm = uTerm + uterms(1) * sumr

        uDeriv(3*i-2:3*i) = uDeriv(3*i-2:3*i) + u * rDeriv(3, :)
        uDeriv(3*j-2:3*j) = uDeriv(3*j-2:3*j) - u * rDeriv(3, :)

        tmp = u * rLapl(3) + &
              uterms(3) * sumr + &
              uterms(4) * sumrd + &
              uterms(1) * sumrn
        uLapli(i) = uLapli(i) + tmp
        uLapli(j) = uLapli(j) + tmp

        if(deriveParams .and. unum > 0) then
          if(curOptMode == 1) then
            ! cutoff derivatives
            p = - (C * spinCoeff(i, j)/r + s)/cutOffEE
            q = p * rij(i, j)
            u = uterms(1) * p + &
                uterms(2) * q - &
                uterms(2) * sumrd - &
                uterms(3) * sumr
            uk(1) = uk(1) - uterms(2) * sumr + uterms(1) * q

            tmpa = u * rDeriv(3, :)
            ukgrad(3*i-2:3*i, 1) = ukgrad(3*i-2:3*i, 1) + tmpa
            ukgrad(3*j-2:3*j, 1) = ukgrad(3*j-2:3*j, 1) - tmpa

            tmp = u * rLapl(3) + &
                  uterms(3) * q + &
                  uterms(4) * p - &
                  uterms(5) * sumr - &
                  uterms(2) * sumrn - &
                  2 * uterms(3) * sumrd

            uklapli(i, 1) = uklapli(i, 1) + tmp
            uklapli(j, 1) = uklapli(j, 1) + tmp
            uklapl(1) = uklapl(1) + 2 * tmp
          endif
          ! alpha(1) derivative
          pos = 1 + offset + unum * (spin - 1)

          w = 1 + v * rij(i, j)
          u = uterms(1) * v + uterms(2) * w
          uk(pos) = uk(pos) + uterms(1) * w

          tmpa = u * rDeriv(3, :)
          ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + tmpa
          ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) - tmpa

          tmp = u * rLapl(3) + &
                uterms(3) * w + &
                uterms(4) * v
          uklapli(i, pos) = uklapli(i, pos) + tmp
          uklapli(j, pos) = uklapli(j, pos) + tmp
          uklapl(pos) = uklapl(pos) + 2 * tmp

          do t = 2, umax
            pos = offset + t + umax * (spin - 1)
            uk(pos) = uk(pos) + uterms(1) * eePowers(i, j, t)

            tmpa(1:3) = (uterms(2) * eePowers(i, j, t) + &
                         uterms(1) * t * eePowers(i, j, t-1)) * rDeriv(3, :)
            ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + tmpa
            ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) - tmpa


            tmp = (uterms(1) * t * eePowers(i, j, t-1) + &
                   uterms(2) * eePowers(i, j, t)) * rLapl(3) + &
                  uterms(3) * eePowers(i, j, t) + &
                  uterms(4) * t * eePowers(i, j, t-1) + &
                  uterms(1) * t * (t-1) * eePowers(i, j, t-2)
            uklapli(i, pos) = uklapli(i, pos) + tmp
            uklapli(j, pos) = uklapli(j, pos) + tmp
            uklapl(pos) = uklapl(pos) + 2 * tmp

          enddo
        endif
      endif
    enddo
  enddo

  !write(iul, *) "e-e corr complete"
  !---------Electron-Nucleus-Correlation-Terms----------------------------

  r = (-cutOffEN)**C
  v = C / cutOffEN
  do a = 1, nclast
    do i = 1, ne
      dist = rai(a, i) - cutOffEN
      if(dist < 0) then
        if(spinX == 1 .or. i <= nalpha) then
          spin = 1
        else
          spin = 2
        endif

        s = beta(spin, a, 1) * v
        p = -nucNum(a)/r + s
        w = (C * nucNum(a)/r - s)/cutOffEN

        sumr = beta(spin, a, 1) + p * rai(a, i)
        sumrd = p
        sumrn = 0

        do t = 2, xmax
          sumr = sumr + beta(spin, a, t) * enPowers(a, i, t)
          sumrd = sumrd + t * beta(spin, a, t) * enPowers(a, i, t-1)
          sumrn = sumrn + t * (t-1) * beta(spin, a, t) * enPowers(a, i, t-2)
        enddo

        rDeriv(1, 1) = xai(a, i) / rai(a, i)
        rDeriv(1, 2) = yai(a, i) / rai(a, i)
        rDeriv(1, 3) = zai(a, i) / rai(a, i)
        rLapl(1) = 2 / rai(a, i)

        distP(1) = dist
        do t = 2, C
          distP(t) = dist * distP(t-1)
        enddo

        xterms(1) = distP(C)
        xterms(2) = C * distP(C-1)
        xterms(3) = C * (C-1) * distP(C-2)
        xterms(4) = 2 * C * distP(C-1)
        xterms(5) = C * (C-1) * (C-2) * distP(C-3)
        u = xterms(1) * sumrd + xterms(2) * sumr

        xTerm = xTerm + xterms(1) * sumr
        xDeriv(3*i-2:3*i) = xDeriv(3*i-2:3*i) + u * rDeriv(1, :)
        xLapli(i) = xLapli(i) + &
                    u * rLapl(1) + &
                    xterms(3) * sumr + &
                    xterms(4) * sumrd + &
                    xterms(1) * sumrn

        if(deriveParams .and. xnum > 0) then
          if(curOptMode == 1) then
            pos = 1
            if(unum > 0) pos = 2
            ! cutoff derivatives
            q = w * rai(a, i)
            u = xterms(1) * w + &
                xterms(2) * q - &
                xterms(2) * sumrd - &
                xterms(3) * sumr

            uk(pos) = uk(pos) - xterms(2) * sumr + xterms(1) * q

            ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + u * rDeriv(1, :)

            tmp = u * rLapl(1) + &
                  xterms(3) * q + &
                  xterms(4) * w - &
                  xterms(5) * sumr - &
                  xterms(2) * sumrn - &
                  2 * xterms(3) * sumrd

            uklapli(i, pos) = uklapli(i, pos) + tmp
            uklapl(pos) = uklapl(pos) + tmp
          endif

          ! beta(1) derivatives
          q = 1 + v * rai(a, i)
          u = xterms(1) * v + xterms(2) * q
          pos = offset + umax * spinU + a + xnum * (spin - 1)

          uk(pos) = uk(pos) + xterms(1) * q

          ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + u * rDeriv(1, :)

          tmp = u * rLapl(1) + &
                xterms(3) * q + &
                xterms(4) * v
          uklapli(i, pos) = uklapli(i, pos) + tmp
          uklapl(pos) = uklapl(pos) + tmp

          do t = 2, xmax
            pos = offset + umax * spinU + (t-1) * nclast + a + xnum * (spin - 1)

            uk(pos) = uk(pos) + xterms(1) * enPowers(a, i, t)

            ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + &
                                     (xterms(2) * enPowers(a, i, t) + &
                                      xterms(1) * t * enPowers(a, i, t-1)) * rDeriv(1, :)

            tmp = (xterms(1) * t * enPowers(a, i, t-1) + &
                   xterms(2) * enPowers(a, i, t)) * rLapl(1) + &
                  xterms(3) * enPowers(a, i, t) + &
                  xterms(4) * t * enPowers(a, i, t-1) + &
                  xterms(1) * t * (t-1) * enPowers(a, i, t-2)
            uklapli(i, pos) = uklapli(i, pos) + tmp
            uklapl(pos) = uklapl(pos) + tmp
          enddo
        endif
      endif
    enddo
  enddo

  !write(iul, *) "e-n corr complete"
  !---------Electron-Electron-Nucleus-Correlation-Terms-------------------
  if(fnum > 0) then
    do a = 1, nclast
      do i = 1, ne
        do j = i + 1, ne
          dist = rai(a, i) - cutOffEEN

          if(dist < 0) then
            distP(1) = dist
            do t = 2, C
              distP(t) = dist * distP(t-1)
            enddo

            distj = rai(a, j) - cutOffEEN

            if(i < j .and. distj < 0) then
              if(spinF == 1 .or. ((i <= nalpha) .eqv. (j <= nalpha))) then
                spin = 1
              else
                spin = 2
              endif

              sumr = 0

              distjP(1) = distj
              do t = 2, C
                distjP(t) = distj * distjP(t-1)
              enddo

              fsum = 0

              do l = 0, fenmax
                do m = 0, fenmax
                  do n = 0, feemax
                    g = gamma(spin, a, l, m, n)

                    sumr = sumr + g * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n)

                    fsum(1) = fsum(1) + g * l * enPowers(a, i, l-1) * enPowers(a, j, m)   * eePowers(i, j, n)
                    fsum(2) = fsum(2) + g * m * enPowers(a, i, l)   * enPowers(a, j, m-1) * eePowers(i, j, n)
                    fsum(3) = fsum(3) + g * n * enPowers(a, i, l)   * enPowers(a, j, m)   * eePowers(i, j, n-1)
                    fsum(4) = fsum(4) + g * l * (l-1) * enPowers(a, i, l-2) * enPowers(a, j, m)   * eePowers(i, j, n)
                    fsum(5) = fsum(5) + g * m * (m-1) * enPowers(a, i, l)   * enPowers(a, j, m-2) * eePowers(i, j, n)
                    fsum(6) = fsum(6) + g * n * (n-1) * enPowers(a, i, l)   * enPowers(a, j, m)   * eePowers(i, j, n-2)
                    fsum(7) = fsum(7) + g * 2 * l * n * enPowers(a, i, l-1) * enPowers(a, j, m)   * eePowers(i, j, n-1)
                    fsum(8) = fsum(8) + g * 2 * m * n * enPowers(a, i, l) * enPowers(a, j, m-1)   * eePowers(i, j, n-1)
                  enddo
                enddo
              enddo

              rDeriv(1, 1) = xai(a, i) / rai(a, i)
              rDeriv(1, 2) = yai(a, i) / rai(a, i)
              rDeriv(1, 3) = zai(a, i) / rai(a, i)
              rLapl(1) = 2 / rai(a, i)

              rDeriv(2, 1) = xai(a, j) / rai(a, j)
              rDeriv(2, 2) = yai(a, j) / rai(a, j)
              rDeriv(2, 3) = zai(a, j) / rai(a, j)
              rLapl(2) = 2 / rai(a, j)

              rDeriv(3, 1) = xij(i, j) / rij(i, j)
              rDeriv(3, 2) = yij(i, j) / rij(i, j)
              rDeriv(3, 3) = zij(i, j) / rij(i, j)
              rLapl(3) = 2 / rij(i, j)

              rComb(1) = sum(rDeriv(1, :) * rDeriv(3, :))
              rComb(2) = sum(rDeriv(2, :) * rDeriv(3, :))

              fterms(1) = distP(C) * distjP(C)
              fterms(2) = C * distP(C-1) * distjP(C)
              fterms(3) = C * distP(C) * distjP(C-1)
              fterms(4) = C * (C-1) * distP(C-2) * distjP(C)
              fterms(5) = C * (C-1) * distP(C) * distjP(C-2)
              fterms(6) = 2 * C * distP(C-1) * distjP(C)
              fterms(7) = 2 * C * distP(C) * distjP(C-1)
              fterms(8) = C * C * distP(C-1) * distjP(C-1)

              fTerm = fTerm + fterms(1) * sumr

              fDeriv(3*i-2:3*i) = fDeriv(3*i-2:3*i) + &
                                  fterms(1) * fsum(1) * rDeriv(1, 1:3) + &
                                  fterms(1) * fsum(3) * rDeriv(3, 1:3) + &
                                  fterms(2) * sumr * rDeriv(1, 1:3)
              fDeriv(3*j-2:3*j) = fDeriv(3*j-2:3*j) + &
                                  fterms(1) * fsum(2) * rDeriv(2, 1:3) - &
                                  fterms(1) * fsum(3) * rDeriv(3, 1:3) + &
                                  fterms(3) * sumr * rDeriv(2, 1:3)

              fsuml(1) = fsum(1) * rLapl(1) + &
                         fsum(3) * rLapl(3) + &
                         fsum(4) + &
                         fsum(6) + &
                         fsum(7) * rComb(1)
              fsuml(2) = fsum(2) * rLapl(2) + &
                         fsum(3) * rLapl(3) + &
                         fsum(5) + &
                         fsum(6) - &
                         fsum(8) * rComb(2)

              fLapli(i) = fLapli(i) + &
                          fterms(1) * fsuml(1) + &
                          fterms(2) * sumr * rLapl(1) + &
                          fterms(4) * sumr + &
                          fterms(6) * (fsum(1) + fsum(3) * rComb(1))
              fLapli(j) = fLapli(j) + &
                          fterms(1) * fsuml(2) + &
                          fterms(3) * sumr * rLapl(2) + &
                          fterms(5) * sumr + &
                          fterms(7) * (fsum(2) - fsum(3) * rComb(2))

              if(deriveParams) then
                if(curOptMode == 1) then
                 pos = 1
                 if(unum > 0) pos = pos + 1
                 if(xnum > 0) pos = pos + 1

                  ! cutoff derivatives
                  uk(pos) = uk(pos) - (fterms(2) + fterms(3)) * sumr

                  ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) - &
                                         (fterms(4) + fterms(8)) * sumr * rDeriv(1, 1:3) - &
                                         (fterms(2) + fterms(3)) * fsum(1) * rDeriv(1, 1:3) - &
                                         (fterms(2) + fterms(3)) * fsum(3) * rDeriv(3, 1:3)
                  ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) - &
                                         (fterms(5) + fterms(8)) * sumr * rDeriv(2, 1:3) - &
                                         (fterms(2) + fterms(3)) * fsum(2) * rDeriv(2, 1:3) + &
                                         (fterms(2) + fterms(3)) * fsum(3) * rDeriv(3, 1:3)

                  tmp = (fterms(2) + fterms(3)) * fsuml(1) + &
                        (fterms(4) + fterms(8)) * sumr * rLapl(1) + &
                        (C * (C-1) * (C-2) * distP(C-3) * distjP(C) + &
                         C * C * (C-1) * distP(C-2) * distjP(C-1)) * sumr + &
                        2 * (fterms(4) + fterms(8)) * (fsum(1) + fsum(3) * rComb(1))
                  uklapli(i, pos) = uklapli(i, pos) - tmp
                  uklapl(pos) = uklapl(pos) - tmp

                  tmp = (fterms(2) + fterms(3)) * fsuml(2) + &
                        (fterms(5) + fterms(8)) * sumr * rLapl(2) + &
                        (C * (C-1) * (C-2) * distP(C) * distjP(C-3) + &
                         C * C * (C-1) * distP(C-1) * distjP(C-2)) * sumr + &
                        2 * (fterms(5) + fterms(8)) * (fsum(2) - fsum(3) * rComb(2))
                  uklapli(j, pos) = uklapli(j, pos) - tmp
                  uklapl(pos) = uklapl(pos) - tmp
                endif

                do l = 0, fenmax
                  do m = 0, fenmax
                    do n = 0, feemax
                      pos = offset + umax * spinU + xnum * spinX + fnum * (spin - 1) + oneDim(l, m, n, a)

                      fsump(1) = l * enPowers(a, i, l-1) * enPowers(a, j, m) * eePowers(i, j, n)
                      fsump(2) = m * enPowers(a, i, l) * enPowers(a, j, m-1) * eePowers(i, j, n)
                      fsump(3) = n * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n-1)
                      fsump(4) = l * (l-1) * enPowers(a, i, l-2) * enPowers(a, j, m) * eePowers(i, j, n)
                      fsump(5) = m * (m-1) * enPowers(a, i, l) * enPowers(a, j, m-2) * eePowers(i, j, n)
                      fsump(6) = n * (n-1) * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n-2)
                      fsump(7) = 2 * l * n * enPowers(a, i, l-1) * enPowers(a, j, m) * eePowers(i, j, n-1)
                      fsump(8) = 2 * m * n * enPowers(a, i, l) * enPowers(a, j, m-1) * eePowers(i, j, n-1)

                      sumpr = enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n)

                      uk(pos) = uk(pos) + fterms(1) * sumpr

                      ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + &
                                               fterms(1) * fsump(1) * rDeriv(1, 1:3) + &
                                               fterms(1) * fsump(3) * rDeriv(3, 1:3) + &
                                               fterms(2) * sumpr * rDeriv(1, 1:3)
                      ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) + &
                                               fterms(1) * fsump(2) * rDeriv(2, 1:3) - &
                                               fterms(1) * fsump(3) * rDeriv(3, 1:3) + &
                                               fterms(3) * sumpr * rDeriv(2, 1:3)

                      fsumpl(1) = fsump(1) * rLapl(1) + &
                                  fsump(3) * rLapl(3) + &
                                  fsump(4) + &
                                  fsump(6) + &
                                  fsump(7) * rComb(1)
                      fsumpl(2) = fsump(2) * rLapl(2) + &
                                  fsump(3) * rLapl(3) + &
                                  fsump(5) + &
                                  fsump(6) - &
                                  fsump(8) * rComb(2)

                      tmp = fterms(1) * fsumpl(1) + &
                            fterms(2) * sumpr * rLapl(1) + &
                            fterms(4) * sumpr + &
                            fterms(6) * (fsump(1) + fsump(3) * rComb(1))
                      uklapli(i, pos) = uklapli(i, pos) + tmp
                      uklapl(pos) = uklapl(pos) + tmp

                      tmp = fterms(1) * fsumpl(2) + &
                            fterms(3) * sumpr * rLapl(2) + &
                            fterms(5) * sumpr + &
                            fterms(7) * (fsump(2) - fsump(3) * rComb(2))
                      uklapli(j, pos) = uklapli(j, pos) + tmp
                      uklapl(pos) = uklapl(pos) + tmp

                    enddo
                  enddo
                enddo
              endif
            endif

          endif
        enddo
      enddo
    enddo
  endif

  !write(iul, *) "e-e-n complete"
  !---------Calculation of U and its Deriviatives-------------------------
  ju = uTerm + xTerm + fTerm
  jud = uDeriv + xDeriv + fDeriv
  julapli = uLapli + xLapli + fLapli
  julapl = sum(julapli(1:ne))

  !write(iul, *) uTerm, xTerm, fTerm
  !write(iul, *) ju, julapl

  !---------Parameter derivatives-----------------------------------------

end subroutine jasdtnall_naiv

!=======================================================================

! calculates U without derivatives
subroutine jasdtnp(init, ie, rai, rij, ju)
  logical, intent(in) ::  init   ! .true. for initialization
  integer, intent(in) ::  ie     ! electron with new position
  real(r8), intent(in)  ::  rai(:, :), rij(:, :)   ! current distances
  real(r8), intent(out) ::  ju     ! returns U rather than exp(U)
  ju = 0
  call assert(.false., "jasdtnp called")
end subroutine jasdtnp

!=======================================================================

subroutine getVectorLenDTN(optMode,npJ1,npJ2,npJnl)
  integer, intent(in)    :: optMode
  integer, intent(inout) :: npJ1     ! one-electron linear
  integer, intent(inout) :: npJ2     ! two-electron linear
  integer, intent(inout) :: npJnl    ! nonlinear

  if(optMode < 5) curOptMode = optMode !XXX
  select case(optMode)
  case(1) ! all parameters
    npJ1 = xnum * spinX 
    npJ2 = unum * spinU + fnum * spinF
    npJnl = getOffset() 
  case(2) ! linear parameters
    npJ1 = xnum * spinX 
    npJ2 = unum * spinU + fnum * spinF
    npJnl = 0
  case(5) ! gamma parameters for cusp
    npJ1 = 0
    npJ2 = fnum * spinF
    npJnl = 0
  case default
    call abortp("getVectorLenDTN: this optMode is not implemented")
  end select
end subroutine getVectorLenDTN

!=======================================================================

subroutine getVectorDTN(optMode, p, spin)
  integer, intent(in) ::  optMode       ! optimization mode, 1 for linear params,
                                        !                    2 for nonlinear,
                                        !                    3 for all params,
                                        !                    5 for gamma only
  real(r8), intent(out) ::  p(:)          ! parameter vector
  integer, intent(in), optional :: spin
  integer :: s, pos, o

  o = 0
  if(optMode < 5) curOptMode = optMode !XXX
  select case(optMode)
  case(1) ! all parameters
    if(unum > 0) then
      o = o + 1
      p(o) = cutOffEE
    endif
    if(xnum > 0) then
      o = o + 1
      p(o) = cutOffEN
    endif
    if(fnum > 0) then
      o = o + 1
      p(o) = cutOffEEN
    endif

    do s = 1, spinU
      pos = o + unum * (s - 1)
      p(pos + 1:pos + unum) = alpha(s, 1:umax)
    enddo
    do s = 1, spinX
      pos = o + unum * spinU + xnum * (s - 1)
      p(pos + 1:pos + xnum) = reshape(beta(s, 1:nclast, 1:xmax), [xnum])
    enddo
    do s = 1, spinF
      pos = o + unum * spinU + xnum * spinX + fnum * (s - 1)
      p(pos + 1:pos + fnum) = reshape(gamma(s, 1:nclast, 0:fenmax, 0:fenmax, 0:feemax), [fnum])
    enddo
  case(2) ! linear parameters
    do s = 1, spinU
      pos = unum * (s - 1)
      p(pos + 1:pos + unum) = alpha(s, 1:umax)
    enddo
    do s = 1, spinX
      pos = unum * spinU + xnum * (s - 1)
      p(pos + 1:pos + xnum) = reshape(beta(s, 1:nclast, 1:xmax), [xnum])
    enddo
    do s = 1, spinF
      pos =  unum * spinU + xnum * spinX + fnum * (s - 1)
      p(pos + 1:pos + fnum) = reshape(gamma(s, 1:nclast, 0:fenmax, 0:fenmax, 0:feemax), [fnum])
    enddo
  case(5) ! gamma parameters for cusp
    do s = 1, spinF
      pos =  fnum * (s - 1)
      p(pos + 1:pos + fnum) = reshape(gamma(spin, 1:nclast, 0:fenmax, 0:fenmax, 0:feemax), [fnum])
    enddo
  case default
    call abortp("getVectorDTN: optMode not implemeted")
  end select
end subroutine getVectorDTN

!=======================================================================

subroutine putVectorDTN(optMode, p)
  integer, intent(in) ::  optMode       ! optimization mode
  real(r8), intent(in)  ::  p(:)          ! parameter vector
  integer :: s, pos, o

  o = 0
  if(optMode < 5) curOptMode = optMode !XXX
  select case(optMode)
  case(1) ! linear and nonlinear parameters
    if(unum > 0) then
      o = o + 1
      cutOffEE = p(o)
    endif
    if(xnum > 0) then
      o = o + 1
      cutOffEN = p(o)
    endif
    if(fnum > 0) then
      o = o + 1
      cutOffEEN = p(o)
    endif

    do s = 1, spinU
      pos = o + unum * (s - 1)
      alpha(s, 1:umax) = p(pos + 1:pos + unum)
    enddo
    do s = 1, spinX
      pos = o + unum * spinU + xnum * (s - 1)
      beta(s, 1:nclast, 1:xmax) = reshape(p(pos + 1:pos + xnum), [nclast, xmax])
    enddo
    do s = 1, spinF
      pos = o + unum * spinU + xnum * spinX + fnum * (s - 1)
      gamma(s, 1:nclast, 0:fenmax, 0:fenmax, 0:feemax) = reshape(p(pos + 1: pos + fnum), &
                                                                 [nclast, fenmax + 1, fenmax + 1, feemax + 1])
    enddo
  case(2) ! linear parameters
    do s = 1, spinU
      pos = unum * (s - 1)
      alpha(s, 1:umax) = p(pos + 1:pos + unum)
    enddo
    do s = 1, spinX
      pos = unum * spinU + xnum * (s - 1)
      beta(s, 1:nclast, 1:xmax) = reshape(p(pos + 1:pos + xnum), [nclast, xmax])
    enddo
    do s = 1, spinF
      pos = unum * spinU + xnum * spinX + fnum * (s - 1)
      gamma(s, 1:nclast, 0:fenmax, 0:fenmax, 0:feemax) = reshape(p(pos + 1: pos + fnum), &
                                                                 [nclast, fenmax + 1, fenmax + 1, feemax + 1])
    enddo
  case(5) ! gamma parameters for cusp
    do s = 1, spinF
      pos =  fnum * (s - 1)
      gamma(s, 1:nclast, 0:fenmax, 0:fenmax, 0:feemax) = reshape(p(pos + 1:pos + fnum), &
                                                                 [nclast, fenmax + 1, fenmax + 1, feemax + 1])
    enddo
  case default
    call abortp("putVectorDTN: optMode not implemeted")
  end select

  !call checkCusp
  !call jasoutput_dtn
end subroutine putVectorDTN

!=======================================================================

subroutine checkCusp()
  real(r8) :: cuspMatrix(0:3 * fenmax + feemax + 1, fnum)
  real(r8) :: paramVector(fnum)
  real(r8) :: resVector(3 * fenmax + feemax + 2)
  real(r8) :: ipiv(fnum)
  integer :: conditions
  integer :: k, l, m, n, a, r, s, info

  if(mytid /= 0) return
  cuspMatrix = 0
  resVector = 0
  conditions = 3 * fenmax + feemax + 2

  if(fnum < 1) return

  !TODO improve if checks
  if(feemax >= 1) then
    do k = 0, 2 * fenmax - 1
      do l = k/2 + 1, k
        m = k - l
        do a = 1, nclast
          cuspMatrix(k, oneDim(l, m, 1, a)) = 2
        enddo
      enddo
      if(mod(k, 2) == 0) then
        l = k/2
        do a = 1, nclast
          cuspMatrix(k, oneDim(l, l, 1, a)) = 1
        enddo
      endif
    enddo
  endif


  do a = 1, nclast
    do k = 0, feemax
      r = k + 2 * fenmax
      cuspMatrix(r, oneDim(0, 0, k, a)) = C
      if(1 <= fenmax) then
        cuspMatrix(r, oneDim(1, 0, k, a)) = -cutOffEEN
      endif
    enddo
    do k = 0, fenmax + feemax
      do l = 1, k
        n = k - l
        if(l <= fenmax .and. n <= feemax) then
          cuspMatrix(r, oneDim(l, 0, n, a)) = C
          if(1 <= fenmax) then
            cuspMatrix(r, oneDim(l, 1, n, a)) = -cutOffEEN
          endif
        endif
      enddo
    enddo
  enddo

  !do a = 1, nclast
  !  do k = 0, fenmax + feemax
  !    r = k + 2 * fenmax
  !    if(k <= feemax) then
  !      cuspMatrix(r, oneDim(0, 0, k, a)) = C
  !      if(1 <= fenmax) then
  !        cuspMatrix(r, oneDim(1, 0, k, a)) = -cutOffEEN
  !      endif
  !    endif
  !    do l = 1, k
  !      n = k - l
  !      if(l <= fenmax .and. n <= feemax) then
  !        cuspMatrix(r, oneDim(l, 0, n, a)) = C
  !        if(1 <= fenmax) then
  !          cuspMatrix(r, oneDim(l, 1, n, a)) = -cutOffEEN
  !        endif
  !      endif
  !    enddo
  !  enddo
  !enddo

  do s = 1, spinF
    call getVectorDTN(5, paramVector, s)
    call dgemv('N', conditions, fnum, 1d0, cuspMatrix(0,1), conditions, paramVector, 1, 0, resVector, 1)
    write(iul,*) "Cusp check for spin function", s
    write(iul,*) resVector
  enddo
end subroutine checkCusp

!=======================================================================

integer function oneDim(l, m, n, a)
  integer :: l, m, n, a
  oneDim = n * (fenmax + 1) * (fenmax + 1) * nclast + &
           m * (fenmax + 1) * nclast + &
           l * nclast + &
           a
end function oneDim

integer function getOffset()
  integer :: o

  o = 0
  if(unum > 0) o = o + 1
  if(xnum > 0) o = o + 1
  if(fnum > 0) o = o + 1

  getOffset = o
end function getOffset

integer function getNumConstraintsDTN(optMode)
  integer, intent(in) :: optMode
  getNumConstraintsDTN = (3 * fenmax + feemax + 2) * spinF
end function getNumConstraintsDTN

subroutine getConstraintsDTN(optMode, nConstraints, constraints, &
                             constraintsGrad, constraintsHessian)
  integer, intent(in) :: optMode
  integer, intent(in) :: nConstraints
  real(r8), intent(inout) :: constraints(:)
  real(r8), intent(inout) :: constraintsGrad(:, :)
  real(r8), intent(inout) :: constraintsHessian(:, :, :)

  integer :: a, k, l, m, n, r, s
  integer :: base, offset, poffset, pos

  ! TODO make gamma symmetric (see appendix A) instead?
  constraints = 0
  constraintsGrad = 0
  constraintsHessian = 0

  base = getOffset() + unum * spinU + xnum * spinX

  do s = 1, spinF
    offset = (3 * fenmax + feemax + 2) * (s - 1)
    poffset = base + fnum * (s - 1)
    do k = 0, 2 * fenmax
      r = k + offset + 1
      do m = 0, k
        l = k - m
        do a = 1, nclast
          constraints(r) = constraints(r) + gamma(s, a, l, m, 1)
          pos = poffset + oneDim(l, m, 1, a)
          constraintsGrad(r, pos) = constraintsGrad(r, pos) + 1
        enddo
      enddo
    enddo

    do k = 0, fenmax + feemax
      r = k + 2 * fenmax + offset + 1
      do m = 0, min(k, fenmax) ! XXX
        n = k - m
        if(n <= feemax) then ! XXX
          do a = 1, nclast
            constraints(r) = constraints(r) + C * gamma(s, a, 0, m, n) &
                                            - cutOffEEN * gamma(s, a, 1, m, n)
            pos = poffset + oneDim(0, m, n, a)
            constraintsGrad(r, pos) = constraintsGrad(r, pos) + C
            pos = poffset + oneDim(1, m, n, a)
            constraintsGrad(r, pos) = constraintsGrad(r, pos) - cutOffEEN
          enddo
        endif
      enddo
    enddo
  enddo

end subroutine getConstraintsDTN


!=======================================================================


subroutine jasgeneric_naiv(x, y, z, rai, rij, optType, ju, jud, julapl, julapli)
  ! main calculations for all electrons moved

  real(r8), intent(in)  ::  x(:), y(:), z(:)     ! cartesian coords
  real(r8), intent(in)  ::  rai(:, :), rij(:, :) ! distances r_ai, r_ij
  character(len=*)    ::  optType              ! parameter optimization
  real(r8), intent(out) ::  ju                   ! U
  real(r8), intent(out) ::  jud(:)               ! nabla U
  real(r8), intent(out) ::  julapl               ! sum of laplacians
  real(r8), intent(out) ::  julapli(:)           ! laplacian U

  ! values of correlation terms
  real(r8) :: uTerm, xTerm, fTerm
  ! nabla of correlation terms
  real(r8) :: uDeriv(3*ne), xDeriv(3*ne), fDeriv(3*ne)
  ! laplacian
  real(r8) :: uLapli(ne), xLapli(ne), fLapli(ne)

  ! parameter derivatives
  logical :: deriveParams

  ! spin coefficient for e-e-term
  real(r8) :: spinCoeff(ne, ne)
  ! numerator in e-n term, used to enforce cusp if necessary
  real(r8) :: nucNum(nclast)

  ! x/y/z distances e-n and e-e
  real(r8) :: xai(nclast, ne), yai(nclast, ne), zai(nclast, ne)
  real(r8) :: xij(ne, ne), yij(ne, ne), zij(ne, ne)

  ! powers of distances
  real(r8) :: eePowers(ne, ne, -2:max(umax, feemax, 1))
  real(r8) :: enPowers(nclast, ne, -2:max(xmax, fenmax, 1))

  ! derivatives of distances
  real(r8) :: rijDeriv(ne, ne, 3), rijDerivj(ne, ne, 3), raiDeriv(nclast, ne, 3)
  ! (dr/dx)^2 + (dr/dy)^2 + (dr/dz)^2
  real(r8) :: rijSquare(ne, ne), rijSquarej(ne, ne), raiSquare(nclast, ne)
  ! laplacian
  real(r8) :: rijLapl(ne, ne), rijLaplj(ne, ne), raiLapl(nclast, ne)
  real(r8) :: rComb(3)

  ! (distance - cut-off length) and powers
  real(r8) :: dist, distj, distP(-1:C), distjP(-1:C)
  ! sum in x, u and f-term, derivative, nabla
  real(r8) :: sumr, sumrd, sumrn, sumpr
  ! commonly used expressions in sum/derivative/nabla/laplacian of u/x/f
  real(r8) :: uterms(4), xterms(4), fterms(12)
  real(r8) :: fsum(8), fsuml(2)
  ! for parameter derivatives
  real(r8) :: fsump(8), fsumpl(2)
  ! position of parameter for param derivatives
  integer :: pos, offset, coOffset
  integer :: spin

  ! commonly used terms
  real(r8) :: g, p, q, r, s, u, v, w, tmp, tmp2, tmpa(3)

  ! loop variables
  integer :: a, i, j, l, m, n, t, ll, mm, nn, xx, yy, zz, np1, np2, npnl

  ! alloc stat
  integer :: error

  !write(iul, *) "start"

  !---------Initialization of arrays-------------------------------
  ju = 0
  jud = 0
  julapl = 0
  julapli = 0

  uTerm = 0
  xTerm = 0
  fTerm = 0
  uDeriv = 0
  xDeriv = 0
  fDeriv = 0
  uLapli = 0
  xLapli = 0
  fLapli = 0

  offset = 0
  deriveParams = .false.
  if(optType == "jastrow") then
    deriveParams = .true.
    call assert(allocated(uk), "jastrowParamData not allocated")
    call getVectorLenDTN(curOptMode,np1,np2,npnl)
    call assert(uk_params == np1+np2+npnl,"jasgeneric_naiv: mismatch of uk_params")
    if(curOptMode == 1) offset = getOffset()
    uk = 0
    ukgrad = 0
    uklapl = 0
    uklapli = 0
  endif

  !---------Precalculation of Terms---------------------------------------

  if(nucCusp) then
    nucNum = atoms(1:nclast)%za
  else
    nucNum = 0
  endif

  distP(-1) = 0
  distP(0) = 1
  distjP(-1) = 0
  distjP(0) = 1
  do i = 1, ne
    do a = 1, nclast
      ! x/y/z distances for e-n
      xai(a, i) = x(i) - atoms(a)%cx
      yai(a, i) = y(i) - atoms(a)%cy
      zai(a, i) = z(i) - atoms(a)%cz

      ! powers of e-n distances
      enPowers(a, i, -2) = 0
      enPowers(a, i, -1) = 0
      enPowers(a, i, 0) = 1
      enPowers(a, i, 1) = rai(a, i)
      do t = 2, max(xmax, fenmax)
        enPowers(a, i, t) = enPowers(a, i, t-1) * enPowers(a, i, 1)
      enddo

      raiDeriv(a, i, 1) = xai(a, i) / rai(a, i)
      raiDeriv(a, i, 2) = yai(a, i) / rai(a, i)
      raiDeriv(a, i, 3) = zai(a, i) / rai(a, i)

      raiSquare(a, i) = 1
      !raiSquare(a, i) = raiDeriv(a, i, 1) ** 2 + raiDeriv(a, i, 2) ** 2 + raiDeriv(a, i, 3) ** 2

      raiLapl(a, i) = 2 / rai(a, i)
    enddo
    do j = i + 1, ne
      ! x/y/z distance for e-e
      xij(i, j) = x(i) - x(j)
      yij(i, j) = y(i) - y(j)
      zij(i, j) = z(i) - z(j)

      ! powers of e-e distances
      eePowers(i, j, -2) = 0
      eePowers(i, j, -1) = 0
      eePowers(i, j, 0) = 1
      eePowers(i, j, 1) = rij(i, j)
      do t = 2, max(umax, feemax)
        eePowers(i, j, t) = eePowers(i, j, t-1) * eePowers(i, j, 1)
      enddo

      rijDeriv(i, j, 1) = xij(i, j) / rij(i, j)
      rijDeriv(i, j, 2) = yij(i, j) / rij(i, j)
      rijDeriv(i, j, 3) = zij(i, j) / rij(i, j)
      rijDerivj(i, j, 1) = -rijDeriv(i, j, 1)
      rijDerivj(i, j, 2) = -rijDeriv(i, j, 2)
      rijDerivj(i, j, 3) = -rijDeriv(i, j, 3)

      rijSquare(i, j) = 1
      rijSquarej(i, j) = 1

      rijLapl(i, j) = 2 / rij(i, j)
      rijLaplj(i, j) = 2 / rij(i, j)
    enddo
  enddo

  ! spin coefficients for e-e correlation
  spinCoeff(1:nalpha, 1:nalpha) = 0.25d0
  spinCoeff(1:nalpha, 1 + nalpha:ne) = 0.5d0
  spinCoeff(1 + nalpha:ne, 1 + nalpha:ne) = 0.25d0


  !---------Electron-Electron-Correlation Terms---------------------------

  r = (-cutOffEE)**C
  v = C / cutOffEE
  do i = 1, ne
    do j = i + 1, ne

      dist = eePowers(i, j, 1) - cutOffEE
      ! only calculate correlation if electron-electron distance is <= cutoff
      if(i < j .and. dist < 0) then
        distP(1) = dist
        do t = 2, C
          distP(t) = dist * distP(t-1)
        enddo

        if(spinU == 1 .or. ((i <= nalpha) .eqv. (j <= nalpha))) then
          spin = 1
        else
          spin = 2
        endif

        s = alpha(spin, 1) * v

        uterms(1) = distP(C)
        uterms(2) = C * distP(C-1)
        uterms(3) = C * (C-1) * distP(C-2)
        uterms(4) = C * (C-1) * (C-2) * distP(C-3)

        sumr = alpha(spin, 1) + (spinCoeff(i, j)/r + s) * eePowers(i, j, 1)
        sumrd = spinCoeff(i, j)/r + s
        sumrn = 0

        do t = 2, umax
          sumr = sumr + alpha(spin, t) * eePowers(i, j, t)
          sumrd = sumrd + t * alpha(spin, t) * eePowers(i, j, t-1)
          sumrn = sumrn + t * (t-1) * alpha(spin, t) * eePowers(i, j, t-2)
        enddo

        uTerm = uTerm + uterms(1) * sumr

        u = uterms(1) * sumrd + uterms(2) * sumr
        uDeriv(3*i-2:3*i) = uDeriv(3*i-2:3*i) + u * rijDeriv(i, j, :)
        uDeriv(3*j-2:3*j) = uDeriv(3*j-2:3*j) + u * rijDerivj(i, j, :)

        tmp = uterms(3) * sumr + uterms(1) * sumrn + 2 * uterms(2) * sumrd
        uLapli(i) = uLapli(i) + u * rijLapl(i, j) + tmp * rijSquare(i, j)
        uLapli(j) = uLapli(j) + u * rijLaplj(i, j) + tmp * rijSquarej(i, j)

        if(deriveParams .and. unum > 0) then
          if(curOptMode == 1) then
            ! cutoff derivatives
            p = - (C * spinCoeff(i, j)/r + s)/cutOffEE
            q = p * eePowers(i, j, 1)
            u = uterms(1) * p + &
                uterms(2) * q - &
                uterms(2) * sumrd - &
                uterms(3) * sumr

            uk(1) = uk(1) - uterms(2) * sumr + uterms(1) * q

            ukgrad(3*i-2:3*i, 1) = ukgrad(3*i-2:3*i, 1) + u * rijDeriv(i, j, :)
            ukgrad(3*j-2:3*j, 1) = ukgrad(3*j-2:3*j, 1) + u * rijDerivj(i, j, :)

            tmp = uterms(3) * q + &
                  2 * uterms(2) * p - &
                  uterms(4) * sumr - &
                  uterms(2) * sumrn - &
                  2 * uterms(3) * sumrd

            uklapli(i, 1) = uklapli(i, 1) + u * rijLapl(i, j) + tmp * rijSquare(i, j)
            uklapli(j, 1) = uklapli(j, 1) + u * rijLaplj(i, j) + tmp * rijSquarej(i, j)
            uklapl(1) = uklapl(1) + u * rijLapl(i, j) + tmp * rijSquare(i, j)
            uklapl(1) = uklapl(1) + u * rijLaplj(i, j) + tmp * rijSquarej(i, j)
          endif

          ! alpha(1) derivative
          pos = 1 + offset + unum * (spin - 1)

          w = 1 + v * eePowers(i, j, 1)
          u = uterms(1) * v + uterms(2) * w
          uk(pos) = uk(pos) + uterms(1) * w

          ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + u * rijDeriv(i, j, :)
          ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) + u * rijDerivj(i ,j, :)

          tmp = uterms(3) * w + 2 * uterms(2) * v
          uklapli(i, pos) = uklapli(i, pos) + u * rijLapl(i, j) + tmp * rijSquare(i, j)
          uklapli(j, pos) = uklapli(j, pos) + u * rijLaplj(i, j) + tmp * rijSquarej(i, j)
          uklapl(pos) = uklapl(pos) + u * rijLapl(i, j) + tmp * rijSquare(i, j)
          uklapl(pos) = uklapl(pos) + u * rijLaplj(i, j) + tmp * rijSquarej(i, j)

          do t = 2, umax
            pos = offset + t + umax * (spin - 1)
            uk(pos) = uk(pos) + uterms(1) * eePowers(i, j, t)

            u = uterms(1) * t * eePowers(i, j, t-1) + &
                uterms(2) * eePowers(i, j, t)
            ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + u * rijDeriv(i, j, :)
            ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) + u * rijDerivj(i, j, :)

            tmp = uterms(3) * eePowers(i, j, t) + &
                  2 * uterms(2) * t * eePowers(i, j, t-1) + &
                  uterms(1) * t * (t-1) * eePowers(i, j, t-2)
            uklapli(i, pos) = uklapli(i, pos) + u * rijLapl(i, j) + tmp * rijSquare(i, j)
            uklapli(j, pos) = uklapli(j, pos) + u * rijLaplj(i, j) + tmp * rijSquarej(i, j)
            uklapl(pos) = uklapl(pos) + u * rijLapl(i, j) + tmp * rijSquare(i, j)
            uklapl(pos) = uklapl(pos) + u * rijLaplj(i, j) + tmp * rijSquarej(i, j)
          enddo
        endif
      endif
    enddo
  enddo

  !write(iul, *) "e-e corr complete"
  !---------Electron-Nucleus-Correlation-Terms----------------------------

  r = (-cutOffEN)**C
  v = C / cutOffEN
  do a = 1, nclast
    do i = 1, ne
      dist = enPowers(a, i, 1) - cutOffEN
      if(dist < 0) then
        distP(1) = dist
        do t = 2, C
          distP(t) = dist * distP(t-1)
        enddo

        if(spinX == 1 .or. i <= nalpha) then
          spin = 1
        else
          spin = 2
        endif

        s = beta(spin, a, 1) * v
        p = -nucNum(a)/r + s
        w = (C * nucNum(a)/r - s)/cutOffEN

        xterms(1) = distP(C)
        xterms(2) = C * distP(C-1)
        xterms(3) = C * (C-1) * distP(C-2)
        xterms(4) = C * (C-1) * (C-2) * distP(C-3)

        sumr = beta(spin, a, 1) + p * enPowers(a, i, 1)
        sumrd = p
        sumrn = 0

        do t = 2, xmax
          sumr = sumr + beta(spin, a, t) * enPowers(a, i, t)
          sumrd = sumrd + t * beta(spin, a, t) * enPowers(a, i, t-1)
          sumrn = sumrn + t * (t-1) * beta(spin, a, t) * enPowers(a, i, t-2)
        enddo

        xTerm = xTerm + xterms(1) * sumr

        u = (xterms(1) * sumrd + xterms(2) * sumr)
        xDeriv(3*i-2:3*i) = xDeriv(3*i-2:3*i) + u * raiDeriv(a, i, :)

        xLapli(i) = xLapli(i) + &
                    u * raiLapl(a, i) + &
                    (xterms(3) * sumr + 2 * xterms(2) * sumrd + xterms(1) * sumrn) * raiSquare(a, i)

        if(deriveParams .and. xnum > 0) then
          if(curOptMode == 1) then
            pos = 1
            if(unum > 0) pos = 2
            ! cutoff derivatives
            q = w * enPowers(a, i, 1)
            u = xterms(1) * w + &
                xterms(2) * q - &
                xterms(2) * sumrd - &
                xterms(3) * sumr

            uk(pos) = uk(pos) - xterms(2) * sumr + xterms(1) * q

            ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + u * raiDeriv(a, i, :)

            tmp = u * raiLapl(a, i) + &
                  (xterms(3) * q + &
                   2 * xterms(2) * w - &
                   xterms(4) * sumr - &
                   xterms(2) * sumrn - &
                   2 * xterms(3) * sumrd) * raiSquare(a, i)

            uklapli(i, pos) = uklapli(i, pos) + tmp
            uklapl(pos) = uklapl(pos) + tmp
          endif

          ! beta(1) derivatives
          q = 1 + v * enPowers(a, i, 1)
          u = xterms(1) * v + xterms(2) * q
          pos = offset + umax * spinU + a + xnum * (spin - 1)

          uk(pos) = uk(pos) + xterms(1) * q

          ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + u * raiDeriv(a, i, :)

          tmp = u * raiLapl(a, i) + &
                (xterms(3) * q + &
                 2 * xterms(2) * v) * raiSquare(a, i)
          uklapli(i, pos) = uklapli(i, pos) + tmp
          uklapl(pos) = uklapl(pos) + tmp

          do t = 2, xmax
            pos = offset + umax * spinU + (t-1) * nclast + a + xnum * (spin - 1)

            uk(pos) = uk(pos) + xterms(1) * enPowers(a, i, t)

            ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + &
                                     (xterms(2) * enPowers(a, i, t) + &
                                      xterms(1) * t * enPowers(a, i, t-1)) * raiDeriv(a, i, :)

            tmp = (xterms(1) * t * enPowers(a, i, t-1) + &
                   xterms(2) * enPowers(a, i, t)) * raiLapl(a, i) + &
                  (xterms(3) * enPowers(a, i, t) + &
                   2 * xterms(2) * t * enPowers(a, i, t-1) + &
                   xterms(1) * t * (t-1) * enPowers(a, i, t-2)) * raiSquare(a, i)
            uklapli(i, pos) = uklapli(i, pos) + tmp
            uklapl(pos) = uklapl(pos) + tmp
          enddo
        endif
      endif
    enddo
  enddo

  !write(iul, *) "e-n corr complete"
  !---------Electron-Electron-Nucleus-Correlation-Terms-------------------
  if(fnum > 0) then
    do a = 1, nclast
      do i = 1, ne
        do j = i + 1, ne
          dist = enPowers(a, i, 1) - cutOffEEN

          if(dist < 0) then
            distP(1) = dist
            do t = 2, C
              distP(t) = dist * distP(t-1)
            enddo

            distj = enPowers(a, j, 1) - cutOffEEN

            if(i < j .and. distj < 0) then
              if(spinF == 1 .or. ((i <= nalpha) .eqv. (j <= nalpha))) then
                spin = 1
              else
                spin = 2
              endif

              sumr = 0

              distjP(1) = distj
              do t = 2, C
                distjP(t) = distj * distjP(t-1)
              enddo

              fsum = 0

              do l = 0, fenmax
                do m = 0, fenmax
                  do n = 0, feemax
                    g = gamma(spin, a, l, m, n)

                    sumr = sumr + g * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n)

                    fsum(1) = fsum(1) + g * l * enPowers(a, i, l-1) * enPowers(a, j, m)   * eePowers(i, j, n)
                    fsum(2) = fsum(2) + g * m * enPowers(a, i, l)   * enPowers(a, j, m-1) * eePowers(i, j, n)
                    fsum(3) = fsum(3) + g * n * enPowers(a, i, l)   * enPowers(a, j, m)   * eePowers(i, j, n-1)
                    fsum(4) = fsum(4) + g * l * (l-1) * enPowers(a, i, l-2) * enPowers(a, j, m)   * eePowers(i, j, n)
                    fsum(5) = fsum(5) + g * m * (m-1) * enPowers(a, i, l)   * enPowers(a, j, m-2) * eePowers(i, j, n)
                    fsum(6) = fsum(6) + g * n * (n-1) * enPowers(a, i, l)   * enPowers(a, j, m)   * eePowers(i, j, n-2)
                    fsum(7) = fsum(7) + g * 2 * l * n * enPowers(a, i, l-1) * enPowers(a, j, m)   * eePowers(i, j, n-1)
                    fsum(8) = fsum(8) + g * 2 * m * n * enPowers(a, i, l)   * enPowers(a, j, m-1) * eePowers(i, j, n-1)
                  enddo
                enddo
              enddo

              rComb(1) = sum(raiDeriv(a, i, :) * rijDeriv(i, j, :))
              rComb(2) = sum(raiDeriv(a, j, :) * rijDerivj(i, j, :))

              fterms(1) = distP(C) * distjP(C)
              fterms(2) = C * distP(C-1) * distjP(C)
              fterms(3) = C * distP(C) * distjP(C-1)
              fterms(4) = C * (C-1) * distP(C-2) * distjP(C)
              fterms(5) = C * (C-1) * distP(C) * distjP(C-2)
              fterms(6) = 2 * C * distP(C-1) * distjP(C)
              fterms(7) = 2 * C * distP(C) * distjP(C-1)
              fterms(8) = C * C * distP(C-1) * distjP(C-1)
              fterms(9) = C * C * (C-1) * distP(C-2) * distjP(C-1)
              fterms(10)= C * C * (C-1) * distP(C-1) * distjP(C-2)
              fterms(11)= C * (C-1) * (C-2) * distP(C-3) * distjP(C)
              fterms(12)= C * (C-1) * (C-2) * distP(C) * distjP(C-3)

              fTerm = fTerm + fterms(1) * sumr

              fDeriv(3*i-2:3*i) = fDeriv(3*i-2:3*i) + &
                                  fterms(1) * fsum(1) * raiDeriv(a, i, :) + &
                                  fterms(1) * fsum(3) * rijDeriv(i, j, :) + &
                                  fterms(2) * sumr * raiDeriv(a, i, :)
              fDeriv(3*j-2:3*j) = fDeriv(3*j-2:3*j) + &
                                  fterms(1) * fsum(2) * raiDeriv(a, j, :) + &
                                  fterms(1) * fsum(3) * rijDerivj(i, j, :) + &
                                  fterms(3) * sumr * raiDeriv(a, j, :)

              fsuml(1) = fsum(1) * raiLapl(a, i) + &
                         fsum(3) * rijLapl(i, j) + &
                         fsum(4) * raiSquare(a, i) + &
                         fsum(6) * rijSquare(i, j) + &
                         fsum(7) * rComb(1)
              fsuml(2) = fsum(2) * raiLapl(a, j) + &
                         fsum(3) * rijLaplj(i, j) + &
                         fsum(5) * raiSquare(a, j) + &
                         fsum(6) * rijSquarej(i, j) + &
                         fsum(8) * rComb(2)

              fLapli(i) = fLapli(i) + &
                          fterms(1) * fsuml(1) + &
                          fterms(2) * sumr * raiLapl(a, i) + &
                          fterms(4) * sumr * raiSquare(a, i) + &
                          fterms(6) * (fsum(1) * raiSquare(a, i) + fsum(3) * rComb(1))
              fLapli(j) = fLapli(j) + &
                          fterms(1) * fsuml(2) + &
                          fterms(3) * sumr * raiLapl(a, j) + &
                          fterms(5) * sumr * raiSquare(a, j) + &
                          fterms(7) * (fsum(2) * raiSquare(a, j) + fsum(3) * rComb(2))

              if(deriveParams) then
                if(curOptMode == 1) then
                 pos = 1
                 if(unum > 0) pos = pos + 1
                 if(xnum > 0) pos = pos + 1

                  ! cutoff derivatives
                  uk(pos) = uk(pos) - (fterms(2) + fterms(3)) * sumr

                  ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) - &
                                           (fterms(4) + fterms(8)) * sumr * raiDeriv(a, i, :) - &
                                           (fterms(2) + fterms(3)) * fsum(1) * raiDeriv(a, i, :) - &
                                           (fterms(2) + fterms(3)) * fsum(3) * rijDeriv(i, j, :)
                  ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) - &
                                           (fterms(5) + fterms(8)) * sumr * raiDeriv(a, j, :) - &
                                           (fterms(2) + fterms(3)) * fsum(2) * raiDeriv(a, j, :) - &
                                           (fterms(2) + fterms(3)) * fsum(3) * rijDerivj(i, j, :)

                  tmp = (fterms(2) + fterms(3)) * fsuml(1) + &
                        (fterms(4) + fterms(8)) * sumr * raiLapl(a, i) + &
                        (fterms(9) + fterms(11)) * sumr * raiSquare(a, i) + &
                        2 * (fterms(4) + fterms(8)) * (fsum(1) * raiSquare(a, i) + fsum(3) * rComb(1))
                  uklapli(i, pos) = uklapli(i, pos) - tmp
                  uklapl(pos) = uklapl(pos) - tmp

                  tmp = (fterms(2) + fterms(3)) * fsuml(2) + &
                        (fterms(5) + fterms(8)) * sumr * raiLapl(a, j) + &
                        (fterms(10) + fterms(12)) * sumr * raiSquare(a, j) + &
                        2 * (fterms(5) + fterms(8)) * (fsum(2) * raiSquare(a, j) + fsum(3) * rComb(2))
                  uklapli(j, pos) = uklapli(j, pos) - tmp
                  uklapl(pos) = uklapl(pos) - tmp
                endif

                ! gamma derivatives

                do l = 0, fenmax
                  do m = 0, fenmax
                    do n = 0, feemax
                      pos = offset + umax * spinU + xnum * spinX + fnum * (spin - 1) + oneDim(l, m, n, a)

                      fsump(1) = l * enPowers(a, i, l-1) * enPowers(a, j, m) * eePowers(i, j, n)
                      fsump(2) = m * enPowers(a, i, l) * enPowers(a, j, m-1) * eePowers(i, j, n)
                      fsump(3) = n * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n-1)
                      fsump(4) = l * (l-1) * enPowers(a, i, l-2) * enPowers(a, j, m) * eePowers(i, j, n)
                      fsump(5) = m * (m-1) * enPowers(a, i, l) * enPowers(a, j, m-2) * eePowers(i, j, n)
                      fsump(6) = n * (n-1) * enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n-2)
                      fsump(7) = 2 * l * n * enPowers(a, i, l-1) * enPowers(a, j, m) * eePowers(i, j, n-1)
                      fsump(8) = 2 * m * n * enPowers(a, i, l) * enPowers(a, j, m-1) * eePowers(i, j, n-1)

                      sumpr = enPowers(a, i, l) * enPowers(a, j, m) * eePowers(i, j, n)

                      uk(pos) = uk(pos) + fterms(1) * sumpr

                      ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + &
                                               fterms(1) * fsump(1) * raiDeriv(a, i, :) + &
                                               fterms(1) * fsump(3) * rijDeriv(i, j, :) + &
                                               fterms(2) * sumpr * raiDeriv(a, i, :)
                      ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) + &
                                               fterms(1) * fsump(2) * raiDeriv(a, j, :) + &
                                               fterms(1) * fsump(3) * rijDerivj(i, j, :) + &
                                               fterms(3) * sumpr * raiDeriv(a, j, :)

                      fsumpl(1) = fsump(1) * raiLapl(a, i) + &
                                  fsump(3) * rijLapl(i, j) + &
                                  fsump(4) * raiSquare(a, i) + &
                                  fsump(6) * rijSquare(i, j) + &
                                  fsump(7) * rComb(1)
                      fsumpl(2) = fsump(2) * raiLapl(a, j) + &
                                  fsump(3) * rijLaplj(i, j) + &
                                  fsump(5) * raiSquare(a, j) + &
                                  fsump(6) * rijSquarej(i, j) + &
                                  fsump(8) * rComb(2)

                      tmp = fterms(1) * fsumpl(1) + &
                            fterms(2) * sumpr * raiLapl(a, i) + &
                            fterms(4) * sumpr * raiSquare(a, i) + &
                            fterms(6) * (fsump(1) * raiSquare(a, i) + fsump(3) * rComb(1))
                      uklapli(i, pos) = uklapli(i, pos) + tmp
                      uklapl(pos) = uklapl(pos) + tmp

                      tmp = fterms(1) * fsumpl(2) + &
                            fterms(3) * sumpr * raiLapl(a, j) + &
                            fterms(5) * sumpr * raiSquare(a, j) + &
                            fterms(7) * (fsump(2) * raiSquare(a, j) + fsump(3) * rComb(2))
                      uklapli(j, pos) = uklapli(j, pos) + tmp
                      uklapl(pos) = uklapl(pos) + tmp

                    enddo
                  enddo
                enddo
              endif
            endif

          endif
        enddo
      enddo
    enddo
  endif

  !write(iul, *) "e-e-n complete"
  !---------Calculation of U and its Deriviatives-------------------------

  ju = uTerm + xTerm + fTerm
  jud = uDeriv + xDeriv + fDeriv
  julapli = uLapli + xLapli + fLapli
  julapl = sum(julapli(1:ne))

  !write(iul, *) uTerm, xTerm, fTerm
  !write(iul, *) ju, julapl

  !---------Parameter derivatives-----------------------------------------

end subroutine jasgeneric_naiv

end module jastrowDTN_m
