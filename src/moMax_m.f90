! Copyright (C) 2019 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module moMax_m
   use kinds_m, only: r8
   use error_m, only: assert, error
   use global_m, only: iul, iull, getNElec, bohr2ang, baseName
   use parsing_m, only: getinta, getdbla, ifinda, getstra
   use posList_m, only: posVal_t, pos_VList_t, isSmaller, isGreater
   use minimizer_w_sing_module, only: minimizer_w_sing
   use minimizer_ws_factory_module, only: create_ws_minimizer
   use initialPositions_m, only: init_LMOSampling, destroy_LMOSampling, lmoPosition, getLMOSize
   use fctn_module, only: Function_t
   use aos_m, only: aocalc, ao1calc, aosplcalc, ao1splcalc
   use mos_m, only: mocalc, mo1calc, mat, mat1x, mat1y, mat1z
   use eConfigs_m, only: eConfigArray, eConfigArray_new, eConfigArray_set, eConfigArray_destroy
   use wfData_m, only: atoms, getNNuc, getNAlpha, getNBeta, spline, getNOrb
   use atom_m, only: atom_getPosition, atom_getDistance, atoms_getPositionMatrix

   implicit none
   private
   public :: moMax_run, moMax_plot, moMax_plot_plane, fctn_mo, mo_f, mo_fg, setMOIdx

   type, extends(Function_t) :: fctn_mo
      integer, private :: moIdx = 0
   contains
      procedure :: eval_fg => mo_fg
      procedure :: eval => mo_f
      procedure :: setMOIdx
   end type fctn_mo

   real(r8) :: mValueThreshold = 0.d0 ! ignore maxima (of mo**2) below this value

contains

   function cross_product(a, b) result(axb)
      real(r8) :: axb(3)
      real(r8), intent(in) :: a(3)
      real(r8), intent(in) :: b(3) 

      axb(1) = a(2)*b(3) - a(3)*b(2)
      axb(2) = a(3)*b(1) - a(1)*b(3)
      axb(3) = a(1)*b(2) - a(2)*b(1)
   end function cross_product  

   function mo_f(this, x) result(f)
      class(fctn_mo), intent(in)  :: this
      real(r8), intent(in) :: x(:)
      real(r8) :: f

      call assert(SIZE(x) == 3, "mo_f: illegal size")

      call mo_mo(x, this%moIdx, f)

      f = - abs(f)
   end function mo_f


   subroutine mo_fg(this, x, f, g, mask)
      class(fctn_mo), intent(in)  :: this
      real(r8), intent(in) :: x(:)
      real(r8), intent(out) :: f
      real(r8), intent(out) :: g(:)
      logical, intent(in), optional :: mask(SIZE(x))

      call assert(SIZE(x) == 3 .and. SIZE(g) == 3, "mo_fg: illegal size")

      call mo_mo(x, this%moIdx, f, g)
      if (f > 0) then
         f = - f
         g = - g
         if (PRESENT(mask)) where (mask) g = 0._r8
      end if
   end subroutine mo_fg


   subroutine setMOIdx(this, idx)
      class(fctn_mo), intent(inout)  :: this
      integer, intent(in) :: idx
      this%moIdx = idx
   end subroutine setMOIdx


   subroutine mo_mo(position, moIdx, value, gradient)
      real(r8), intent(in) :: position(3)
      integer, intent(in) :: moIdx
      real(r8), intent(out) :: value
      real(r8), optional, intent(out) :: gradient(3)
      real(r8) :: rrai(getNNuc(), getNElec(), 1)
      type(eConfigArray)  :: eca
      integer :: a

      real(r8) :: x(getNElec()), y(getNElec()), z(getNElec())

      ! filling rrai (electron-nucleus distances) only for electron 1
      rrai = 0._r8
      do a = 1, getNNuc()
         rrai(a, 1, 1) = NORM2(position - atom_getPosition(atoms(a)))
      enddo

      ! filling x, y, and z only for electron 1
      x = 0._r8; y = 0._r8; z = 0._r8;
      x(1) = position(1)
      y(1) = position(2)
      z(1) = position(3)

      ! calculating aos and mos for electron 1
      if (PRESENT(gradient)) then
         call eConfigArray_new(eca, getNElec(), 1)
         call eConfigArray_set(eca, 1, x, y, z)

         if (spline) then
            call aosplcalc(1, eca, rrai)
         else
            call aocalc(1, eca, rrai)
         end if
         call mocalc(1) ! calculates mat, mat1x, mat1y, mat1z

         gradient(1) = mat1x(moIdx, 1, 1) 
         gradient(2) = mat1y(moIdx, 1, 1) 
         gradient(3) = mat1z(moIdx, 1, 1) 

         call eConfigArray_destroy(eca)
      else
         if (spline) then
            call ao1splcalc(1, x, y, z, rrai(:,:,1))
         else
            call ao1calc(1, x, y, z, rrai(:,:,1))
         end if
         call mo1calc(1) ! calculates mat
      end if

      value = mat(moIdx, 1, 1)
   end subroutine mo_mo


   subroutine moMax_run(lines, nl)
      ! sample LMOs using ellipsoids with scale factor 'centroidScale' followed by minimizing -phi**2
      character(len=120), intent(in) :: lines(:)
      integer, intent(in) :: nl
      integer :: i, ii, mo, iflag, maxCount, n, io, refThreshold
      character(len=10) :: mode
      type(posVal_t), pointer :: pvp
      type(pos_VList_t), allocatable :: posList(:)
      integer, parameter :: iu = 13

      if (getNOrb() <= 1) call error("$maximize_mos: no orbitals available")
      allocate(posList(getNOrb()))

      mValueThreshold = 1.0d-4
      call getdbla(lines, nl, 'value_threshold=', mValueThreshold, iflag)

      refThreshold = 0
      call getinta(lines, nl, 'ref_threshold=', refThreshold, iflag)

      mode = 'lmo_min'
      call getstra(lines, nl, "mode=", mode, iflag)

      if (mode == 'lmo_min') then
         call moMax_lmo_min(lines, nl, posList)
      else if (mode == 'grid') then
         call moMax_grid(lines, nl, posList)
      else
         call error("$maximize_mos: illegal mode (lmo_min|grid)")
      end if

      ! print maxima from posList
      maxCount = 0
      do mo = 1, size(posList)
         if (posList(mo)%size() == 0) cycle
         maxCount = max(maxCount, posList(mo)%size())
         call posList(mo)%sort(isGreater)
         write(iul,'(/a,i5)') ' maxima of squared mo #', mo
         do i = 1, posList(mo)%size()
            pvp => posList(mo)%elem(i)
            write(iul, '(i5,g16.6,a,3g15.4,a,i8)') i, -pvp%value, ' pos:', bohr2ang * pvp%pos, &
               ' count:', pvp%count
         end do
      end do

      ! write mo maxima in ref-file

      open(iu,file=trim(baseName)//'.ref', iostat=io)
      if (io/=0) call error('(maximize_mos): opening file '//trim(baseName)//'.ref failed')
      write(iu,'(i5,a)') getNNuc(), " nuclei:"
      do i = 1, getNNuc()
         write(iu,'(i4,1X,a2,1X,3f14.7)') i, atoms(i)%elem, atoms(i)%cx, atoms(i)%cy, atoms(i)%cz
      enddo

      write(iu,'(i5,1x,a)') (maxCount + 1) / 2, "MOMAX"

      do n = 1, (maxCount + 1) / 2
         write(iu,'(a,i5,a,i4,a,i4)') "MOMAX:", n, "  1 F(MOMAX):  0  found:  0  0", &
            getNAlpha(), "    0   ", getNelec()
         write(iu,'(i5)') 2 * size(posList)
         do mo = 1, size(posList)
            ii = 0
            do i = 1, min(posList(mo)%size(), 2 * n - 1)
               pvp => posList(mo)%elem(i)
               if (pvp%count < refThreshold) then
                  ii = i - 1
                  exit
               else
                  ii = i
               end if 
               ii = max(ii, 1)
            end do
            pvp => posList(mo)%elem(ii)
            write(iu,'(3f13.6)') pvp%pos(1), pvp%pos(2), pvp%pos(3)
         end do
         do mo = 1, size(posList)
            do i = 1, min(posList(mo)%size(), 2 * n)
               pvp => posList(mo)%elem(i)
               if (pvp%count < refThreshold) then
                  ii = i - 1
                  exit
               else
                  ii = i
               end if 
               ii = max(ii, 1)
            end do
            pvp => posList(mo)%elem(ii)
            write(iu,'(3f13.6)') pvp%pos(1), pvp%pos(2), pvp%pos(3)
         end do
      end do
      close(iu)


      if (allocated(posList)) then
         do mo = 1, size(posList)
            call posList(mo)%destroy()
         end do      
         deallocate(posList)
      end if

   end subroutine moMax_run

   subroutine moMax_lmo_min(lines, nl, posList)
      ! sample LMOs using ellipsoids with scale factor 'centroidScale' followed by minimizing -phi**2
      character(len=120), intent(in) :: lines(:)
      integer, intent(in) :: nl
      type(pos_VList_t), allocatable, intent(inout) :: posList(:)
      real(r8) :: centroidScale, r(3), distThreshold
      integer :: i, iflag, m, gridSize, sampleSize, verbose, mo
      logical :: found
      class(minimizer_w_sing), pointer :: minimizer_p => null()
      type(fctn_mo)  :: fg
      type(posVal_t) :: pv
      type(posVal_t), pointer :: pvp

      fg = fctn_mo()
      minimizer_p => create_ws_minimizer(lines)

      verbose = 0
      call getinta(lines, nl, "verbose=", verbose, iflag)

      gridSize = 100
      call getinta(lines, nl, 'grid_size=', gridSize, iflag)
      sampleSize = 5
      call getinta(lines, nl, 'attempts_per_mo=', sampleSize, iflag)      
      centroidScale = 1.d0
      call getdbla(lines, nl, 'scale=', centroidScale, iflag)
      distThreshold = 1.d-3
      call getdbla(lines, nl, 'distance_threshold=', distThreshold, iflag)

      call minimizer_p%set_verbose(verbose)
      call minimizer_p%set_verbose_unit(iull)

      call minimizer_p%set_singularities(atoms_getPositionMatrix(atoms), &
                                         scalings=REAL(atoms%Get_atomic_number(), r8))

      ! create ellipsoids around lmo centroids
      call init_lmoSampling(gridSize, centroidScale)

      if (getLMOSize() /= getNOrb()) call error("maximize_mos: # LMOs <> # MOs (norb)")
      if (getLMOSize() /= size(posList)) call error("maximize_mos: # LMOs <> # MOs (posList)")

      do mo = 1, size(posList)

         call posList(mo)%create()
         call fg%setMOIdx(mo)

         do m = 1, sampleSize

            call lmoPosition(mo, r)
            call minimizer_p%reset()
            call minimizer_p%minimize(fg, r)

            if (minimizer_p%is_converged()) then
               found = .false.
               do i = 1, posList(mo)%size()
                  pvp => posList(mo)%elem(i)
                  if ( NORM2(r - pvp%pos) < distThreshold ) then 
                     found = .true.
                     exit
                  end if
               end do
               if (found .and. abs(minimizer_p%value()) > mValueThreshold) then
                  pvp%count = pvp%count + 1
               else if (.not.found .and. abs(minimizer_p%value()) > mValueThreshold) then 
                  pv = posVal_t(r, minimizer_p%value(), 1)
                  call posList(mo)%append(pv)
               end if
            end if
            if (verbose >= 1) then
               write(iul, '(2i6,l5,3g16.5,g16.6,i6)') mo, m, minimizer_p%is_converged(), r, &
                  minimizer_p%value(), minimizer_p%iterations()
            end if

         end do
      end do 

      call destroy_lmoSampling()
      nullify(minimizer_p)
      
   end subroutine moMax_lmo_min

   subroutine moMax_grid(lines, nl, posList)
      ! alternative simple 3d search code for testing
      character(len=120), intent(in) :: lines(:)
      integer, intent(in) :: nl
      type(pos_VList_t), allocatable, intent(inout) :: posList(:) 
      real(r8) :: r(3), value
      real(r8) :: xStart, yStart, zStart, xStep, yStep, zStep, gridStep
      real(r8), allocatable :: vox(:,:,:)
      integer :: ix, iy, iz, mo, verbose, xGridSize, yGridSize, zGridSize, iflag
      integer, parameter :: iu = 99
      type(posVal_t) :: maxP

      verbose = 0
      call getinta(lines, nl, "verbose=", verbose, iflag)

      xGridSize = 100
      call getinta(lines, nl, 'x_grid_size=', xGridSize, iflag)
      yGridSize = 100
      call getinta(lines, nl, 'y_grid_size=', yGridSize, iflag)
      zGridSize = 100
      call getinta(lines, nl, 'z_grid_size=', zGridSize, iflag)

      gridStep = 0.02d0
      call getdbla(lines, nl, 'grid_step=', gridStep, iflag)

      xStep = gridStep / bohr2ang
      yStep = gridStep / bohr2ang
      zStep = gridStep / bohr2ang

      xStart = -(xGridSize/2) * gridStep / bohr2ang
      yStart = -(yGridSize/2) * gridStep / bohr2ang
      zStart = -(zGridSize/2) * gridStep / bohr2ang

      write(iul,'(a,g15.5,a,g15.5,a,i5)') 'xStart =', xStart*bohr2ang, ' xStep =', xStep*bohr2ang, ' xGridSize =', xGridSize
      write(iul,'(a,g15.5,a,g15.5,a,i5)') 'yStart =', yStart*bohr2ang, ' yStep =', yStep*bohr2ang, ' yGridSize =', yGridSize
      write(iul,'(a,g15.5,a,g15.5,a,i5)') 'zStart =', zStart*bohr2ang, ' zStep =', zStep*bohr2ang, ' zGridSize =', zGridSize

      allocate(vox(0:xGridSize, 0:yGridSize, 0:zGridSize))

      if (getNOrb() /= size(posList)) call error("maximize_mos: # MOs <> # MOs (posList)")

      do mo = 1, size(posList)
         call posList(mo)%create()
         vox = 0.d0
         do ix = 0, xGridSize
            r(1) = xStart + ix * xStep
            do iy = 0, yGridSize
               r(2) = yStart + iy * yStep
               do iz = 0, zGridSize
                  r(3) = zStart + iz * zStep
                  call mo_mo(r, mo, value)
                  !!!write(iull, '(i5, 3g15.4, g16.6)') mo, bohr2ang*r, value
                  vox(ix, iy, iz) = -value
               end do
            end do
         end do

         do ix = 1, xGridSize - 1
            do iy = 1, yGridSize - 1
               do iz = 1, zGridSize - 1

                  if ( &
                      vox(ix,iy,iz) < vox(ix-1, iy, iz) .and. &
                      vox(ix,iy,iz) < vox(ix+1, iy, iz) .and. &
                      vox(ix,iy,iz) < vox(ix, iy-1, iz) .and. &
                      vox(ix,iy,iz) < vox(ix, iy+1, iz) .and. &
                      vox(ix,iy,iz) < vox(ix, iy, iz-1) .and. &
                      vox(ix,iy,iz) < vox(ix, iy, iz+1)       &
                      ) then

                     if (abs(vox(ix, iy, iz)) > mValueThreshold) then
                        r = [ xStart + ix * xStep, yStart + iy * yStep, zStart + iz * zStep ]
                        maxP = posVal_t(r, vox(ix, iy, iz), 0)
                        call posList(mo)%append(maxP)
                     end if
                  end if
               end do
            end do
         end do

      end do

      deallocate(vox)
      
   end subroutine moMax_grid

   subroutine moMax_plot(lines, nl)
      ! write MO function values for plotting
      character(len=120), intent(in) :: lines(:)
      integer, intent(in) :: nl
      real(r8) :: r(3), value
      real(r8) :: xStart, yStart, zStart, xStep, yStep, zStep, gridStep
      integer :: io, ix, iy, iz, mo, xGridSize, yGridSize, zGridSize, iflag
      integer, parameter :: iu = 99

      xGridSize = 100
      call getinta(lines, nl, 'x_grid_size=', xGridSize, iflag)
      yGridSize = 100
      call getinta(lines, nl, 'y_grid_size=', yGridSize, iflag)
      zGridSize = 100
      call getinta(lines, nl, 'z_grid_size=', zGridSize, iflag)

      gridStep = 0.02d0
      call getdbla(lines, nl, 'grid_step=', gridStep, iflag)

      xStep = gridStep / bohr2ang
      yStep = gridStep / bohr2ang
      zStep = gridStep / bohr2ang

      xStart = -(xGridSize/2) * gridStep / bohr2ang
      yStart = -(yGridSize/2) * gridStep / bohr2ang
      zStart = -(zGridSize/2) * gridStep / bohr2ang

      open(iu,file=trim(baseName)//'.plt', iostat=io)
      if (io/=0) call error('($plot_mos): opening file '//trim(baseName)//'.plt failed')

      write(iul,'(/a,i5,a/)') 'printing plot data for ', getNOrb(), ' MOs to file '//trim(baseName)//'.plt'

      write(iu,'(i5,2g15.5)') xGridSize, xStart, xStep
      write(iu,'(i5,2g15.5)') yGridSize, yStart, yStep
      write(iu,'(i5,2g15.5)') zGridSize, zStart, zStep
      write(iu,'(i5)') getNOrb()

      do mo = 1, getNOrb()
         write(iu,'(i5)') mo
         do ix = 0, xGridSize
            r(1) = xStart + ix * xStep
            do iy = 0, yGridSize
               r(2) = yStart + iy * yStep
               do iz = 0, zGridSize
                  r(3) = zStart + iz * zStep
                  call mo_mo(r, mo, value)
                  write(iu,'(ES18.8)')  value
               end do
            end do
         end do
      end do

      close(iu)

   end subroutine moMax_plot

   subroutine moMax_plot_plane(lines, nl)
      ! write MO function values for plotting a given plane
      character(len=120), intent(in) :: lines(:)
      integer, intent(in) :: nl
      real(r8) :: r(3), value
      real(r8) :: xStart, yStart, xStep, yStep, gridStep, nucInPlaneThresh
      real(r8) :: ex(3), ey(3), ex0(3), ey0(3), rStart(3), normalOfPlane(3), A(3,3)
      real(r8) :: origin(3), point1(3), point2(3)
      real(r8), allocatable :: B(:, :), extraPoints(:,:)
      integer :: io, idx, ix, iy, mo, xGridSize, yGridSize, iflag, nExtraPoints
      integer :: nRHS, ipiv(getNNuc()), info, countInPlane, ndim, n, n0
      integer, parameter :: iu = 99

      idx = ifinda(lines, nl, 'grid_data=')
      if (idx == 0) call error('$plot_mo_in_plane: grid_data required')

      ! note: construct x direction ex with point1 and y directions with point2
      ! use negative xGridSize for reversed direction
      read(lines(idx+1),*) nExtraPoints, origin
      read(lines(idx+2),*) xGridSize, point1
      read(lines(idx+3),*) yGridSize, point2

      allocate(B(3, getNNuc() + nExtraPoints), extraPoints(3, nExtraPoints))
      ! read further points e.g. ref points to be projected into plot plane
      do n = 1, nExtraPoints
         read(lines(idx+3+n),*) extraPoints(:, n)
      end do

      ex = point1 - origin
      if (xGridSize > 0) ex = - ex
      ey = point2 - origin
      if (yGridSize < 0) ey = - ey

      xGridSize = abs(xGridSize)
      yGridSize = abs(yGridSize)

      ex0 = ex / NORM2(ex)
      ey = ey - dot_product(ex0, ey) * ex0
      ey0 = ey / NORM2(ey)

      write(iul,'(a,i5,3g16.6)') 'ex:', xGridSize, ex0
      write(iul,'(a,i5,3g16.6)') 'ey:', yGridSize, ey0
      write(iul,'(a,3g16.6)') 'orthonormal?', NORM2(ex0), NORM2(ey0), dot_product(ex0, ey0)

      nucInPlaneThresh = 0.01d0
      call getdbla(lines, nl, 'nuc_in_plane_thresh=', nucInPlaneThresh, iflag)

      gridStep = 0.02d0
      call getdbla(lines, nl, 'grid_step=', gridStep, iflag)

      xStep = gridStep / bohr2ang
      yStep = gridStep / bohr2ang
      origin = origin / bohr2ang

      xStart = -(xGridSize/2) * xStep
      yStart = -(yGridSize/2) * yStep
      rStart = origin + xStart * ex0 + yStart * ey0
      write(iul,'(a,3g16.6)') 'rStart:', rStart

      open(iu,file=trim(baseName)//'.plt', iostat=io)
      if (io/=0) call error('($plot_mos): opening file '//trim(baseName)//'.plt failed')

      write(iul,'(/a,i5,a/)') 'printing plot data for ', getNOrb(), ' MOs to file '//trim(baseName)//'.plt'

      ! plane spanned by ex0, ey0 from origin
      ! normal = ex0 x ey0
      ! line thru nuc normal to plane cuts plane: p_nuc + x1*normal = origin + x2*ex0 + x3*ey0
      ! abs(x1) is distance to plane, x2 and x3 are coordinates in plane coord system
      normalOfPlane = cross_product(ex0, ey0)
      A(:, 1) = -normalOfPlane
      A(:, 2) = ex0
      A(:, 3) = ey0
      ndim = 3
      do n = 1, getNNuc()
         B(:, n) = atom_getPosition(atoms(n)) - origin
      end do
      do n = 1, nExtraPoints
         B(:, n + getNNuc()) = extraPoints(:, n) - origin
      end do
      nRHS = getNNuc() + nExtraPoints
      call DGESV(ndim, nRHS, A, ndim, ipiv, B, ndim, info)   ! LAPACK solve A*X = B
      if (info /= 0) write(iul,'(a,i5)') 'moMax_plot_plane: DGESV failed with info=', info 
      countInPlane = 0
      do n = 1, getNNuc()
         if (abs(B(1, n)) < nucInPlaneThresh) countInPlane = countInPlane + 1
      end do
      write(iu,'(i5,a)') countInPlane, ' nuclei in plane:'
      do n = 1, getNNuc()
         write(iul,'(i5,3g15.5)') n, B(:, n) * bohr2ang
         if (abs(B(1, n)) < nucInPlaneThresh) write(iu,'(2g15.5)') B(2, n), B(3, n)
      end do

      countInPlane = 0
      n0 = getNNuc()
      do n = 1, nExtraPoints
         if (abs(B(1, n + n0)) < nucInPlaneThresh) countInPlane = countInPlane + 1
      end do
      write(iu,'(i5,a)') countInPlane, ' extra points in plane:'
      do n = 1, nExtraPoints
         write(iul,'(i5,3g15.5)') n, B(:, n + n0) * bohr2ang
         if (abs(B(1, n + n0)) < nucInPlaneThresh) write(iu,'(2g15.5)') B(2, n + n0), B(3, n + n0)
      end do

      write(iu,'(i5,2g15.5)') xGridSize, xStart, xStep
      write(iu,'(i5,2g15.5)') yGridSize, yStart, yStep
      write(iu,'(i5)') getNOrb()

      do mo = 1, getNOrb()
         write(iu,'(i5)') mo
         do ix = 0, xGridSize
            do iy = 0, yGridSize
               r = rStart + ix * xStep * ex0 + iy * yStep * ey0
               call mo_mo(r, mo, value)
               write(iu,'(ES18.8)')  value
            end do
         end do
      end do

      close(iu)

   end subroutine moMax_plot_plane

end module moMax_m
