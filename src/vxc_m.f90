! Copyright (C) 2019-2020 Vladimir Terzi

module vxc_m
   use kinds_m, only: r8
   use parsing_m, only: getinta, getstra
   use globalUtils_m, only: iul
   use global_m, only: getNElec
   use error_m, only: asserts, assert
   use utils_m, only: printHeader, getFormattedHeader
   use wfData_m, only: getNNuc, atoms
   use rhoMax_m, only: rhoMax_t
   use rwSample_m, only: rwSample, getSampleSize, getFirst, getNext
   use randomWalker_m, only: RandomWalker, pos
   use rhoGrid_m, only: RhoGrid_t

   implicit none

   private
   public :: CalcVxc

   ! TODO: separate initialization of the array of electrons and assigning them to basins from the Vxc calculation
   ! TODO: check which is right: Vxc=Vee-Vcl or Vxc=Vee-(getNElec()-1)/getNElec()*Vcl
   ! TODO: check if VxcAA value needs to be divided by 2

   real(r8), parameter :: Eh = 4.3597447222071e-18_r8  ! https://physics.nist.gov/cgi-bin/cuu/Value?hr
   real(r8), parameter :: NA = 6.02214076e23_r8  ! https://physics.nist.gov/cgi-bin/cuu/Value?na
   integer, parameter :: warnOutUnit = 7  ! TODO: choose better unit for outputing warnings
   character(len = *), parameter :: defaultBasinSearch = 'rho'
   logical :: vxcIsInited = .false., vxcRhoMaxIsInited = .false.
   integer :: smplSize, validSmplSize, rhoMaxValidSmplSize
   real(r8), allocatable :: elecs(:, :)
   integer, allocatable :: basins(:), rhoMaxBasins(:), nSelectedElecsInPoint(:, :)

contains
   subroutine InitVxc(smpl)
      type(rwSample), intent(inout) :: smpl

      real(r8) :: smplPoint(3 * getNElec())
      type(RandomWalker), pointer :: rwp
      integer :: x, i, j

      smplSize = getSampleSize(smpl)
      allocate(elecs(3, smplSize * getNElec()), basins(smplSize * getNElec()), nSelectedElecsInPoint(2, smplSize))
      rwp => getFirst(smpl)
      i = 0  ! iterator for the array of electrons
      do x = 1, smplSize
         call pos(rwp, smplPoint)
         do j = 1, 3 * getNElec(), 3
            i = i + 1
            elecs(:, i) = smplPoint(j : j + 2)
         end do
         rwp => getNext(smpl)
      end do
      vxcIsInited = .true.
   end subroutine InitVxc


   subroutine InitVxcRhoMax(rhoMax)
      type(rhoMax_t), intent(inout) :: rhoMax

      real(r8) :: rMax(3)
      logical :: converged
      integer :: x, i, k, n

      n = 0
      do i = 1, SIZE(elecs, 2)
         call rhoMax%getBasin(elecs(:, i), basins(i), rMax, converged)
         if (basins(i) == 0) then
            n = n + 1
            write(warnOutUnit, '(6(a, f0.3), a)', advance = 'no') 'WARNING: (', elecs(1, i), ',', elecs(2, i), ',', &
               elecs(3, i), ') -> (', rMax(1), ',', rMax(2), ',', rMax(3), ') '
            if (converged) then
               write(warnOutUnit, '(a)') 'was assigned to an undefined basin.'
            else
               write(warnOutUnit, '(a)') 'did not converge.'
            end if
         end if
      end do

      validSmplSize = smplSize
      i = 1
      do x = 1, smplSize
         k = i + getNElec() - 1
         if (ANY(basins(i : k) == 0)) then
            basins(i : k) = 0
            validSmplSize = validSmplSize - 1
         end if
         i = k + 1
      end do

      write(iul, '(2(i0, a))') n, ' electrons out of ', SIZE(elecs, 2), ' did not converge to a defined basin.'
      write(iul, '(2(i0, a))') smplSize - validSmplSize, ' sample points out of ', smplSize, ' were discarded.'
      write(iul, *)

      vxcRhoMaxIsInited = .true.

      ! Saving a copy of the basins array and the number of valid sample points, in case it will be overwritten by
      ! InitVxcRhoGrid.
      allocate(rhoMaxBasins(SIZE(basins)))
      rhoMaxBasins = basins
      rhoMaxValidSmplSize = validSmplSize
   end subroutine InitVxcRhoMax


   subroutine InitVxcRhoGrid(rhoGrid)
      type(RhoGrid_t), intent(in) :: rhoGrid

      real(r8) :: rMax(3)
      logical :: inBounds
      integer :: x, i, k, n

      n = 0
      do i = 1, SIZE(elecs, 2)
         call rhoGrid%GetBasin(elecs(:, i), basins(i), rMax, inBounds)
         if (basins(i) == 0) then
            n = n + 1
             write(warnOutUnit, '(6(a, f0.3), a)', advance = 'no') 'WARNING: (', elecs(1, i), ',', elecs(2, i), ',', &
                elecs(3, i), ') -> (', rMax(1), ',', rMax(2), ',', rMax(3), ') '
             if (.not. inBounds) then
                write(warnOutUnit, '(a)') 'is out of grid bounds.'
             else
                write(warnOutUnit, '(a)') 'did not converge.'
             end if
         end if
      end do

      validSmplSize = smplSize
      i = 1
      do x = 1, smplSize
         k = i + getNElec() - 1
         if (ANY(basins(i : k) == 0)) then
            basins(i : k) = 0
            validSmplSize = validSmplSize - 1
         end if
         i = k + 1
      end do

      write(iul, '(2(i0, a))') n, ' electrons out of ', SIZE(elecs, 2), ' did not converge to a defined basin.'
      write(iul, '(2(i0, a))') smplSize - validSmplSize, ' sample points out of ', smplSize, ' were discarded.'
      write(iul, *)
   end subroutine InitVxcRhoGrid


   subroutine InitElecsAA(a, elecsA)
      integer, intent(in) :: a
      real(r8), dimension(:, :), allocatable, intent(out) :: elecsA

      integer :: x, i, j, iA

      ! counting electrons in the selected basin for each sample point
      i = 0  ! iterator for the array of electrons
      do x = 1, smplSize
         iA = 0
         do j = 1, getNElec()
            i = i + 1
            if (basins(i) == a) then
               iA = iA + 1
            end if
         end do
         nSelectedElecsInPoint(1, x) = iA
      end do

      ! saving selected electrons in the array
      allocate(elecsA(3, SUM(nSelectedElecsInPoint(1, :))))
      iA = 0
      do i = 1, SIZE(basins)
         if (basins(i) == a) then
            iA = iA + 1
            elecsA(:, iA) = elecs(:, i)
         end if
      end do
   end subroutine InitElecsAA


   subroutine InitElecsAB(a, b, elecsA, elecsB)
      integer, intent(in) :: a, b
      real(r8), dimension(:, :), allocatable, intent(out) :: elecsA, elecsB

      integer :: x, i, j, iA, iB

      ! counting electrons in the selected basins for each sample point
      i = 0  ! iterator for the array of electrons
      do x = 1, smplSize
         iA = 0
         iB = 0
         do j = 1, getNElec()
            i = i + 1
            if (basins(i) == a) then
               iA = iA + 1
            else if (basins(i) == b) then
               iB = iB + 1
            end if
         end do
         nSelectedElecsInPoint(1, x) = iA
         nSelectedElecsInPoint(2, x) = iB
      end do

      ! saving selected electrons in the corresponding arrays
      allocate(elecsA(3, SUM(nSelectedElecsInPoint(1, :))), elecsB(3, SUM(nSelectedElecsInPoint(2, :))))
      iA = 0
      iB = 0
      do i = 1, SIZE(basins)
         if (basins(i) == a) then
            iA = iA + 1
            elecsA(:, iA) = elecs(:, i)
         else if (basins(i) == b) then
            iB = iB + 1
            elecsB(:, iB) = elecs(:, i)
         end if
      end do
   end subroutine InitElecsAB


   subroutine CalcVeeAA(elecsA, vee, intee)
      real(r8), dimension(:, :), intent(in) :: elecsA
      real(r8), intent(out) :: vee, intee

      real(r8) :: r(3)
      integer :: x, i, k, n, j1, j2

      ! The selected electrons are saved in the same order in the array as the sample points. For each sample point
      ! initial (i) and final (k) electron indices are calculated using the saved numbers. The electrons
      ! between these indices are read from the arrays consequentially to build all possible electron pairs.
      vee = 0._r8
      n = 0
      i = 1
      do x = 1, smplSize
         k = i + nSelectedElecsInPoint(1, x) - 1
         do j1 = i, k - 1
            r = elecsA(:, j1)
            do j2 = j1 + 1, k
               vee = vee + 1 / NORM2(r - elecsA(:, j2))
               n = n + 1
            end do
         end do
         i = k + 1
      end do
      vee = vee / validSmplSize
      intee = 1._r8 * n / validSmplSize
   end subroutine CalcVeeAA


   subroutine CalcVeeAB(elecsA, elecsB, vee, intee)
      real(r8), dimension(:, :), intent(in) :: elecsA, elecsB
      real(r8), intent(out) :: vee, intee

      real(r8) :: r(3)
      integer :: x, n, iA, iB, jA, jB, kA, kB

      ! The selected electrons are saved in the same order in the corresponding arrays as the sample points. For each
      ! sample point initial (iA / iB) and final (kA / kB) electron indices are calculated using the saved numbers. The
      ! electrons between these indices are read from the arrays consequentially to build all possible electron pairs.
      vee = 0._r8
      n = 0
      iA = 1
      iB = 1
      do x = 1, smplSize
         kA = iA + nSelectedElecsInPoint(1, x) - 1
         kB = iB + nSelectedElecsInPoint(2, x) - 1
         do jA = iA, kA
            r = elecsA(:, jA)
            do jB = iB, kB
               vee = vee + 1 / NORM2(r - elecsB(:, jB))
               n = n + 1
            end do
         end do
         iA = kA + 1
         iB = kB + 1
      end do
      vee = vee / validSmplSize
      intee = 1._r8 * n / validSmplSize
   end subroutine CalcVeeAB


   subroutine CalcVclAA(elecsA, vcl, intcl)
      real(r8), dimension(:, :), intent(inout) :: elecsA
      real(r8), intent(out) :: vcl, intcl

      integer :: i, j
      real(r8) :: r(3)

      vcl = 0._r8
      do i = 1, SIZE(elecsA, 2) - 1
         r = elecsA(:, i)
         do j = i + 1, SIZE(elecsA, 2)
            vcl = vcl + 1 / NORM2(r - elecsA(:, j))
         end do
      end do
      ! 1._r8 in denominator is needed to avoid integer(4) overflow
      vcl = vcl * getNElec() / (validSmplSize * (1._r8 * validSmplSize * getNElec() - 1))
      intcl = .5_r8 * SIZE(elecsA, 2) * (SIZE(elecsA, 2) - 1) * getNElec() &
         / (validSmplSize * (1._r8 * validSmplSize * getNElec() - 1))
   end subroutine CalcVclAA


   subroutine CalcVclAB(elecsA, elecsB, vcl, intcl)
      real(r8), dimension(:, :), intent(inout) :: elecsA, elecsB
      real(r8), intent(out) :: vcl, intcl

      integer :: iA, iB
      real(r8) :: r(3)
      
      vcl = 0._r8
      do iA = 1, SIZE(elecsA, 2)
         r = elecsA(:, iA)
         do iB = 1, SIZE(elecsB, 2)
            vcl = vcl + 1 / NORM2(r - elecsB(:, iB))
         end do
      end do
      ! 1._r8 in denominator is needed to avoid integer(4) overflow
      vcl = vcl * getNElec() / (validSmplSize * (1._r8 * validSmplSize * getNElec() - 1))
      intcl = 1._r8 * SIZE(elecsA, 2) * SIZE(elecsB, 2) * getNElec() &
         / (validSmplSize * (1._r8 * validSmplSize * getNElec() - 1))
   end subroutine CalcVclAB


   subroutine CalcVxcAA(a)
      integer, intent(in) :: a

      real(r8), dimension(:, :), allocatable :: elecsA
      real(r8) :: vee, vcl, vxc
      real(r8) :: intee, intcl, intxc

      write(iul, '(2a, i0)') 'Vxc for the atom ', TRIM(atoms(a)%elem), a
      write(iul, *)

      call InitElecsAA(a, elecsA)

      call CalcVeeAA(elecsA, vee, intee)
      write(iul, '(2(a, f0.16), a)') 'Vee = ', vee, ' E_h = ', vee * Eh * NA / 1e3_r8, ' kJ/mol'
      write(iul, '(a, f0.16)') 'int rho2(1,2)/2 = ', intee
      write(iul, *)

      call CalcVclAA(elecsA, vcl, intcl)
      write(iul, '(2(a, f0.16), a)') 'Vcl = ', vcl, ' E_h = ', vcl * Eh * NA / 1e3_r8, ' kJ/mol'
      write(iul, '(a, f0.16)') 'int rho(1)*rho(2)/2 = ', intcl
      write(iul, *)

      deallocate(elecsA)

      vxc = vee - vcl
      intxc = intee - intcl
      write(iul, '(2(a, f0.16), a)') 'Vxc = ', vxc, ' E_h = ', vxc * Eh * NA / 1e3_r8, ' kJ/mol'
      write(iul, '(a, f0.16)') 'int rhoxc(1,2)/2 = ', intxc
      write(iul, *)
   end subroutine CalcVxcAA


   subroutine CalcVxcAB(a, b)
      integer, intent(in) :: a, b

      real(r8), dimension(:, :), allocatable :: elecsA, elecsB
      real(r8) :: vee, vcl, vxc
      real(r8) :: intee, intcl, intxc

      write(iul, '(2(2a, i0))') 'Vxc for the bond between ', TRIM(atoms(a)%elem), a, ' and ', TRIM(atoms(b)%elem), b
      write(iul, *)

      call InitElecsAB(a, b, elecsA, elecsB)

      call CalcVeeAB(elecsA, elecsB, vee, intee)
      write(iul, '(2(a, f0.16), a)') 'Vee = ', vee, ' E_h = ', vee * Eh * NA / 1e3_r8, ' kJ/mol'
      write(iul, '(a, f0.16)') 'int rho2(1,2) = ', intee
      write(iul, *)

      call CalcVclAB(elecsA, elecsB, vcl, intcl)
      write(iul, '(2(a, f0.16), a)') 'Vcl = ', vcl, ' E_h = ', vcl * Eh * NA / 1e3_r8, ' kJ/mol'
      write(iul, '(a, f0.16)') 'int rho(1)*rho(2) = ', intcl
      write(iul, *)

      deallocate(elecsA, elecsB)

      vxc = vee - vcl
      intxc = intee - intcl
      write(iul, '(2(a, f0.16), a)') 'Vxc = ', vxc, ' E_h = ', vxc * Eh * NA / 1e3_r8, ' kJ/mol'
      write(iul, '(a, f0.16)') 'int rhoxc(1,2) = ', intxc
      write(iul, *)
   end subroutine CalcVxcAB


   subroutine CalcVxc(lines, nl, smpl, rhoMax, rhoGrid)
      integer, intent(in) :: nl
      character(len = 120), intent(in) :: lines(nl)
      type(rwSample), intent(inout) :: smpl
      type(rhoMax_t), intent(inout) :: rhoMax
      type(RhoGrid_t), intent(in) :: rhoGrid

      integer :: a, b
      character(len = 80) :: basinSearch
      integer :: iflagA, iflagB, iflag
      real(r8) :: start, finish

      call getinta(lines, nl, 'A=', a, iflagA)
      call getinta(lines, nl, 'B=', b, iflagB)
      if (iflagA /= 0 .and. iflagB /= 0) then
         if (asserts) call assert(.false., "CalcVxc: no atom indices are specified.")
      else if (iflagA == 0 .and. iflagB /= 0) then
         b = a
      else if (iflagA /= 0 .and. iflagB == 0) then
         a = b
      end if
      if (asserts) call assert(a >= 1 .and. a <= getNNuc() .and. b >= 1 .and. b <= getNNuc(), &
         "CalcVxc: the atom indeces are invalid.")

      call getstra(lines, nl, 'basin_search=', basinSearch, iflag)
      if (iflag /= 0) basinSearch = defaultBasinSearch

      if (.not. vxcIsInited) then
         call InitVxc(smpl)
      end if
      if (basinSearch == 'rho') then
         if (.not. vxcRhoMaxIsInited) then
            call printHeader(iul, getFormattedHeader('vxc', 'assigning electrons to basins'), '-')
            write(iul, *)
            call CPU_TIME(start)

            if (asserts) call assert(rhoMax%isInitialized(), "CalcVxc: rhoMax is not initialized.")
            call InitVxcRhoMax(rhoMax)

            call CPU_TIME(finish)
            write(iul, '(a, f0.3, a)') 'cpu time for initializing Vxc: ', finish - start, ' s'
            write(iul, *)
         else
            basins = rhoMaxBasins
            validSmplSize = rhoMaxValidSmplSize
         end if
      else if (basinSearch == 'rho_grid') then
         call printHeader(iul, getFormattedHeader('vxc', 'assigning electrons to basins'), '-')
         write(iul, *)
         call CPU_TIME(start)

         if (asserts) call assert(rhoGrid%IsInitialized(), 'CalcVxc: electronic density grid is not initialized.')
         if (asserts) call assert(rhoGrid%AtomInGrid(a) .and. rhoGrid%AtomInGrid(b), 'CalcVxc: the specified atoms ' &
            // 'are not in the defined electronic density grid.')
         call InitVxcRhoGrid(rhoGrid)

         call CPU_TIME(finish)
         write(iul, '(a, f0.3, a)') 'cpu time for initializing Vxc: ', finish - start, ' s'
         write(iul, *)
      else
         if (asserts) call assert(.false., 'CalcVxc: the specified basin search method is invalid.')
      end if
      if (asserts) call assert(validSmplSize > 0, 'CalcVxc: no valid sample points were found.')

      call CPU_TIME(start)
      
      if (a == b) then
         call CalcVxcAA(a)
      else
         call CalcVxcAB(a, b)
      end if

      call CPU_TIME(finish)
      write(iul, '(a, f0.3, a)') 'cpu time for calculating Vxc: ', finish - start, ' s'
   end subroutine CalcVxc
end module vxc_m
