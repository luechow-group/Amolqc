! Copyright (C) 2011-2013 Alexander Sturm
! Copyright (C) 2012-2013, 2015 Arne Luechow
! Copyright (C) 2013, 2015-2016 Kaveh Haghighi Mood
! Copyright (C) 2014, 2016, 2018 Christoph Schulte
!
! SPDX-License-Identifier: GPL-3.0-or-later

module jastrowIC_m

  use kinds_m, only: r8
  use wfData_m
  use jastrowParamData_m
  use aosData_m
  use rdataUpdate_m
  use jastrowAniso_m

  implicit none


  public :: jasinput_ic, jasoutput_ic, jasoutput_ic_new, jas_shortoutput_ic, jasChangeType_ic, &
            jasicall, jasicp, getVectorLenIC, getVectorIC, putVectorIC, jas_diffeecusp_ic, &
            jasicInit, jasicUpdate, jasicInitWithUk, jasicUpdateWithUk

  ! max index of parameter
  integer :: umax, xmax, fmax
  ! number of parameters per term
  integer :: unum=0, xnum=0, fnum=0, gnum=0
  ! number of e-n/e-e-n-terms per core
  integer :: xpnum, fpnum
  ! total number of parameters
  integer :: numParams
  ! coefficients for power expansions
  real(r8) :: alpha(kmax) = 0d0
  real(r8) :: beta(kmax, amax) = 0d0
  real(r8) :: gamma(kmax, amax) = 0d0
  real(r8), allocatable :: gl(:)
  ! last nucleus with correlation
  integer :: nclast, ncdiff
  ! whether electron nucleus should be enforced; whether different cusp condition for same and opposite spin ee
  logical :: nucCusp, diffeecusp
  ! use generic or hardcoded terms
  logical :: useGeneric = .false.
  logical :: useAOJasTerms = .false.

  !anisotropic Jastrow Object
  type(JasAnisoType) :: JasAniso
  !count of one and two body anisotropic terms
  integer :: anisoJ1_count, anisoJ2_count

  ! terms used in calculations
  ! reusable arrays allocated
  logical :: allocatedTerms = .false.
  ! powers of distances
  real(r8), allocatable :: eePowers(:, :, :), enPowers(:, :, :)
  ! derivatives of distances
  real(r8), allocatable :: rijDeriv(:, :,  :), raiDeriv(:, :, :)
  ! (dr/dx)^2 + (dr/dy)^2 + (dr/dz)^2
  real(r8), allocatable :: rijSquare(:, :), raiSquare(:, :)
  ! laplacian
  real(r8), allocatable :: rijLapl(:, :), raiLapl(:, :)

  ! scaled distance type
  ! SM distance function 1/(1+a*r)
  integer, parameter :: DIST_SM = 1
  ! double exponential function 1 - exp(-a*r)
  integer, parameter :: DIST_DOUBLEEXP = 2
  ! power of fractions r/(a+r^b)
  ! by Lopez Rios, Seth, Drummond, Needs DOI: 10.1103/PhysRevE.86.036703
  integer, parameter :: DIST_NEEDS = 3

  integer :: distType

  ! sm/doubleexp params
  real(r8) :: scaleEE = 1 ! scale factor for ee distance
  real(r8) :: scaleEN(amax) = 1 ! scale factor for en distance

  ! additional params for DIST_NEEDS
  real(r8) :: powerEE = 1
  real(r8) :: powerEN(amax) = 1

  !! different optmodes
  !! only linear parameters, default
  !integer, parameter :: OPT_LIN = 1
  !! only linear parameters, numerical
  !integer, parameter :: OPT_NUM_LIN = 2
  !! only non-linear parameter (always numerical)
  !integer, parameter :: OPT_NONLIN = 3
  !! all parameters (linear analytical, nonlinear numerical)
  !integer, parameter :: OPT_ALL = 4
  !! all parameters (numerical)
  !integer, parameter :: OPT_NUM_ALL = 5

  ! current optMode, used to identify which parameter derivatives need
  ! to be returned in jastrow calculation
  ! XXX should be removed, pass optMode to jastrow calc?
  integer :: curOptMode

contains

!=======================================================================


subroutine jasinput_ic(lines,nl)
!------------------------------!
  ! read Jastrow related input from 'lines'
  character(len=*), intent(in) :: lines(:)! lines array
  integer, intent(in)          :: nl      ! actual # of lines
  integer :: a, i, m, stat, idx, nWords
  integer :: offset
  character :: t
  character(len=10) :: word(20)

  diffeecusp=.FALSE.
  idx = 2


  read(lines(idx), *) nclast
  if(nclast == 0) nclast = atoms_getNCLast(atoms)
  call assert(nclast <= atoms_getNCLast(atoms), "nclast is too large!")
  ncdiff = atoms(nclast)%sa
  idx = idx+1
  read(lines(idx),"(2L2)", iostat=stat) nucCusp, diffeecusp
  idx = idx+1
  call tokenize(lines(idx),word,nWords)
  read(lines(idx), *) umax, xmax, fmax   ! ee, en max powert, een max degree
  select case (nWords)
  case (3)
    ! only ic ee, en, een
    ! already read
  case (6:)
    ! additional terms; currently only anisotropic AO terms for en, een or eenn
    call JasAniso%init(lines,idx,useAOJasTerms)
  case default
    call abortp("ic jastrow input: wrong format in 3rd line")
  end select

  call setNumParams()

  idx = idx + 1
  read(lines(idx), "(A)") t
  idx = idx + 1
  if(t /= 'x') then
    idx = idx - 1
    if(unum > 0) then
      read(lines(idx), *) (alpha(m), m = 2, umax)
      idx = idx + 1
    endif

    if(xnum > 0) then
      do i = 1, ncdiff
        read(lines(idx), *) (beta(m, i), m = 2, xmax)
        idx = idx + 1
      enddo
    endif

    if(fnum > 0) then
      do m = 1, fpnum
        read(lines(idx), *) (gamma(m, i), i = 1, ncdiff)
        idx = idx + 1
      enddo
    endif

    if (useAOJasTerms) then
      call JasAniso%read_coefs(lines,idx)
    endif

  else
    alpha = 0d0
    beta = 0d0
    gamma = 0d0
    if (useAOJasTerms) then
      call JasAniso%zero_coefs()
    endif

  endif

  ! for old input format (no dist type specified), default to SM type distance
  ! unless jastype deXYZ is used
  distType = DIST_SM
  if(jastype(1:2) == 'de') distType = DIST_DOUBLEEXP

  scaleEE = 1d0
  scaleEN(1:ncdiff) = 1d0

  read(lines(idx), "(I4)", iostat=stat) distType
  if(stat == 0) then
    select case(distType)
    case(DIST_SM, DIST_DOUBLEEXP)
      offset = 1
      a=1
      read(lines(idx+offset), *, iostat=stat) scaleEE, scaleEN(a)
      if (stat /= 0) call abortp("jastrowIC_m: Error while reading distance parameters")
      offset = offset +1
      do
        if (ncdiff-a == 0 ) then
          !all parameters read
          exit
        else if (ncdiff-a == 1) then
          read(lines(idx+offset), *, iostat=stat) scaleEN(a+1)
          if(stat /= 0) then
            ! use same parameter for all nuclei types if read fails
            scaleEN(2:ncdiff) = scaleEN(1)
            exit
          endif
          a = a+1
          offset = offset+1
        else if (ncdiff-a >= 2) then
          read(lines(idx+offset), *, iostat=stat) scaleEN(a+1),scaleEN(a+2)
          if (stat /= 0) then
            ! use same parameter for all nuclei types if read fails
            scaleEN(2:ncdiff) = scaleEN(1)
            exit
          endif
          a = a+2
          offset = offset+1
        else
          call abortp("jastrowIC_m: Error while reading distance parameters")
        endif
      enddo
    case(DIST_NEEDS)
      !read scaleEE,scaleEN first
      offset = 1
      a=1
      read(lines(idx+offset), *, iostat=stat) scaleEE, scaleEN(a)
      if (stat /= 0) call abortp("jastrowIC_m: Error while reading Needs-distance parameters")
      offset = offset +1
      do
        if (ncdiff-a == 0 ) then
          !all parameters read
          exit
        else if (ncdiff-a == 1) then
          read(lines(idx+offset), *, iostat=stat) scaleEN(a+1)
          if(stat /= 0) then
            ! use same parameter for all nuclei types if read fails
            scaleEN(2:ncdiff) = scaleEN(1)
            exit
          endif
          a = a+1
          offset = offset+1
        else if (ncdiff-a >= 2) then
          read(lines(idx+offset), *, iostat=stat) scaleEN(a+1),scaleEN(a+2)
          if (stat /= 0) then
            ! use same parameter for all nuclei types if read fails
            scaleEN(2:ncdiff) = scaleEN(1)
            exit
          endif
          a = a+2
          offset = offset+1
        else
          call abortp("jastrowIC_m: Error while reading Needs-distance parameters")
        endif
      enddo
      !then read powerEE
      a=0
      do
        if (ncdiff-a == 0 ) then
          !all parameters read
          exit
        else if (ncdiff-a == 1) then
          read(lines(idx+offset), *, iostat=stat) powerEN(a+1)
          if(stat /= 0) then
            ! use same parameter for all nuclei types if read fails
            powerEN(2:ncdiff) = powerEN(1)
            exit
          endif
          a = a+1
          offset = offset+1
        else if (ncdiff-a >= 2) then
          read(lines(idx+offset), *, iostat=stat) powerEN(a+1),powerEN(a+2)
          if (stat /= 0) then
            ! use same parameter for all nuclei types if read fails
            powerEN(2:ncdiff) = powerEN(1)
            exit
          endif
          a = a+2
          offset = offset+1
        else
          call abortp("jastrowIC_m: Error while reading Needs-distance parameters")
        endif
      enddo
    case default
      call abortp("Unknown scaled distance type in jastrow input")
    end select
    idx = idx+offset
  endif

  if(jastype(1:3) == 'icg' .or. jastype(1:3) == 'gmg') then
    useGeneric = .true.
  elseif(fmax > 6) then
    useGeneric = .true.
    if (logmode>=2) write(iul, *) "Using generic jastrow"
  else
    useGeneric = .false.
  endif

  ! normalize saved jastype value to properly handle ic/icg/de/...
  jastype = "ic"

end subroutine jasinput_ic


subroutine jas_addAnisoTerms_ic(lines,nl)
!---------------------------------------!
  ! add anisotropic AO terms to current Jastrow: read input from 'lines'
  character(len=*), intent(in) :: lines(:)! lines array
  integer, intent(in)          :: nl      ! actual # of lines

  call JasAniso%add(lines,nl,useAOJasTerms)
  call setNumParams()
end subroutine jas_addAnisoTerms_ic

!=======================================================================

subroutine jasoutput_ic()
!---------------------!
! write Jastrow related terms in human readable format
  integer :: a, i

  if (.not. MASTER) return

  write(iul, "(/A)") "Jastrow part:"
  if(umax > 0) then
    write(iul, "(/A,I4)") "No. of el-el Jastrow parameters = ", unum
    do i = 2, umax
      write(iul, "(F22.15)", advance="no") alpha(i)
    enddo
    write(iul, *)
  endif

  write(iul, *)

  if(xmax > 0) then
    write(iul, "(/A,I4)") "No. of el-nuc Jastrow parameters = ", xnum
    do i = 2, xmax
      do a = 1, ncdiff
        write(iul, "(F22.15)", advance="no") beta(i, a)
      enddo
      write(iul, *)
    enddo
  endif

  if (diffeecusp) then
      write(iul,"(a)") "Different ee cusp condition for spin like electrons used"
  endif

  write(iul, *)

  if(fmax > 0) then
    write(iul, "(/A,I4)") "No. of el-el-nuc Jastrow terms = ", fnum
    do i = 1, fpnum
      do a = 1, ncdiff
        write(iul, "(F22.15)", advance="no") gamma(i, a)
      enddo
      write(iul, *)
    enddo
  endif

  write(iul, *)
  if (useAOJasTerms) call JasAniso%output()

  write(iul, "(/A,I4)") "Jastrow distance type: ", distType

  select case(distType)
  case(DIST_SM, DIST_DOUBLEEXP)
    write(iul, "(2ES15.7)") scaleEE, (scaleEN(a),a=1,ncdiff)
  case(DIST_NEEDS)
    write(iul, "(2ES15.7)") scaleEE, (scaleEN(a),a=1,ncdiff)
    write(iul, "(2ES15.7)") powerEE, (powerEN(a),a=1,ncdiff)
  case default
    call abortp("Unknown Jastrow distance type in jasoutput")
  end select

  write(iul, *)

end subroutine jasoutput_ic


subroutine jasoutput_ic_new(iu)
!-----------------------------!
! writes jastrow data to log unit iu in input format
! used for writeWF
  integer, intent(in) :: iu
  integer :: a, i, m
  character(len=120) :: tline

  write(iu, "(I4)") nclast
  write(iu, "(2L2)") nucCusp, diffeecusp

  write(tline, "(4I4)") umax, xmax, fmax
  if (useAOJasTerms) call JasAniso%output_header_new(tline)
  write(iu, "(a)") trim(tline)
  if (useAOJasTerms) call JasAniso%output_header_new_lines(iu)

  if(unum > 0) then
    do m = 2, umax
      write(iu, "(ES15.7)", advance="no") alpha(m)
    enddo
    write(iu, *)
  endif

  if(xnum > 0) then
    do i = 1, ncdiff
      do m = 2, xmax
        write(iu, "(ES15.7)", advance="no") beta(m, i)
      enddo
      write(iu, *)
    enddo
  endif

  if(fnum > 0) then
    do m = 1, fpnum
      do i = 1, ncdiff
        write(iu, "(ES15.7)", advance="no") gamma(m, i)
      enddo
      write(iu, *)
    enddo
  endif

  if (useAOJasTerms) call JasAniso%output_param_new(iu)

  write(iu, "(I4)") distType

  select case(distType)
  case(DIST_SM, DIST_DOUBLEEXP)
    write(iu, "(2ES15.7)") scaleEE, (scaleEN(a),a=1,ncdiff)
  case(DIST_NEEDS)
    write(iu, "(2ES15.7)") scaleEE, (scaleEN(a),a=1,ncdiff)
    write(iu, "(2ES15.7)") powerEE, (powerEN(a),a=1,ncdiff)
  case default
    call abortp("Unknown Jastrow distance type in jasoutput")
  end select

end subroutine jasoutput_ic_new

!==============================================================

subroutine jas_shortoutput_ic(iu)
!-------------------------------!
  ! writes jastrow data to log unit iu in input format
  integer, intent(in) :: iu

  select case(distType)
  case (DIST_SM)
    write(iu,'(a)') '  Schmidt-Moskowitz radial function'
  case (DIST_DOUBLEEXP)
    write(iu,'(a)') '  Double exponential radial function'
  case (DIST_NEEDS)
    write(iu,'(a)') '  Drummond-Needs radial function'
  end select
  write(iu,'(3(i4,a))') unum,' ee terms, ',xnum,' en terms, ',fnum,' een terms'
  if (useAOJasTerms) then
     call JasAniso%shortoutput(iu)
  end if
  if (diffeecusp) then
    write(iu,'(a)')"with different cusp term for spin like and spin unlike e-"
  endif
end subroutine jas_shortoutput_ic

!=======================================================================

subroutine jasChangeType_ic(jt)
!-----------------------------!
  ! reset Jastrow parameters (to start e.g. the optimization)
  character(len=*), intent(in) :: jt

  call deallocateTerms()

  if(jastype(1:2) /= jt(1:2)) then
    if ( .not.( (jastype(1:2)== 'ic' .and. distType==DIST_DOUBLEEXP) .and. jt(1:2)=='de') ) then
      ! don't preserve old parameters
      diffeecusp = .false.
      nucCusp = .false.
      nclast = atoms_getNCLast(atoms)
      ncdiff = atoms_getNSCenter(atoms)
      alpha = 0d0
      beta = 0d0
      gamma = 0d0

      ! default distance type is SM
      distType = DIST_SM
      scaleEE = 1d0
      scaleEN(1:ncdiff) = 1d0
    endif
  endif

  if(jt(1:2) == 'de') then
    distType = DIST_DOUBLEEXP
  endif

  if(jt(1:3) == 'icg') then
    useGeneric = .true.
    read(jt, "(3X,I1,I1,I1)") umax, xmax, fmax
  else
    read(jt, "(2X,I1,I1,I1)") umax, xmax, fmax
    if(fmax > 6) then
      useGeneric = .true.
      if(MASTER .and. logmode >= 2) then
        write(iul, *) "Using generic jastrow"
      endif
    else
      useGeneric = .false.
    endif
  endif

  jastype = 'ic'

  if (.not. (useAOJasTerms)) then
    !deallocate anisotropic terms if needed
     call JasAniso%destroy()
  endif
  call setNumParams()

end subroutine jasChangeType_ic

!=======================================================================

subroutine jasicp(init, ie, rai, rij, ju)
!---------------------------------------!
! calculates U without derivatives
  logical, intent(in) ::  init                   ! .true. for initialization
  integer, intent(in) ::  ie                     ! electron with new position
  real(r8), intent(in)  ::  rai(:, :), rij(:, :)   ! current distances
  real(r8), intent(out) ::  ju                     ! returns U rather than exp(U)

  ! dummy variables
  real(r8) :: x(ne), y(ne), z(ne), jud(3*ne), julapl, julapli(ne)

  !call assert(ie == 0 .and. .not. init, "jasicp only implemented for all electron move")
  call jasicall(x, y, z, rai, rij, "none", ju, jud, julapl, julapli,withDerivs= 0)
end subroutine jasicp

!=======================================================================
!=======================================================================

subroutine jasicpWithUk(init, ie, rai, rij, ju, uuk)
!---------------------------------------!
! calculates U without derivatives
  logical, intent(in) ::  init                   ! .true. for initialization
  integer, intent(in) ::  ie                     ! electron with new position
  real(r8), intent(in)  ::  rai(:, :), rij(:, :)   ! current distances
  real(r8), intent(out) ::  ju                     ! returns U rather than exp(U)
  real(r8), intent(out) ::  uuk(:)
  ! dummy variables
  real(r8) :: jud(3*ne), julapl, julapli(ne)
  real(r8) ::x(ne), y(ne), z(ne)
  !call assert(ie == 0 .and. .not. init, "jasicp only implemented for all electron move")
  call jasicall(x, y, z, rai, rij, "jastrow", ju, jud, julapl, julapli, uuk=uuk, withDerivs = 2)
end subroutine jasicpWithUk

!=======================================================================


subroutine getVectorLenIC(optmode,npJ1,npJ2,npJnl)
  integer, intent(in) :: optMode ! optimization mode
  integer, intent(inout) :: npJ1     ! one-electron linear
  integer, intent(inout) :: npJ2     ! two-electron linear
  integer, intent(inout) :: npJnl    ! nonlinear
  integer :: aniso_J1, aniso_J2

  aniso_J1 = 0
  aniso_J2 = 0
  if (useAOJasTerms) call JasAniso%getParamCount(optMode,aniso_J1,aniso_J2)

  npJ1 = xnum + aniso_J1
  npJ2 = unum + fnum + aniso_J2

  select case(optMode)
  case (OPT_LIN, OPT_NUM_LIN)
    npJ1 = xnum + aniso_J1
    npJ2 = unum + fnum + aniso_J2
    npJnl = 0
  case (OPT_NONLIN)
    npJ1 = 0
    npJ2 = 0
    select case(distType)
    case (DIST_SM, DIST_DOUBLEEXP)
      npJnl = 1 + ncdiff
    case (DIST_NEEDS)
      npJnl = 2 + 2*ncdiff
    case default
      call abortp("getVectorLenIC: distType not implemented")
    end select
  case (OPT_ALL, OPT_NUM_ALL)
    npJ1 = xnum + aniso_J1
    npJ2 = unum + fnum + aniso_J2
    select case(distType)
    case (DIST_SM, DIST_DOUBLEEXP)
      npJnl = 1 + ncdiff
    case (DIST_NEEDS)
      npJnl = 2 + 2*ncdiff
    case default
      call abortp("getVectorLenIC: distType not implemented")
    end select
  case default
    call abortp("getVectorLenIC: optMode not implemented")
  end select

  curOptMode = optMode
end subroutine getVectorLenIC


!=======================================================================

recursive subroutine getVectorIC(optMode, p, offset)
  integer, intent(in)   ::  optMode       ! optimization mode, 1 for linear params,
  real(r8), intent(inout) ::  p(:)          ! parameter vector
  integer, optional     ::  offset
  integer               ::  start, aostart, aniso_count
  real(r8), allocatable :: aniso_p(:)

  start = 1
  aostart = 1
  if(present(offset)) start = offset

  select case(optMode)
  case(OPT_LIN, OPT_NUM_LIN) ! linear parameters
    p(start:unum + start - 1) = alpha(2:umax)
    start = start + unum
    p(start:start + xnum - 1) = reshape(beta(2:xmax, 1:ncdiff), [xnum])
    start = start + xnum
    p(start:start + fnum - 1) = reshape(gamma(1:fpnum, 1:ncdiff), [fnum])
    start = start + fnum
    if (useAOJasTerms) then
      call JasAniso%getParamVector(optMode,aniso_p,aniso_count)
      p(start:start+aniso_count-1) = aniso_p(1:aniso_count)
      start = start + aniso_count
    endif
  case(OPT_NONLIN)
    select case(distType)
    case(DIST_SM, DIST_DOUBLEEXP)
      p(start) = scaleEE
      p(start + 1:start + ncdiff) = scaleEN(1:ncdiff)
      start = start + ncdiff
    case(DIST_NEEDS)
      p(start) = scaleEE
      p(start + 1) = powerEE
      p(start + 2:start + 1 + ncdiff) = scaleEN(1:ncdiff)
      p(start + 2 + ncdiff:start + 1 + 2*ncdiff) = powerEN(1:ncdiff)
      start = start + 1 + 2 * ncdiff
    case default
      call abortp("getVectorIC: distType not implemented")
    end select

  case(OPT_ALL, OPT_NUM_ALL)
    call getVectorIC(OPT_LIN, p)
    call getVectorIC(OPT_NONLIN, p, numParams+1)
  case default
    call abortp("getVectorIC: optMode not implemented")
  end select

end subroutine getVectorIC

!=======================================================================

recursive subroutine putVectorIC(optMode, p)
!------------------------------------------!
  integer, intent(in) ::  optMode       ! optimization mode
  real(r8), intent(in)  ::  p(:)          ! parameter vector
  integer             ::  start, aostart, aniso_count
  real(r8), allocatable :: aniso_p(:)

  start = 1
  aostart = 1
  select case(optMode)
  case(OPT_LIN, OPT_NUM_LIN) ! linear parameters
    alpha(2:umax) = p(start:start+unum-1)
    start = start + unum
    beta(2:xmax, 1:ncdiff) = reshape(p(start:start + xnum-1), [xpnum, ncdiff])
    start = start + xnum
    gamma(1:fpnum, 1:ncdiff) = reshape(p(start:start+fnum-1), [fpnum, ncdiff])
    start = start + fnum
    if (useAOJasTerms) then
      aniso_count = JasAniso%getNumberofParams()
      allocate(aniso_p(aniso_count))
      aniso_p = p(start:start+aniso_count-1)
      call JasAniso%setParamVector(optmode,aniso_p)
    endif
  case(OPT_NONLIN)
    select case(distType)
    case(DIST_SM, DIST_DOUBLEEXP)
      scaleEE = p(1)
      scaleEN(1:ncdiff) = p(2:ncdiff + 1)
      start = ncdiff+2
    case(DIST_NEEDS)
      scaleEE = p(1)
      powerEE = p(2)
      scaleEN(1:ncdiff) = p(3:ncdiff + 2)
      powerEN(1:ncdiff) = p(ncdiff + 3:2*ncdiff + 2)
      start = 2*ncdiff+3
    case default
      call abortp("putVectorIC: distType not implemented")
    end select
    !!!if (useAOJasTerms) gl(gnum) = p(start)
    !print*, p(:)
  case(OPT_ALL, OPT_NUM_ALL)
    call putVectorIC(OPT_LIN, p(1:numParams))
    call putVectorIC(OPT_NONLIN, p(numParams+1:))
  case default
    call abortp("putVectorIC: optMode not implemented")
  end select

  !call jasoutput_ic

end subroutine putVectorIC


!==================================  jasicInit  ======================================================


subroutine jasicInit(Rdu,mode)
!-----------------------------

! ic generic Jastrow: one electron updates only, no derivatives
! This version initializes the arrays for one electron updates
! Rdu keeps the auxiliary data for the update

   type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
   integer, intent(in)              :: mode

   ! values of correlation terms
   real(r8) :: uTerm, xTerm, fTerm, gTerm

   ! sum in x, u and f-term
   real(r8) :: sumr
   ! parameter value to satisfy e-e cusp condition
   real(r8) :: eeCusp

   ! commonly used terms
   real(r8) :: b, r, scale, power, Fsum, Gsum

   ! loop variables
   integer :: a, c, i, j, t

   !!!write(iul,'(a)') 'DBG:jasicInit:start'

   call Rdu%initENSize(nclast,fmax)

   uTerm = 0d0
   xTerm = 0d0
   fTerm = 0d0
   gTerm = 0d0

   Rdu%Fij = 0.d0   ! -> u+f
   Rdu%Gi = 0.d0    ! -> x+X

   if(.not.allocated(eePowers)) call abortp("JasicInit requires allocated module vars")

   call internal_precalculationOfTerms()
   ! if (useAOJasTerms) call JasAniso%getParamCount(curOptMode,anisoJ1_count,anisoJ2_count)

   ! value for e-e cusp
   eeCusp = 0.5d0
   select case(distType)
   case(DIST_SM, DIST_DOUBLEEXP)
      eeCusp = 0.5d0 / scaleEE
   case(DIST_NEEDS)
      eeCusp = 0.5d0 * scaleEE
   end select


   call internal_electronElectronTerms()

   call internal_electronNucleusTerms()

   !---------Electron-Electron-Nucleus-Correlation-Terms-------------------
   if (fmax > 2) then
      if (useGeneric) then
         call abortp("generic ic update not yet implemented")
      else
         call eenHardcodedInit(Rdu,fTerm)
      end if
   end if

   if (useAOJasTerms) call JasAniso%CalcOnlyU(Rdu,gTerm)

  !---------Calculation of U -------------------------
  Rdu%U = uTerm + xTerm + fTerm + gTerm

  !!!write(iul,'(a,6g20.10)') 'DBG:jasicInit:',Rdu%U,uTerm,xTerm,fTerm,uTerm+fTerm,gTerm
!    Fsum = 0; Gsum = 0
!    do i=1,ne
!       Gsum = Gsum + Rdu%Gi(i)
!       do j=i+1,ne
!          Fsum = Fsum + Rdu%Fij(j,i)
!          write(iul,'(a,2i3,g20.10)') 'DBG:jasicInit3:',i,j,Rdu%Fij(j,i)
!       enddo
!    enddo
   !!!write(iul,'(a,4g20.10)') 'DBG:jasicInit3:', Fsum,Gsum


  Rdu%U0 = Rdu%U
  Rdu%ieJasold = 0

  call Rdu%markJastrowValid()

contains

   subroutine internal_precalculationOfTerms()

      do i = 1, ne
         do a = 1, nclast
            c = atoms(a)%sa
            scale = scaleEN(c)
            power = powerEN(c)

            ! powers of e-n distances
            enPowers(a, i, -2) = 0d0
            enPowers(a, i, -1) = 0d0
            enPowers(a, i, 0) = 1d0

            select case(distType)
            case(DIST_SM)
               r = 1 / (1 + scale * Rdu%rai(a, i))
               enPowers(a, i, 1) = scale * Rdu%rai(a, i) * r

            case(DIST_DOUBLEEXP)
               r = exp(-scale * Rdu%rai(a, i))
               enPowers(a, i, 1) = 1 - r

            case(DIST_NEEDS)
               b = Rdu%rai(a, i) ** power
               r = 1 / (b + scale)
               enPowers(a, i, 1) = Rdu%rai(a, i) * r

            case default
               call abortp("Unknown Jastrow distance type in jasicall")
            end select

            do t = 2, max(xmax, fmax)
               enPowers(a, i, t) = enPowers(a, i, t-1) * enPowers(a, i, 1)
            enddo
         enddo

         do j = i + 1, ne
            ! powers of e-e distances
            eePowers(i, j, -2) = 0d0
            eePowers(i, j, -1) = 0d0
            eePowers(i, j, 0) = 1d0

            select case(distType)
            case(DIST_SM)
               r = 1 / (1 + scaleEE * Rdu%rij(i, j))
               eePowers(i, j, 1) = scaleEE * Rdu%rij(i, j) * r

            case(DIST_DOUBLEEXP)
               r = exp(-scaleEE * Rdu%rij(i, j))
               eePowers(i, j, 1) = 1 - r

            case(DIST_NEEDS)
               b = Rdu%rij(i, j) ** powerEE
               r = 1 / (b + scaleEE)
               eePowers(i, j, 1) = Rdu%rij(i, j) * r

            case default
               call abortp("Unknown Jastrow distance type in jasicall")
            end select

            do t = 2, max(umax, fmax)
               eePowers(i, j, t) = eePowers(i, j, t-1) * eePowers(i, j, 1)
            enddo
         enddo
      enddo

   end subroutine internal_precalculationOfTerms


   subroutine internal_electronElectronTerms()
      do i = 1, ne
         do j = i + 1, ne
            sumr  = eeCusp * eePowers(i, j, 1)
            if (diffeecusp .and. ((i<=nalpha) .eqv. (j<=nalpha))) then
               sumr  = 0.5d0*eeCusp * eePowers(i, j, 1)
            endif

            ! start at 2 to satisfy cusp
            do t = 2, umax
               sumr  = sumr  + alpha(t) * eePowers(i, j, t)
            enddo
            Rdu%Fij(j,i) = Rdu%Fij(j,i) + sumr
            uTerm = uTerm + sumr
         enddo
      enddo
   end subroutine internal_electronElectronTerms

   subroutine internal_electronNucleusTerms()
      do a = 1, nclast
         c = atoms(a)%sa
         do i = 1, ne
            sumr  = 0d0

            if(nucCusp .and. .not. useAOJasTerms) then
               sumr  = atoms(a)%za * enPowers(a, i, 1)
            endif

            ! start at 2 to satisfy cusp
            do t = 2, xmax
               sumr  = sumr  + beta(t, c) * enPowers(a, i, t)
            enddo
            Rdu%Gi(i) = Rdu%Gi(i) + sumr
            xTerm = xTerm + sumr
         enddo
      enddo
   end subroutine internal_electronNucleusTerms

end subroutine jasicInit


subroutine eenHardcodedInit(Rdu,fTerm)
  ! value of f term + derivs
   type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
   real(r8), intent(inout) :: fTerm

   real(r8) :: rComb(2)   !!! ???
   real(r8) :: tmp
   real(r8) :: terms(fpnum), termsn(6, fpnum), termsl(2, fpnum)
   real(r8) :: nfterms(fmax), nftermsn(6, fmax), nftermsl(2, fmax)
   integer :: a, c, g, i, j, t

   do a = 1, nclast
      c = atoms(a)%sa

      do i = 1, ne
         do j = i + 1, ne
            ! rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
            ! rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

            if(fmax > 2) then
               ! r_i^2 r_j + r_i r_j^2
               nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                          enPowers(a, i, 1) * enPowers(a, j, 2)

               ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
               tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                   2 * enPowers(a, i, 1) * enPowers(a, j, 1)
               nfterms(2) = -0.5d0 * eePowers(i, j, 1) * tmp

               ! r_ij^2 (r_i + r_j)
               tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
               nfterms(3) = eePowers(i, j, 2) * tmp

               do t = 1, 2
                  terms(t) = nfterms(t) - nfterms(3)
               enddo
            endif

            if(fmax > 3) then
               ! 2 ri^2 rj^2
               terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)

               ! rij^2 (ri^2 + rj^2)
               tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
               terms(4) = eePowers(i, j, 2) * tmp

               ! ri^3 rj + ri rj^3
               nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                          enPowers(a, i, 1) * enPowers(a, j, 3)

               ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
               tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                   enPowers(a, i, 2) * enPowers(a, j, 1) - &
                   enPowers(a, i, 1) * enPowers(a, j, 2)
               nfterms(2) = - eePowers(i, j, 1) * tmp

               ! rij^3 (ri + rj)
               tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
               nfterms(3) = eePowers(i, j, 3) * tmp

               ! rij^2 ri rj
               tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
               nfterms(4) = eePowers(i, j, 2) * tmp

               do t = 1, 3
                  terms(4+t) = nfterms(t) - nfterms(4)
               enddo
            endif

            if(fmax > 4) then
               ! ri^3 rj^2 + ri^2 rj^3
               terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                        enPowers(a, i, 2) * enPowers(a, j, 3)

               ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
               tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                   2 * enPowers(a, i, 2) * enPowers(a, j, 2)
               terms(9)  = eePowers(i, j, 1) * tmp

               ! rij^3 (ri^2 + rj^2)
               tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
               terms(10) = eePowers(i, j, 3) * tmp

               ! rij^2 (ri^3 + rj^3)
               tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
               terms(11) = eePowers(i, j, 2) * tmp

               ! ri^4 rj + ri rj^4
               nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                          enPowers(a, i, 1) * enPowers(a, j, 4)
               ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
               tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                   enPowers(a, i, 3) * enPowers(a, j, 1) - &
                   enPowers(a, i, 1) * enPowers(a, j, 3)
               nfterms(2) = -eePowers(i, j, 1) * tmp

               ! rij^4 (ri + rj)
               tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
               nfterms(3) = eePowers(i, j, 4) * tmp

               ! rij^2 (ri^2 rj + ri rj^2)
               tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                   enPowers(a, i, 1) * enPowers(a, j, 2)
               nfterms(4) = eePowers(i, j, 2) * tmp

               ! rij^3 ri rj
               tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
               nfterms(5) = eePowers(i, j, 3) * tmp

               do t = 1, 4
                  terms(11+t) = nfterms(t) - nfterms(5)
                  termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
                  termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
               enddo
            endif

            if(fmax > 5) then
               ! ri^2 rj^4 + ri^4 rj^2
               terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                         enPowers(a, i, 2) * enPowers(a, j, 4)

               ! 2 ri^3 rj^3
               terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

               ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
               tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                   enPowers(a, i, 3) * enPowers(a, j, 2) - &
                   enPowers(a, i, 2) * enPowers(a, j, 3)
               terms(18) = eePowers(i, j, 1) * tmp

               ! rij^4 (ri^2 + rj^2)
               tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
               terms(19) = eePowers(i, j, 4) * tmp

               ! 2 rij^2 ri^2 rj^2
               tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
               terms(20) = 2 * eePowers(i, j, 2) * tmp

               ! rij^3 (ri^3 + rj^3)
               tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
               terms(21) = eePowers(i, j, 3) * tmp

               ! rij^2 (ri^4 + rj^4)
               tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
               terms(22) = eePowers(i, j, 2) * tmp

               ! ri^5 rj + ri rj^5
               nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                          enPowers(a, i, 1) * enPowers(a, j, 5)

               ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
               tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                   enPowers(a, i, 4) * enPowers(a, j, 1) - &
                   enPowers(a, i, 1) * enPowers(a, j, 4)
               nfterms(2) = - eePowers(i, j, 1) * tmp

               ! rij^5 (ri + rj)
               tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
               nfterms(3) = eePowers(i, j, 5) * tmp

               ! rij^2 (ri rj^3 + ri^3 rj)
               tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                   enPowers(a, i, 1) * enPowers(a, j, 3)
               nfterms(4) = eePowers(i, j, 2) * tmp

               ! rij^3 (ri rj^2 + ri^2 rj)
               tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                   enPowers(a, i, 1) * enPowers(a, j, 2)
               nfterms(5) = eePowers(i, j, 3) * tmp

               ! rij^4 ri rj
               tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
               nfterms(6) = eePowers(i, j, 4) * tmp

               do t = 1, 5
                  terms(22+t) = nfterms(t) - nfterms(6)
               enddo
            endif

            do g = 1, fpnum
               fTerm = fTerm + gamma(g, c) * terms(g)
               Rdu%Fij(j,i) = Rdu%Fij(j,i) + gamma(g, c) * terms(g)
            enddo
         enddo  ! j
      enddo  ! i
   enddo  ! a

end subroutine eenHardcodedInit




subroutine jasicInitWithUk(Rdu)
!------------------------------

! ic generic Jastrow: one electron updates only, no derivatives
! This version initializes the arrays for one electron updates
! This version also for param derivs
! Rdu keeps the auxiliary data for the update

   type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations

   ! values of correlation terms
   real(r8) :: uTerm, xTerm, fTerm, gTerm

   ! sum in x, u and f-term
   real(r8) :: sumr
   ! parameter value to satisfy e-e cusp condition
   real(r8) :: eeCusp

   ! commonly used terms
   real(r8) :: b, r, scale, power, Fsum, Gsum

   ! loop variables
   integer :: a, c, i, j, t, k
   integer :: offsetJ1, offsetJ2

   !!!write(iul,'(a)') 'DBG:jasicInit:start'

   call Rdu%initENSize(nclast,fmax)

   uTerm = 0d0
   xTerm = 0d0
   fTerm = 0d0
   gTerm = 0d0

   Rdu%Fij = 0.d0   ! -> u+f
   Rdu%Gi = 0.d0    ! -> x+X
   Rdu%Fijk = 0.d0  ! parameter derivs
   Rdu%Gki = 0.d0

   if(.not.allocated(eePowers)) call abortp("JasicInit requires allocated module vars")

   call internal_precalculationOfTerms()

   ! value for e-e cusp
   eeCusp = 0.5d0
   select case(distType)
   case(DIST_SM, DIST_DOUBLEEXP)
      eeCusp = 0.5d0 / scaleEE
   case(DIST_NEEDS)
      eeCusp = 0.5d0 * scaleEE
   end select


   call internal_electronElectronTerms()

   call internal_electronNucleusTerms()

   !---------Electron-Electron-Nucleus-Correlation-Terms-------------------
   if (fmax > 2) then
      if (useGeneric) then
         call abortp("generic ic update not yet implemented")
      else
         call eenHardcodedInitWithUk(Rdu,fTerm)
      end if
   end if

   ! AO term contributions
   offsetJ1 = xnum
   offsetJ2 = unum+fnum
   if (useAOJasTerms)  call JasAniso%CalcOnlyUandUK(Rdu,gTerm,offsetJ1,offsetJ2)

   !---------Calculation of U -------------------------
   Rdu%U = uTerm + xTerm + fTerm + gTerm

   !!!write(iul,'(a,6g20.10)') 'DBG:jasicInit:',Rdu%U,uTerm,xTerm,fTerm,uTerm+fTerm,gTerm
!    Fsum = 0; Gsum = 0
!    do i=1,ne
!       Gsum = Gsum + Rdu%Gi(i)
!       do j=i+1,ne
!          Fsum = Fsum + Rdu%Fij(j,i)
!          write(iul,'(a,2i3,g20.10)') 'DBG:jasicInit3:',i,j,Rdu%Fij(j,i)
!       enddo
!    enddo
   !!!write(iul,'(a,4g20.10)') 'DBG:jasicInit3:', Fsum,Gsum

   !call JasAniso%getParamCount(curOptMode,anisoJ1_count,anisoJ2_count)
   do i=1,ne
      do j=i+1,ne
        !isotropic ee
         Rdu%Uk(1:unum) = Rdu%Uk(1:unum) + Rdu%Fijk(j,i,1:unum)
        !isotropic een
         Rdu%Uk(unum+xnum+1:unum+xnum+fnum) = Rdu%Uk(unum+xnum+1:unum+xnum+fnum) +&
             Rdu%Fijk(j,i,unum+1:unum+fnum)
        !anisotropic two body terms
         Rdu%Uk(unum+xnum+fnum+anisoJ1_count+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) = &
            Rdu%Uk(unum+xnum+fnum+anisoJ1_count+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) + &
                Rdu%Fijk(j,i,unum+fnum+1:unum+fnum+anisoJ2_count)

      enddo
   enddo
   do i=1,ne
      !isotropic en
      Rdu%Uk(unum+1:unum+xnum) = Rdu%Uk(unum+1:unum+xnum) + Rdu%Gki(1:xnum,i)
      !anisotropic en
      Rdu%Uk(unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count) = &
         Rdu%Uk(unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count) + Rdu%Gki(xnum+1:xnum+anisoJ1_count,i)
   enddo
   !!!write(iul,'(a,g20.10)') 'DBG:jasicInit:U,nums,Uk',Rdu%U
   !!!write(iul,'(4i4)') unum,xnum,fnum,gnum
   !!!write(iul,'(5g20.10)') Rdu%Uk

   Rdu%U0 = Rdu%U
   Rdu%Uk0 = Rdu%Uk
   Rdu%ieJasold = 0

   call Rdu%markJastrowValid()

contains

   subroutine internal_precalculationOfTerms()

      do i = 1, ne
         do a = 1, nclast
            c = atoms(a)%sa
            scale = scaleEN(c)
            power = powerEN(c)

            ! powers of e-n distances
            enPowers(a, i, -2) = 0d0
            enPowers(a, i, -1) = 0d0
            enPowers(a, i, 0) = 1d0

            select case(distType)
            case(DIST_SM)
               r = 1 / (1 + scale * Rdu%rai(a, i))
               enPowers(a, i, 1) = scale * Rdu%rai(a, i) * r

            case(DIST_DOUBLEEXP)
               r = exp(-scale * Rdu%rai(a, i))
               enPowers(a, i, 1) = 1 - r

            case(DIST_NEEDS)
               b = Rdu%rai(a, i) ** power
               r = 1 / (b + scale)
               enPowers(a, i, 1) = Rdu%rai(a, i) * r

            case default
               call abortp("Unknown Jastrow distance type in jasicall")
            end select

            do t = 2, max(xmax, fmax)
               enPowers(a, i, t) = enPowers(a, i, t-1) * enPowers(a, i, 1)
            enddo
         enddo

         do j = i + 1, ne
            ! powers of e-e distances
            eePowers(i, j, -2) = 0d0
            eePowers(i, j, -1) = 0d0
            eePowers(i, j, 0) = 1d0

            select case(distType)
            case(DIST_SM)
               r = 1 / (1 + scaleEE * Rdu%rij(i, j))
               eePowers(i, j, 1) = scaleEE * Rdu%rij(i, j) * r

            case(DIST_DOUBLEEXP)
               r = exp(-scaleEE * Rdu%rij(i, j))
               eePowers(i, j, 1) = 1 - r

            case(DIST_NEEDS)
               b = Rdu%rij(i, j) ** powerEE
               r = 1 / (b + scaleEE)
               eePowers(i, j, 1) = Rdu%rij(i, j) * r

            case default
               call abortp("Unknown Jastrow distance type in jasicall")
            end select

            do t = 2, max(umax, fmax)
               eePowers(i, j, t) = eePowers(i, j, t-1) * eePowers(i, j, 1)
            enddo
         enddo
      enddo

   end subroutine internal_precalculationOfTerms


   subroutine internal_electronElectronTerms()
      do i = 1, ne
         do j = i + 1, ne
            sumr  = eeCusp * eePowers(i, j, 1)
            if (diffeecusp .and. ((i<=nalpha) .eqv. (j<=nalpha))) then
               sumr  = 0.5d0*eeCusp * eePowers(i, j, 1)
            endif

            ! start at 2 to satisfy cusp
            do t = 2, umax
               sumr  = sumr  + alpha(t) * eePowers(i, j, t)
               Rdu%Fijk(j,i,t-1) = eePowers(i,j,t)
            enddo
            Rdu%Fij(j,i) = Rdu%Fij(j,i) + sumr
            uTerm = uTerm + sumr
         enddo
      enddo
   end subroutine internal_electronElectronTerms


   subroutine internal_electronNucleusTerms()
      do a = 1, nclast
         c = atoms(a)%sa
         do i = 1, ne
            sumr  = 0d0

            if(nucCusp .and. .not. useAOJasTerms) then
               sumr  = atoms(a)%za * enPowers(a, i, 1)
            endif

            ! start at 2 to satisfy cusp
            do t = 2, xmax
               sumr  = sumr  + beta(t, c) * enPowers(a, i, t)
               k = (c-1)*xpnum + t-1        ! see setNumTerms() below
               Rdu%Gki(k,i) = Rdu%Gki(k,i) + enPowers(a, i, t)
            enddo
            Rdu%Gi(i) = Rdu%Gi(i) + sumr
            xTerm = xTerm + sumr
         enddo
      enddo
   end subroutine internal_electronNucleusTerms

end subroutine jasicInitWithUk


subroutine eenHardcodedInitWithUk(Rdu,fTerm)
  ! value of f term + derivs
   type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
   real(r8), intent(inout) :: fTerm

   real(r8) :: rComb(2)   !!! ???
   real(r8) :: tmp
   real(r8) :: terms(fpnum), termsn(6, fpnum), termsl(2, fpnum)
   real(r8) :: nfterms(fmax), nftermsn(6, fmax), nftermsl(2, fmax)
   integer :: a, c, g, i, j, t, k

   do a = 1, nclast
      c = atoms(a)%sa

      do i = 1, ne
         do j = i + 1, ne
            ! rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
            ! rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

            if(fmax > 2) then
               ! r_i^2 r_j + r_i r_j^2
               nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                          enPowers(a, i, 1) * enPowers(a, j, 2)

               ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
               tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                   2 * enPowers(a, i, 1) * enPowers(a, j, 1)
               nfterms(2) = -0.5d0 * eePowers(i, j, 1) * tmp

               ! r_ij^2 (r_i + r_j)
               tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
               nfterms(3) = eePowers(i, j, 2) * tmp

               do t = 1, 2
                  terms(t) = nfterms(t) - nfterms(3)
               enddo
            endif

            if(fmax > 3) then
               ! 2 ri^2 rj^2
               terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)

               ! rij^2 (ri^2 + rj^2)
               tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
               terms(4) = eePowers(i, j, 2) * tmp

               ! ri^3 rj + ri rj^3
               nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                          enPowers(a, i, 1) * enPowers(a, j, 3)

               ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
               tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                   enPowers(a, i, 2) * enPowers(a, j, 1) - &
                   enPowers(a, i, 1) * enPowers(a, j, 2)
               nfterms(2) = - eePowers(i, j, 1) * tmp

               ! rij^3 (ri + rj)
               tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
               nfterms(3) = eePowers(i, j, 3) * tmp

               ! rij^2 ri rj
               tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
               nfterms(4) = eePowers(i, j, 2) * tmp

               do t = 1, 3
                  terms(4+t) = nfterms(t) - nfterms(4)
               enddo
            endif

            if(fmax > 4) then
               ! ri^3 rj^2 + ri^2 rj^3
               terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                        enPowers(a, i, 2) * enPowers(a, j, 3)

               ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
               tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                   2 * enPowers(a, i, 2) * enPowers(a, j, 2)
               terms(9)  = eePowers(i, j, 1) * tmp

               ! rij^3 (ri^2 + rj^2)
               tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
               terms(10) = eePowers(i, j, 3) * tmp

               ! rij^2 (ri^3 + rj^3)
               tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
               terms(11) = eePowers(i, j, 2) * tmp

               ! ri^4 rj + ri rj^4
               nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                          enPowers(a, i, 1) * enPowers(a, j, 4)
               ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
               tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                   enPowers(a, i, 3) * enPowers(a, j, 1) - &
                   enPowers(a, i, 1) * enPowers(a, j, 3)
               nfterms(2) = -eePowers(i, j, 1) * tmp

               ! rij^4 (ri + rj)
               tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
               nfterms(3) = eePowers(i, j, 4) * tmp

               ! rij^2 (ri^2 rj + ri rj^2)
               tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                   enPowers(a, i, 1) * enPowers(a, j, 2)
               nfterms(4) = eePowers(i, j, 2) * tmp

               ! rij^3 ri rj
               tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
               nfterms(5) = eePowers(i, j, 3) * tmp

               do t = 1, 4
                  terms(11+t) = nfterms(t) - nfterms(5)
                  termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
                  termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
               enddo
            endif

            if(fmax > 5) then
               ! ri^2 rj^4 + ri^4 rj^2
               terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                         enPowers(a, i, 2) * enPowers(a, j, 4)

               ! 2 ri^3 rj^3
               terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

               ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
               tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                   enPowers(a, i, 3) * enPowers(a, j, 2) - &
                   enPowers(a, i, 2) * enPowers(a, j, 3)
               terms(18) = eePowers(i, j, 1) * tmp

               ! rij^4 (ri^2 + rj^2)
               tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
               terms(19) = eePowers(i, j, 4) * tmp

               ! 2 rij^2 ri^2 rj^2
               tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
               terms(20) = 2 * eePowers(i, j, 2) * tmp

               ! rij^3 (ri^3 + rj^3)
               tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
               terms(21) = eePowers(i, j, 3) * tmp

               ! rij^2 (ri^4 + rj^4)
               tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
               terms(22) = eePowers(i, j, 2) * tmp

               ! ri^5 rj + ri rj^5
               nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                          enPowers(a, i, 1) * enPowers(a, j, 5)

               ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
               tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                   enPowers(a, i, 4) * enPowers(a, j, 1) - &
                   enPowers(a, i, 1) * enPowers(a, j, 4)
               nfterms(2) = - eePowers(i, j, 1) * tmp

               ! rij^5 (ri + rj)
               tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
               nfterms(3) = eePowers(i, j, 5) * tmp

               ! rij^2 (ri rj^3 + ri^3 rj)
               tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                   enPowers(a, i, 1) * enPowers(a, j, 3)
               nfterms(4) = eePowers(i, j, 2) * tmp

               ! rij^3 (ri rj^2 + ri^2 rj)
               tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                   enPowers(a, i, 1) * enPowers(a, j, 2)
               nfterms(5) = eePowers(i, j, 3) * tmp

               ! rij^4 ri rj
               tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
               nfterms(6) = eePowers(i, j, 4) * tmp

               do t = 1, 5
                  terms(22+t) = nfterms(t) - nfterms(6)
               enddo
            endif

            do g = 1, fpnum
               fTerm = fTerm + gamma(g, c) * terms(g)
               Rdu%Fij(j,i) = Rdu%Fij(j,i) + gamma(g, c) * terms(g)
               k = unum + (c-1)*fpnum + g        ! see setNumTerms() below
               !!!write(iul,'(a,5i4,g20.10)') 'DBG:eenhcwuk:',k,unum,ncdiff,c,g,terms(g)
               Rdu%Fijk(j,i,k) = Rdu%Fijk(j,i,k) + terms(g)
            enddo
         enddo  ! j
      enddo  ! i
   enddo  ! a

end subroutine eenHardcodedInitWithUk





!==================================  jasicUpdate  ======================================================


subroutine jasicUpdate(Rdu, ie)
!------------------------------

! ic generic Jastrow: one electron updates only, no derivatives
! This routine calculates U by one electron updates
! Rdu keeps the auxiliary data for the update

   type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
   integer, intent(in)              :: ie     ! update for electron ie

   ! values of correlation terms
   real(r8) :: uTerm, xTerm, fTerm

   ! sum in x, u and f-term
   real(r8) :: sumr
   ! parameter value to satisfy e-e cusp condition
   real(r8) :: eeCusp

   ! commonly used terms
   real(r8) :: b, r, scale, power, Fsum, Gsum

   ! loop variables
   integer :: a, c, i, j, t

   !uTerm = 0d0   ! -> F
   !xTerm = 0d0   ! -> G
   !fTerm = 0d0   ! -> H
   !gTerm = 0d0   ! -> X

   !!!write(iul,'(a,g20.10,2i3)') 'DBG:jasicUpdate:start', Rdu%U,Rdu%ieJasold,ie

   if (ie /= Rdu%ieJasold) then
      if (Rdu%ieJasold /= 0) then
        ! necessary to restore one-electron terms for previous electron ieold
        ! also restore U to original U
        Rdu%Fij(Rdu%ieJasold,1:Rdu%ieJasold-1) = Rdu%Fijold(1:Rdu%ieJasold-1)
        Rdu%Fij(Rdu%ieJasold+1:ne,Rdu%ieJasold) = Rdu%Fijold(Rdu%ieJasold+1:ne)
        Rdu%Gi(Rdu%ieJasold) = Rdu%Giold
        enPowers(1:Rdu%enDim1,Rdu%ieJasold,1:Rdu%enDim2) = Rdu%raibarOld(:,:)
        Rdu%U = Rdu%U0
      endif
      Rdu%Fijold(1:ie-1) = Rdu%Fij(ie,1:ie-1)
      Rdu%Fijold(ie+1:ne) = Rdu%Fij(ie+1:ne,ie)
      Rdu%Giold = Rdu%Gi(ie)
      Rdu%raibarOld(:,:) = enPowers(1:Rdu%enDim1,ie,1:Rdu%enDim2)
      Rdu%ieJasold = ie
   ! else same electron as before, simply update ie
   endif

   ! subtract current=previous values the three Jastrow terms
   Fsum = 0.d0
   do i=1,ie-1
      Fsum = Fsum + Rdu%Fij(ie,i)
   enddo
   do j=ie+1,ne
      Fsum = Fsum + Rdu%Fij(j,ie)
   enddo
   Rdu%U = Rdu%U - Fsum - Rdu%Gi(ie)

   !!!write(iul,'(a,4g20.10)') 'DBG:jasicUpdate1:', Rdu%U,Fsum,Rdu%Gi(ie),Rdu%Xi(ie)

   if (.not.allocated(eePowers)) call abortp("JasicUpdate requires allocated module vars")

   call internal_precalculationOfTerms()   ! ie only !

   ! value for e-e cusp
   eeCusp = 0.5d0
   select case(distType)
   case(DIST_SM, DIST_DOUBLEEXP)
      eeCusp = 0.5d0 / scaleEE
   case(DIST_NEEDS)
      eeCusp = 0.5d0 * scaleEE
   end select

   ! only for the current electron ie !
   i = ie

   ! recalculate the ie terms
   Rdu%Fij(i,1:i-1) = 0.d0
   Rdu%Fij(i+1:ne,i) = 0.d0
   Rdu%Gi(i) = 0.d0

   !---------Electron-Electron-Correlation-Terms-------------------
   call internal_electronElectronTerms()   ! ie only !

   !---------Electron-Nucleus-Correlation-Terms-------------------
   call internal_electronNucleusTerms()    ! ie only !

   !---------Electron-Electron-Nucleus-Correlation-Terms-------------------
   if (fmax > 2) then
      if (useGeneric) then
         call abortp("generic ic update not yet implemented")
      else
         call eenHardcodedUpdate(Rdu,ie)
      end if
   end if

   !-------- AO term contributions
   if (useAOJasTerms) call JasAniso%UpdateOnlyU(Rdu,ie)

   !--------- update of U -------------------------
   ! add updated values to U
   Fsum = 0.d0
   do i=1,ie-1
      Fsum = Fsum + Rdu%Fij(ie,i)
   enddo
   do j=ie+1,ne
      Fsum = Fsum + Rdu%Fij(j,ie)
   enddo
   Rdu%U = Rdu%U + Fsum + Rdu%Gi(ie)

   !!!write(iul,'(a,4g20.10)') 'DBG:jasicUpdate2:', Rdu%U,Fsum,Rdu%Gi(ie),Rdu%Xi(ie)

!    Fsum = 0; Gsum = 0
!    do i=1,ne
!       Gsum = Gsum + Rdu%Gi(i)
!       do j=i+1,ne
!          Fsum = Fsum + Rdu%Fij(j,i)
!          !!!write(iul,'(a,2i3,g20.10)') 'DBG:jasicUpdate3:',i,j,Rdu%Fij(j,i)

!       enddo
!    enddo
!    write(iul,'(a,4g20.10)') 'DBG:jasicUpdate3:', Fsum,Gsum

   call Rdu%markJastrowValid()

contains

   subroutine internal_precalculationOfTerms()

      do a = 1, nclast
         c = atoms(a)%sa
         scale = scaleEN(c)
         power = powerEN(c)

         ! powers of e-n distances
         enPowers(a, i, -2) = 0d0
         enPowers(a, i, -1) = 0d0
         enPowers(a, i, 0) = 1d0

         select case(distType)
         case(DIST_SM)
            r = 1 / (1 + scale * Rdu%rai(a, i))
            enPowers(a, i, 1) = scale * Rdu%rai(a, i) * r

         case(DIST_DOUBLEEXP)
            r = exp(-scale * Rdu%rai(a, i))
            enPowers(a, i, 1) = 1 - r

         case(DIST_NEEDS)
            b = Rdu%rai(a, i) ** power
            r = 1 / (b + scale)
            enPowers(a, i, 1) = Rdu%rai(a, i) * r

         case default
            call abortp("Unknown Jastrow distance type in jasicall")
         end select

         do t = 2, max(xmax, fmax)
            enPowers(a, i, t) = enPowers(a, i, t-1) * enPowers(a, i, 1)
         enddo
      enddo

      ! note: ie in r_ij can be i or j ! Here case j < ie
      do j = 1, i-1
         ! powers of e-e distances
         eePowers(j, i, -2) = 0d0
         eePowers(j, i, -1) = 0d0
         eePowers(j, i, 0) = 1d0

         select case(distType)
         case(DIST_SM)
            r = 1 / (1 + scaleEE * Rdu%rij(j, i))
            eePowers(j, i, 1) = scaleEE * Rdu%rij(j, i) * r

         case(DIST_DOUBLEEXP)
            r = exp(-scaleEE * Rdu%rij(j, i))
            eePowers(j, i, 1) = 1 - r

         case(DIST_NEEDS)
            b = Rdu%rij(j, i) ** powerEE
            r = 1 / (b + scaleEE)
            eePowers(j, i, 1) = Rdu%rij(j, i) * r

         case default
            call abortp("Unknown Jastrow distance type in jasicall")
         end select

         do t = 2, max(umax, fmax)
            eePowers(j, i, t) = eePowers(j, i, t-1) * eePowers(j, i, 1)
         enddo
      enddo
      ! now j > ie
      do j = i + 1, ne
         ! powers of e-e distances
         eePowers(i, j, -2) = 0d0
         eePowers(i, j, -1) = 0d0
         eePowers(i, j, 0) = 1d0

         select case(distType)
         case(DIST_SM)
            r = 1 / (1 + scaleEE * Rdu%rij(i, j))
            eePowers(i, j, 1) = scaleEE * Rdu%rij(i, j) * r

         case(DIST_DOUBLEEXP)
            r = exp(-scaleEE * Rdu%rij(i, j))
            eePowers(i, j, 1) = 1 - r

         case(DIST_NEEDS)
            b = Rdu%rij(i, j) ** powerEE
            r = 1 / (b + scaleEE)
            eePowers(i, j, 1) = Rdu%rij(i, j) * r

         case default
            call abortp("Unknown Jastrow distance type in jasicall")
         end select

         do t = 2, max(umax, fmax)
            eePowers(i, j, t) = eePowers(i, j, t-1) * eePowers(i, j, 1)
         enddo
      enddo

   end subroutine internal_precalculationOfTerms


   subroutine internal_electronElectronTerms()
      do j = 1, i-1
         sumr  = eeCusp * eePowers(j, i, 1)
         if (diffeecusp .and. ((i<=nalpha) .eqv. (j<=nalpha))) then
            sumr  = 0.5d0*eeCusp * eePowers(j, i, 1)
         endif

         ! start at 2 to satisfy cusp
         do t = 2, umax
            sumr  = sumr  + alpha(t) * eePowers(j, i, t)
         enddo
         Rdu%Fij(i,j) = Rdu%Fij(i,j) + sumr
         !uTerm = uTerm + sumr
      enddo
      do j = i+1, ne
         sumr  = eeCusp * eePowers(i, j, 1)
         if (diffeecusp .and. ((i<=nalpha) .eqv. (j<=nalpha))) then
            sumr  = 0.5d0*eeCusp * eePowers(i, j, 1)
         endif

         ! start at 2 to satisfy cusp
         do t = 2, umax
            sumr  = sumr  + alpha(t) * eePowers(i, j, t)
         enddo
         Rdu%Fij(j,i) = Rdu%Fij(j,i) + sumr
         !uTerm = uTerm + sumr
      enddo
   end subroutine internal_electronElectronTerms


   subroutine internal_electronNucleusTerms()
      do a = 1, nclast
         c = atoms(a)%sa
         sumr  = 0d0

         if(nucCusp .and. .not. useAOJasTerms) then
            sumr  = atoms(a)%za * enPowers(a, i, 1)
         endif

         ! start at 2 to satisfy cusp
         do t = 2, xmax
            sumr  = sumr  + beta(t, c) * enPowers(a, i, t)
         enddo
         Rdu%Gi(i) = Rdu%Gi(i) + sumr
         !xTerm = xTerm + sumr
      enddo
   end subroutine internal_electronNucleusTerms

end subroutine jasicUpdate


subroutine eenHardcodedUpdate(Rdu,ie)
  ! value of f term
   type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
   integer, intent(in) :: ie                  ! electron to update

   real(r8) :: rComb(2)   !!! ???
   real(r8) :: tmp
   real(r8) :: terms(fpnum), termsn(6, fpnum), termsl(2, fpnum)
   real(r8) :: nfterms(fmax), nftermsn(6, fmax), nftermsl(2, fmax)
   integer :: a, c, g, i, j, t

   i = ie ! ie only !

   do a = 1, nclast
      c = atoms(a)%sa

      do j = 1, i-1
         ! rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
         ! rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

         if(fmax > 2) then
            ! r_i^2 r_j + r_i r_j^2
            nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                       enPowers(a, i, 1) * enPowers(a, j, 2)

            ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                2 * enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(2) = -0.5d0 * eePowers(j, i, 1) * tmp

            ! r_ij^2 (r_i + r_j)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(j, i, 2) * tmp

            do t = 1, 2
               terms(t) = nfterms(t) - nfterms(3)
            enddo
         endif

         if(fmax > 3) then
            ! 2 ri^2 rj^2
            terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)

            ! rij^2 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(4) = eePowers(j, i, 2) * tmp

            ! ri^3 rj + ri rj^3
            nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 3)

            ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                enPowers(a, i, 2) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(2) = - eePowers(j, i, 1) * tmp

            ! rij^3 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(j, i, 3) * tmp

            ! rij^2 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(4) = eePowers(j, i, 2) * tmp

            do t = 1, 3
               terms(4+t) = nfterms(t) - nfterms(4)
            enddo
         endif

         if(fmax > 4) then
            ! ri^3 rj^2 + ri^2 rj^3
            terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                     enPowers(a, i, 2) * enPowers(a, j, 3)

            ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
            tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                2 * enPowers(a, i, 2) * enPowers(a, j, 2)
            terms(9)  = eePowers(j, i, 1) * tmp

            ! rij^3 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(10) = eePowers(j, i, 3) * tmp

            ! rij^2 (ri^3 + rj^3)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
            terms(11) = eePowers(j, i, 2) * tmp

            ! ri^4 rj + ri rj^4
            nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 4)
            ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
            tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                enPowers(a, i, 3) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 3)
            nfterms(2) = -eePowers(j, i, 1) * tmp

            ! rij^4 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(j, i, 4) * tmp

            ! rij^2 (ri^2 rj + ri rj^2)
            tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(4) = eePowers(j, i, 2) * tmp

            ! rij^3 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(5) = eePowers(j, i, 3) * tmp

            do t = 1, 4
               terms(11+t) = nfterms(t) - nfterms(5)
               termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
               termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
            enddo
         endif

         if(fmax > 5) then
            ! ri^2 rj^4 + ri^4 rj^2
            terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                      enPowers(a, i, 2) * enPowers(a, j, 4)

            ! 2 ri^3 rj^3
            terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

            ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
            tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 3) * enPowers(a, j, 2) - &
                enPowers(a, i, 2) * enPowers(a, j, 3)
            terms(18) = eePowers(j, i, 1) * tmp

            ! rij^4 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(19) = eePowers(j, i, 4) * tmp

            ! 2 rij^2 ri^2 rj^2
            tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
            terms(20) = 2 * eePowers(j, i, 2) * tmp

            ! rij^3 (ri^3 + rj^3)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
            terms(21) = eePowers(j, i, 3) * tmp

            ! rij^2 (ri^4 + rj^4)
            tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
            terms(22) = eePowers(j, i, 2) * tmp

            ! ri^5 rj + ri rj^5
            nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 5)

            ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
            tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 4) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 4)
            nfterms(2) = - eePowers(j, i, 1) * tmp

            ! rij^5 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(j, i, 5) * tmp

            ! rij^2 (ri rj^3 + ri^3 rj)
            tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 3)
            nfterms(4) = eePowers(j, i, 2) * tmp

            ! rij^3 (ri rj^2 + ri^2 rj)
            tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(5) = eePowers(j, i, 3) * tmp

            ! rij^4 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(6) = eePowers(j, i, 4) * tmp

            do t = 1, 5
               terms(22+t) = nfterms(t) - nfterms(6)
            enddo
         endif

         do g = 1, fpnum
            !fTerm = fTerm + gamma(g, c) * terms(g)
            Rdu%Fij(i,j) = Rdu%Fij(i,j) + gamma(g, c) * terms(g)
         enddo
      enddo  ! j
      do j = i+1, ne
         ! rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
         ! rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

         if(fmax > 2) then
            ! r_i^2 r_j + r_i r_j^2
            nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                       enPowers(a, i, 1) * enPowers(a, j, 2)

            ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                2 * enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(2) = -0.5d0 * eePowers(i, j, 1) * tmp

            ! r_ij^2 (r_i + r_j)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(i, j, 2) * tmp

            do t = 1, 2
               terms(t) = nfterms(t) - nfterms(3)
            enddo
         endif

         if(fmax > 3) then
            ! 2 ri^2 rj^2
            terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)

            ! rij^2 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(4) = eePowers(i, j, 2) * tmp

            ! ri^3 rj + ri rj^3
            nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 3)

            ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                enPowers(a, i, 2) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(2) = - eePowers(i, j, 1) * tmp

            ! rij^3 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(i, j, 3) * tmp

            ! rij^2 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(4) = eePowers(i, j, 2) * tmp

            do t = 1, 3
               terms(4+t) = nfterms(t) - nfterms(4)
            enddo
         endif

         if(fmax > 4) then
            ! ri^3 rj^2 + ri^2 rj^3
            terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                     enPowers(a, i, 2) * enPowers(a, j, 3)

            ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
            tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                2 * enPowers(a, i, 2) * enPowers(a, j, 2)
            terms(9)  = eePowers(i, j, 1) * tmp

            ! rij^3 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(10) = eePowers(i, j, 3) * tmp

            ! rij^2 (ri^3 + rj^3)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
            terms(11) = eePowers(i, j, 2) * tmp

            ! ri^4 rj + ri rj^4
            nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 4)
            ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
            tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                enPowers(a, i, 3) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 3)
            nfterms(2) = -eePowers(i, j, 1) * tmp

            ! rij^4 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(i, j, 4) * tmp

            ! rij^2 (ri^2 rj + ri rj^2)
            tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(4) = eePowers(i, j, 2) * tmp

            ! rij^3 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(5) = eePowers(i, j, 3) * tmp

            do t = 1, 4
               terms(11+t) = nfterms(t) - nfterms(5)
               termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
               termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
            enddo
         endif

         if(fmax > 5) then
            ! ri^2 rj^4 + ri^4 rj^2
            terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                      enPowers(a, i, 2) * enPowers(a, j, 4)

            ! 2 ri^3 rj^3
            terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

            ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
            tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 3) * enPowers(a, j, 2) - &
                enPowers(a, i, 2) * enPowers(a, j, 3)
            terms(18) = eePowers(i, j, 1) * tmp

            ! rij^4 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(19) = eePowers(i, j, 4) * tmp

            ! 2 rij^2 ri^2 rj^2
            tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
            terms(20) = 2 * eePowers(i, j, 2) * tmp

            ! rij^3 (ri^3 + rj^3)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
            terms(21) = eePowers(i, j, 3) * tmp

            ! rij^2 (ri^4 + rj^4)
            tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
            terms(22) = eePowers(i, j, 2) * tmp

            ! ri^5 rj + ri rj^5
            nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 5)

            ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
            tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 4) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 4)
            nfterms(2) = - eePowers(i, j, 1) * tmp

            ! rij^5 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(i, j, 5) * tmp

            ! rij^2 (ri rj^3 + ri^3 rj)
            tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 3)
            nfterms(4) = eePowers(i, j, 2) * tmp

            ! rij^3 (ri rj^2 + ri^2 rj)
            tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(5) = eePowers(i, j, 3) * tmp

            ! rij^4 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(6) = eePowers(i, j, 4) * tmp

            do t = 1, 5
               terms(22+t) = nfterms(t) - nfterms(6)
            enddo
         endif

         do g = 1, fpnum
            !fTerm = fTerm + gamma(g, c) * terms(g)
            Rdu%Fij(j,i) = Rdu%Fij(j,i) + gamma(g, c) * terms(g)
         enddo
      enddo  ! j
   enddo  ! a

end subroutine eenHardcodedUpdate


subroutine jasicUpdateWithUk(Rdu, ie)
!------------------------------------

! ic generic Jastrow: one electron updates only, no derivatives
! This routine calculates U by one electron updates
! Rdu keeps the auxiliary data for the update

   type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
   integer, intent(in)              :: ie     ! update for electron ie

   ! values of correlation terms
   real(r8) :: uTerm, xTerm, fTerm

   ! sum in x, u and f-term
   real(r8) :: sumr
   ! parameter value to satisfy e-e cusp condition
   real(r8) :: eeCusp

   ! commonly used terms
   !real(r8) :: b, r, s, v, tmpn(3), tmpl, scale, power, f, za, Fsum, Gsum, Fksum(unum+fnum+geennum+geennnum)
   real(r8) :: b, r, scale, power
   real(r8) :: Fsum, Gsum, Fksum(unum+fnum+anisoJ2_count)

   ! loop variables
   integer :: a, c, i, j, t, k
   integer :: offsetJ1, offsetJ2

   !uTerm = 0d0   ! -> F
   !xTerm = 0d0   ! -> G
   !fTerm = 0d0   ! -> H
   !gTerm = 0d0   ! -> X

   !!!write(iul,'(a,g20.10,2i3)') 'DBG:jasicUpdate:start', Rdu%U,Rdu%ieJasold,ie

   if (ie /= Rdu%ieJasold) then
      if (Rdu%ieJasold /= 0) then
        ! necessary to restore one-electron terms for previous electron ieold
        ! also restore U to original U
        Rdu%Fij(Rdu%ieJasold,1:Rdu%ieJasold-1) = Rdu%Fijold(1:Rdu%ieJasold-1)
        Rdu%Fij(Rdu%ieJasold+1:ne,Rdu%ieJasold) = Rdu%Fijold(Rdu%ieJasold+1:ne)
        do k=1,unum+fnum+anisoJ2_count
           Rdu%Fijk(Rdu%ieJasold,1:Rdu%ieJasold-1,k) = Rdu%Fijkold(1:Rdu%ieJasold-1,k)
           Rdu%Fijk(Rdu%ieJasold+1:ne,Rdu%ieJasold,k) = Rdu%Fijkold(Rdu%ieJasold+1:ne,k)
        enddo
        Rdu%Gi(Rdu%ieJasold) = Rdu%Giold
        Rdu%Gki(:,Rdu%ieJasold) = Rdu%Gkiold(:)
        enPowers(1:Rdu%enDim1,Rdu%ieJasold,1:Rdu%enDim2) = Rdu%raibarOld(:,:)
        Rdu%U = Rdu%U0
        Rdu%Uk = Rdu%Uk0
      endif
      Rdu%Fijold(1:ie-1) = Rdu%Fij(ie,1:ie-1)
      Rdu%Fijold(ie+1:ne) = Rdu%Fij(ie+1:ne,ie)
      do k=1,unum+fnum+anisoJ2_count
        Rdu%Fijkold(1:ie-1,k) = Rdu%Fijk(ie,1:ie-1,k)
        Rdu%Fijkold(ie+1:ne,k) = Rdu%Fijk(ie+1:ne,ie,k)
      enddo
      Rdu%Giold = Rdu%Gi(ie)
      Rdu%Gkiold(:) = Rdu%Gki(:,ie)
      Rdu%raibarOld(:,:) = enPowers(1:Rdu%enDim1,ie,1:Rdu%enDim2)
      Rdu%ieJasold = ie
   ! else same electron as before, simply update ie
   endif

   ! subtract current=previous values the three Jastrow terms
   Fsum = 0.d0
   do i=1,ie-1
      Fsum = Fsum + Rdu%Fij(ie,i)
   enddo
   do j=ie+1,ne
      Fsum = Fsum + Rdu%Fij(j,ie)
   enddo

   Fksum = 0.d0
   do i=1,ie-1
      Fksum(:) = Fksum(:) + Rdu%Fijk(ie,i,1:unum+fnum+anisoJ2_count)
   enddo
   do j=ie+1,ne
      Fksum(:) = Fksum(:) + Rdu%Fijk(j,ie,1:unum+fnum+anisoJ2_count)
   enddo


   Rdu%U  = Rdu%U  - Rdu%Gi(ie) - Fsum
   Rdu%Uk(1:unum) = Rdu%Uk(1:unum) - Fksum(1:unum)
   Rdu%Uk(unum+1:unum+xnum) = Rdu%Uk(unum+1:unum+xnum) - Rdu%Gki(1:xnum,ie)
   Rdu%Uk(unum+xnum+1:unum+xnum+fnum) = Rdu%Uk(unum+xnum+1:unum+xnum+fnum) - Fksum(unum+1:unum+fnum)
   !en AO
   Rdu%Uk(unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count) = &
      Rdu%Uk(unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count) - Rdu%Gki(xnum+1:xnum+anisoJ1_count,ie)
   !een AO + eenn AO
   Rdu%Uk(unum+xnum+fnum+anisoJ1_count+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) = &
      Rdu%Uk(unum+xnum+fnum+anisoJ1_count+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) - &
      Fksum(unum+fnum+1:unum+fnum+anisoJ2_count)
   !!!write(iul,'(a,4g20.10)') 'DBG:jasicUpdate1:', Rdu%U,Fsum,Rdu%Gi(ie),Rdu%Xi(ie)

   if (.not.allocated(eePowers)) call abortp("JasicUpdate requires allocated module vars")

   call internal_precalculationOfTerms()   ! ie only !

   ! value for e-e cusp
   eeCusp = 0.5d0
   select case(distType)
   case(DIST_SM, DIST_DOUBLEEXP)
      eeCusp = 0.5d0 / scaleEE
   case(DIST_NEEDS)
      eeCusp = 0.5d0 * scaleEE
   end select

   ! only for the current electron ie !
   i = ie

   ! recalculate the ie terms
   Rdu%Fij(i,1:i-1) = 0.d0
   Rdu%Fij(i+1:ne,i) = 0.d0
   Rdu%Gi(i) = 0.d0
   Rdu%Fijk(i,1:i-1,:) = 0.d0
   Rdu%Fijk(i+1:ne,i,:) = 0.d0
   Rdu%Gki(:,i) = 0.d0


   !---------Electron-Electron-Correlation-Terms-------------------
   call internal_electronElectronTerms()   ! ie only !

   !---------Electron-Nucleus-Correlation-Terms-------------------
   call internal_electronNucleusTerms()    ! ie only !

   !---------Electron-Electron-Nucleus-Correlation-Terms-------------------
   if (fmax > 2) then
      if (useGeneric) then
         call abortp("generic ic update not yet implemented")
      else
         call eenHardcodedUpdateWithUk(Rdu,ie)
      end if
   end if

   !-------- AO term contributions
   offsetJ1 = xnum
   offsetJ2 = unum+fnum
   if (useAOJasTerms) call JasAniso%UpdateOnlyUandUK(Rdu,offsetJ1,offsetJ2,ie)

   !--------- update of U -------------------------
   ! add updated values to U
   Fsum = 0.d0
   do i=1,ie-1
      Fsum = Fsum + Rdu%Fij(ie,i)
   enddo
   do j=ie+1,ne
      Fsum = Fsum + Rdu%Fij(j,ie)
   enddo
   Rdu%U = Rdu%U + Fsum + Rdu%Gi(ie)


   Fksum = 0.d0
   do i=1,ie-1
      !Fksum(:) = Fksum(:) + Rdu%Fijk(ie,i,1:unum+fnum+geennum+geennnum)
      Fksum(:) = Fksum(:) + Rdu%Fijk(ie,i,1:unum+fnum+anisoJ2_count)
   enddo
   do j=ie+1,ne
      !Fksum(:) = Fksum(:) + Rdu%Fijk(j,ie,1:unum+fnum+geennum+geennnum)
      Fksum(:) = Fksum(:) + Rdu%Fijk(j,ie,1:unum+fnum+anisoJ2_count)
   enddo

   Rdu%Uk(1:unum) = Rdu%Uk(1:unum) + Fksum(1:unum)
   Rdu%Uk(unum+1:unum+xnum) = Rdu%Uk(unum+1:unum+xnum) + Rdu%Gki(1:xnum,ie)
   Rdu%Uk(unum+xnum+1:unum+xnum+fnum) = Rdu%Uk(unum+xnum+1:unum+xnum+fnum) + Fksum(unum+1:unum+fnum)
   if (useAOJasTerms) then
     !en AO
     Rdu%Uk(unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count) = &
         Rdu%Uk(unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count) + Rdu%Gki(xnum+1:xnum+anisoJ1_count,ie)
     !two body AO terms
     Rdu%Uk(unum+xnum+fnum+anisoJ1_count+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) = &
              Rdu%Uk(unum+xnum+fnum+anisoJ1_count+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) +&
                    Fksum(unum+fnum+1:unum+fnum+anisoJ2_count)
  endif

   !!!write(iul,'(a,g20.10)') 'DBG:jasicUpdatewpd:U,Uk',Rdu%U
   !!!write(iul,'(5g20.10)') Rdu%Uk

   !!!write(iul,'(a,4g20.10)') 'DBG:jasicUpdate2:', Rdu%U,Fsum,Rdu%Gi(ie),Rdu%Xi(ie)

!    Fsum = 0; Gsum = 0
!    do i=1,ne
!       Gsum = Gsum + Rdu%Gi(i)
!       do j=i+1,ne
!          Fsum = Fsum + Rdu%Fij(j,i)
!          !!!write(iul,'(a,2i3,g20.10)') 'DBG:jasicUpdate3:',i,j,Rdu%Fij(j,i)

!       enddo
!    enddo
!    write(iul,'(a,4g20.10)') 'DBG:jasicUpdate3:', Fsum,Gsum


   call Rdu%markJastrowValid()

contains

   subroutine internal_precalculationOfTerms()

      do a = 1, nclast
         c = atoms(a)%sa
         scale = scaleEN(c)
         power = powerEN(c)

         ! powers of e-n distances
         enPowers(a, i, -2) = 0d0
         enPowers(a, i, -1) = 0d0
         enPowers(a, i, 0) = 1d0

         select case(distType)
         case(DIST_SM)
            r = 1 / (1 + scale * Rdu%rai(a, i))
            enPowers(a, i, 1) = scale * Rdu%rai(a, i) * r

         case(DIST_DOUBLEEXP)
            r = exp(-scale * Rdu%rai(a, i))
            enPowers(a, i, 1) = 1 - r

         case(DIST_NEEDS)
            b = Rdu%rai(a, i) ** power
            r = 1 / (b + scale)
            enPowers(a, i, 1) = Rdu%rai(a, i) * r

         case default
            call abortp("Unknown Jastrow distance type in jasicall")
         end select

         do t = 2, max(xmax, fmax)
            enPowers(a, i, t) = enPowers(a, i, t-1) * enPowers(a, i, 1)
         enddo
      enddo

      ! note: ie in r_ij can be i or j ! Here case j < ie
      do j = 1, i-1
         ! powers of e-e distances
         eePowers(j, i, -2) = 0d0
         eePowers(j, i, -1) = 0d0
         eePowers(j, i, 0) = 1d0

         select case(distType)
         case(DIST_SM)
            r = 1 / (1 + scaleEE * Rdu%rij(j, i))
            eePowers(j, i, 1) = scaleEE * Rdu%rij(j, i) * r

         case(DIST_DOUBLEEXP)
            r = exp(-scaleEE * Rdu%rij(j, i))
            eePowers(j, i, 1) = 1 - r

         case(DIST_NEEDS)
            b = Rdu%rij(j, i) ** powerEE
            r = 1 / (b + scaleEE)
            eePowers(j, i, 1) = Rdu%rij(j, i) * r

         case default
            call abortp("Unknown Jastrow distance type in jasicall")
         end select

         do t = 2, max(umax, fmax)
            eePowers(j, i, t) = eePowers(j, i, t-1) * eePowers(j, i, 1)
         enddo
      enddo
      ! now j > ie
      do j = i + 1, ne
         ! powers of e-e distances
         eePowers(i, j, -2) = 0d0
         eePowers(i, j, -1) = 0d0
         eePowers(i, j, 0) = 1d0

         select case(distType)
         case(DIST_SM)
            r = 1 / (1 + scaleEE * Rdu%rij(i, j))
            eePowers(i, j, 1) = scaleEE * Rdu%rij(i, j) * r

         case(DIST_DOUBLEEXP)
            r = exp(-scaleEE * Rdu%rij(i, j))
            eePowers(i, j, 1) = 1 - r

         case(DIST_NEEDS)
            b = Rdu%rij(i, j) ** powerEE
            r = 1 / (b + scaleEE)
            eePowers(i, j, 1) = Rdu%rij(i, j) * r

         case default
            call abortp("Unknown Jastrow distance type in jasicall")
         end select

         do t = 2, max(umax, fmax)
            eePowers(i, j, t) = eePowers(i, j, t-1) * eePowers(i, j, 1)
         enddo
      enddo

   end subroutine internal_precalculationOfTerms

   subroutine internal_electronElectronTerms()
      do j = 1, i-1
         sumr  = eeCusp * eePowers(j, i, 1)
         if (diffeecusp .and. ((i<=nalpha) .eqv. (j<=nalpha))) then
            sumr  = 0.5d0*eeCusp * eePowers(j, i, 1)
         endif

         ! start at 2 to satisfy cusp
         do t = 2, umax
            sumr  = sumr  + alpha(t) * eePowers(j, i, t)
            Rdu%Fijk(i,j,t-1) = Rdu%Fijk(i,j,t-1) + eePowers(j,i,t)
         enddo
         Rdu%Fij(i,j) = Rdu%Fij(i,j) + sumr
         !uTerm = uTerm + sumr
      enddo
      do j = i+1, ne
         sumr  = eeCusp * eePowers(i, j, 1)
         if (diffeecusp .and. ((i<=nalpha) .eqv. (j<=nalpha))) then
            sumr  = 0.5d0*eeCusp * eePowers(i, j, 1)
         endif

         ! start at 2 to satisfy cusp
         do t = 2, umax
            sumr  = sumr  + alpha(t) * eePowers(i, j, t)
            Rdu%Fijk(j,i,t-1) = Rdu%Fijk(j,i,t-1) + eePowers(i,j,t)
         enddo
         Rdu%Fij(j,i) = Rdu%Fij(j,i) + sumr
         !uTerm = uTerm + sumr
      enddo
   end subroutine internal_electronElectronTerms


   subroutine internal_electronNucleusTerms()
      do a = 1, nclast
         c = atoms(a)%sa
         sumr  = 0d0

         if(nucCusp .and. .not. useAOJasTerms) then
            sumr  = atoms(a)%za * enPowers(a, i, 1)
         endif

         ! start at 2 to satisfy cusp
         do t = 2, xmax
            sumr  = sumr  + beta(t, c) * enPowers(a, i, t)
            k = (c-1)*xpnum + t-1        ! see setNumTerms() below
            Rdu%Gki(k,i) = Rdu%Gki(k,i) + enPowers(a, i, t)
         enddo
         Rdu%Gi(i) = Rdu%Gi(i) + sumr
         !xTerm = xTerm + sumr
      enddo
   end subroutine internal_electronNucleusTerms

end subroutine jasicUpdateWithUk


subroutine eenHardcodedUpdateWithUk(Rdu,ie)
  ! value of f term
   type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
   integer, intent(in) :: ie                  ! electron to update

   real(r8) :: rComb(2)   !!! ???
   real(r8) :: tmp
   real(r8) :: terms(fpnum), termsn(6, fpnum), termsl(2, fpnum)
   real(r8) :: nfterms(fmax), nftermsn(6, fmax), nftermsl(2, fmax)
   integer :: a, c, g, i, j, t, k

   i = ie ! ie only !

   do a = 1, nclast
      c = atoms(a)%sa

      do j = 1, i-1
         ! rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
         ! rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

         if(fmax > 2) then
            ! r_i^2 r_j + r_i r_j^2
            nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                       enPowers(a, i, 1) * enPowers(a, j, 2)

            ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                2 * enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(2) = -0.5d0 * eePowers(j, i, 1) * tmp

            ! r_ij^2 (r_i + r_j)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(j, i, 2) * tmp

            do t = 1, 2
               terms(t) = nfterms(t) - nfterms(3)
            enddo
         endif

         if(fmax > 3) then
            ! 2 ri^2 rj^2
            terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)

            ! rij^2 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(4) = eePowers(j, i, 2) * tmp

            ! ri^3 rj + ri rj^3
            nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 3)

            ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                enPowers(a, i, 2) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(2) = - eePowers(j, i, 1) * tmp

            ! rij^3 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(j, i, 3) * tmp

            ! rij^2 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(4) = eePowers(j, i, 2) * tmp

            do t = 1, 3
               terms(4+t) = nfterms(t) - nfterms(4)
            enddo
         endif

         if(fmax > 4) then
            ! ri^3 rj^2 + ri^2 rj^3
            terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                     enPowers(a, i, 2) * enPowers(a, j, 3)

            ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
            tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                2 * enPowers(a, i, 2) * enPowers(a, j, 2)
            terms(9)  = eePowers(j, i, 1) * tmp

            ! rij^3 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(10) = eePowers(j, i, 3) * tmp

            ! rij^2 (ri^3 + rj^3)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
            terms(11) = eePowers(j, i, 2) * tmp

            ! ri^4 rj + ri rj^4
            nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 4)
            ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
            tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                enPowers(a, i, 3) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 3)
            nfterms(2) = -eePowers(j, i, 1) * tmp

            ! rij^4 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(j, i, 4) * tmp

            ! rij^2 (ri^2 rj + ri rj^2)
            tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(4) = eePowers(j, i, 2) * tmp

            ! rij^3 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(5) = eePowers(j, i, 3) * tmp

            do t = 1, 4
               terms(11+t) = nfterms(t) - nfterms(5)
               termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
               termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
            enddo
         endif

         if(fmax > 5) then
            ! ri^2 rj^4 + ri^4 rj^2
            terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                      enPowers(a, i, 2) * enPowers(a, j, 4)

            ! 2 ri^3 rj^3
            terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

            ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
            tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 3) * enPowers(a, j, 2) - &
                enPowers(a, i, 2) * enPowers(a, j, 3)
            terms(18) = eePowers(j, i, 1) * tmp

            ! rij^4 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(19) = eePowers(j, i, 4) * tmp

            ! 2 rij^2 ri^2 rj^2
            tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
            terms(20) = 2 * eePowers(j, i, 2) * tmp

            ! rij^3 (ri^3 + rj^3)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
            terms(21) = eePowers(j, i, 3) * tmp

            ! rij^2 (ri^4 + rj^4)
            tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
            terms(22) = eePowers(j, i, 2) * tmp

            ! ri^5 rj + ri rj^5
            nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 5)

            ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
            tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 4) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 4)
            nfterms(2) = - eePowers(j, i, 1) * tmp

            ! rij^5 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(j, i, 5) * tmp

            ! rij^2 (ri rj^3 + ri^3 rj)
            tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 3)
            nfterms(4) = eePowers(j, i, 2) * tmp

            ! rij^3 (ri rj^2 + ri^2 rj)
            tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(5) = eePowers(j, i, 3) * tmp

            ! rij^4 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(6) = eePowers(j, i, 4) * tmp

            do t = 1, 5
               terms(22+t) = nfterms(t) - nfterms(6)
            enddo
         endif

         do g = 1, fpnum
            !fTerm = fTerm + gamma(g, c) * terms(g)
            Rdu%Fij(i,j) = Rdu%Fij(i,j) + gamma(g, c) * terms(g)
            k = unum + (c-1)*fpnum + g        ! see setNumTerms() below
            !!!write(iul,'(a,5i4,g20.10)') 'DBG:up1eenhcwuk:',k,unum,ncdiff,c,g,terms(g)
            Rdu%Fijk(i,j,k) = Rdu%Fijk(i,j,k) + terms(g)
         enddo
      enddo  ! j
      do j = i+1, ne
         ! rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
         ! rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

         if(fmax > 2) then
            ! r_i^2 r_j + r_i r_j^2
            nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                       enPowers(a, i, 1) * enPowers(a, j, 2)

            ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                2 * enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(2) = -0.5d0 * eePowers(i, j, 1) * tmp

            ! r_ij^2 (r_i + r_j)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(i, j, 2) * tmp

            do t = 1, 2
               terms(t) = nfterms(t) - nfterms(3)
            enddo
         endif

         if(fmax > 3) then
            ! 2 ri^2 rj^2
            terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)

            ! rij^2 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(4) = eePowers(i, j, 2) * tmp

            ! ri^3 rj + ri rj^3
            nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 3)

            ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                enPowers(a, i, 2) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(2) = - eePowers(i, j, 1) * tmp

            ! rij^3 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(i, j, 3) * tmp

            ! rij^2 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(4) = eePowers(i, j, 2) * tmp

            do t = 1, 3
               terms(4+t) = nfterms(t) - nfterms(4)
            enddo
         endif

         if(fmax > 4) then
            ! ri^3 rj^2 + ri^2 rj^3
            terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                     enPowers(a, i, 2) * enPowers(a, j, 3)

            ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
            tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                2 * enPowers(a, i, 2) * enPowers(a, j, 2)
            terms(9)  = eePowers(i, j, 1) * tmp

            ! rij^3 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(10) = eePowers(i, j, 3) * tmp

            ! rij^2 (ri^3 + rj^3)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
            terms(11) = eePowers(i, j, 2) * tmp

            ! ri^4 rj + ri rj^4
            nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 4)
            ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
            tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                enPowers(a, i, 3) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 3)
            nfterms(2) = -eePowers(i, j, 1) * tmp

            ! rij^4 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(i, j, 4) * tmp

            ! rij^2 (ri^2 rj + ri rj^2)
            tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(4) = eePowers(i, j, 2) * tmp

            ! rij^3 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(5) = eePowers(i, j, 3) * tmp

            do t = 1, 4
               terms(11+t) = nfterms(t) - nfterms(5)
               termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
               termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
            enddo
         endif

         if(fmax > 5) then
            ! ri^2 rj^4 + ri^4 rj^2
            terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                      enPowers(a, i, 2) * enPowers(a, j, 4)

            ! 2 ri^3 rj^3
            terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

            ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
            tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 3) * enPowers(a, j, 2) - &
                enPowers(a, i, 2) * enPowers(a, j, 3)
            terms(18) = eePowers(i, j, 1) * tmp

            ! rij^4 (ri^2 + rj^2)
            tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
            terms(19) = eePowers(i, j, 4) * tmp

            ! 2 rij^2 ri^2 rj^2
            tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
            terms(20) = 2 * eePowers(i, j, 2) * tmp

            ! rij^3 (ri^3 + rj^3)
            tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
            terms(21) = eePowers(i, j, 3) * tmp

            ! rij^2 (ri^4 + rj^4)
            tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
            terms(22) = eePowers(i, j, 2) * tmp

            ! ri^5 rj + ri rj^5
            nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 5)

            ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
            tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 4) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 4)
            nfterms(2) = - eePowers(i, j, 1) * tmp

            ! rij^5 (ri + rj)
            tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
            nfterms(3) = eePowers(i, j, 5) * tmp

            ! rij^2 (ri rj^3 + ri^3 rj)
            tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 3)
            nfterms(4) = eePowers(i, j, 2) * tmp

            ! rij^3 (ri rj^2 + ri^2 rj)
            tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
            nfterms(5) = eePowers(i, j, 3) * tmp

            ! rij^4 ri rj
            tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
            nfterms(6) = eePowers(i, j, 4) * tmp

            do t = 1, 5
               terms(22+t) = nfterms(t) - nfterms(6)
            enddo
         endif

         do g = 1, fpnum
            !fTerm = fTerm + gamma(g, c) * terms(g)
            Rdu%Fij(j,i) = Rdu%Fij(j,i) + gamma(g, c) * terms(g)
            k = unum + (c-1)*fpnum + g        ! see setNumTerms() below
            !!!write(iul,'(a,5i4,g20.10)') 'DBG:up2eenhcwuk:',k,unum,ncdiff,c,g,terms(g)
            Rdu%Fijk(j,i,k) = Rdu%Fijk(j,i,k) + terms(g)
         enddo
      enddo  ! j
   enddo  ! a

end subroutine eenHardcodedUpdateWithUk





!=======================================================================

subroutine numDerivs(x, y, z, rai, rij, optType, optMode, nn)
  !use wfParameters_m

  real(r8), intent(in)            :: x(:), y(:), z(:)
  real(r8), intent(in)            :: rai(:, :), rij(:, :)
  integer, intent(in), optional :: nn
  character(len=*),intent(in)   :: optType
  integer, intent(in)           :: optMode

  integer :: offset, np, stat, i, opt, np1, np2, npnl
  real(r8), allocatable :: up(:), upgrad(:, :), uplapl(:), uplapli(:, :)
  integer, parameter :: points = 7
  logical :: resetOptMode

  if(optMode == OPT_ALL) then
    ! only generate numerical derivatives for nonlinear params
    opt = OPT_NONLIN
    resetOptMode = .true.
    offset = numParams
  else
    opt = optMode
    resetOptMode = .false.
    offset = 0
  endif

  call abortp("jastrowIC_m: numDerivs currently disabled")

  ! since optMode is being passed as intent(in), it's actually a pointer
  ! to the originally passed variable - which is curOptMode
  ! since getVectorLenIC overwrites curOptMode, the local optMode variable
  ! here is also overwritten. making a copy doesn't help, as this is "optimized"
  ! away from the compiler
  call getVectorLenIC(opt,np1,np2,npnl)
  np = np1+np2+npnl

  allocate(up(np), upgrad(3*ne, np), uplapl(np), uplapli(ne, np), stat=stat)
  call assert(stat == 0, "jastrowIC_m: numDerivs allocation failed")

  !!!call wfparams_numDerivs(0, x, y, z, rai, rij, &
  !!!                        up, upgrad, uplapl, uplapli, &
  !!!                        optType, optMode, points, nn)

  uk(offset + 1:offset + np) = up(1:np)
  uklapl(offset + 1:offset + np) = uplapl(1:np)
  do i = 1, 3*ne
    ukgrad(i, offset + 1:offset + np) = upgrad(i, 1:np)
  enddo
  do i = 1, ne
    uklapli(i, offset + 1:offset + np) = uplapli(i, 1:np)
  enddo

  deallocate(up, upgrad, uplapl, uplapli)

  ! reset saved optMode to original one
  ! WARNING: DO NOT USE optMode, as it's being overwritten by the call to
  ! getVectorLen, see above
  if(resetOptMode) curOptMode = OPT_ALL
end subroutine numDerivs

!=======================================================================


subroutine jasicall(x, y, z, rai, rij, optType, ju, jud, julapl, julapli, uuk, withDerivs, nn)
!---------------------------------------------------------------------------------------!
! main calculations for all electrons moved

  real(r8), intent(in)  ::  x(:), y(:), z(:)     ! cartesian coords
  real(r8), intent(in)  ::  rai(:, :), rij(:, :) ! distances r_ai, r_ij
  character(len=*)    ::  optType              ! parameter optimization
  real(r8), intent(out) ::  ju                   ! U
  real(r8), intent(out) ::  jud(:)               ! nabla U
  real(r8), intent(out) ::  julapl               ! sum of laplacians
  real(r8), intent(out) ::  julapli(:)           ! laplacian U
  integer, intent(in), optional :: withDerivs  ! 0 for no drv 1 for all and 2 for only uk(used in one electron move in ECP)
  integer, intent(in), optional ::  nn                   ! current index of nElecConfigs
  real(r8), intent(out), optional ::  uuk(:)

  ! values of correlation terms
  real(r8) :: uTerm, xTerm, fTerm, gTerm
  ! nabla of correlation terms
  real(r8) :: uDeriv(3*ne), xDeriv(3*ne), fDeriv(3*ne), gDeriv(3*ne)
  ! laplacian
  real(r8) :: uLapli(ne), xLapli(ne), fLapli(ne), gLapli(ne)

  ! jastrow derivatives
  integer :: calcDerivs ! 0 for no drv 1 for all and 2 for only uk(used in one electron move in ECP)
  ! parameter derivatives
  logical :: deriveLinParams
  logical :: deriveNumericalParams

  ! x/y/z distances e-n and e-e
  real(r8) :: xai(nclast, ne), yai(nclast, ne), zai(nclast, ne)
  real(r8) :: xij(ne, ne), yij(ne, ne), zij(ne, ne)

  ! sum in x, u and f-term, derivative, nabla
  real(r8) :: sumr, sumrd, sumrn
  ! position of parameter for param derivatives
  integer :: pos
  ! parameter value to satisfy e-e cusp condition
  real(r8) :: eeCusp

  ! commonly used terms
  real(r8) :: b, r, s, v, tmpn(3), tmpl, scale, power

  ! loop variables
  integer :: a, c, i, j, t, m, nntemp

  !aniso
  real(r8), allocatable :: gk(:),gkgrad(:,:),gklapl(:),gklapli(:,:)

  call initArrays()        ! inner subroutine: initialization of arrays

  call precalculationOfTerms()  ! inner subroutine

  ! value for e-e cusp
  eeCusp = 0.5d0
  select case(distType)
  case(DIST_SM, DIST_DOUBLEEXP)
    eeCusp = 0.5d0 / scaleEE
  case(DIST_NEEDS)
    eeCusp = 0.5d0 * scaleEE
  end select
  !!!if (useAOJasTerms) eeCusp = 0.d0


  if (calcDerivs == 1) then
    call electronElectronTermsAll()
    call electronNucleusTermsAll()
  else
    call electronElectronTerms()
    call electronNucleusTerms()
  endif


  ! what happens if this is all in a geminal loop?
  !??? Very strange, why here uklapl for fAOIdx%nums ???
!   if (deriveLinParams) then
!     begin = unum + xnum
!     end   = unum + xnum + fnum
!     if (useAOJasTerms) then
!       begin = begin + fnum! + 1
!       end   = end + gnum - 1
!     end if
!     do m = begin, end
!       uklapl(m) = sum(uklapli(:,m))
!     end do
!     if (useAOJasTerms) then
!       begin = begin + gnum-1
!       end = end + (gnum-1)*(gnum-1)
!       do m = begin, end
!         uklapl(m) = sum(uklapli(:,m))
!       end do
!     end if
!   end if


  !---------Electron-Electron-Nucleus-Correlation-Terms-------------------
  if (fmax > 2) then
    if (useGeneric) then
      call eenGeneric(fTerm, fDeriv, fLapli, deriveLinParams, calcDerivs)
    else
      if (calcDerivs == 1) then
        call eenHardcodedAll(fTerm, fDeriv, fLapli, deriveLinParams)
      elseif(calcDerivs == 0) then
        call eenHardcodedNoDerivs(fTerm)
      else
        call eenHardcodedWithUk(fTerm,deriveLinParams)
      end if
    end if
  end if


  if (.not.(present(nn))) then
    nntemp=1
  else
    nntemp = nn
  endif
  !  AO terms contributions
  if (useAOJasTerms) then
    if(deriveLinParams) then
      if (allocated(gk)) then
        if (.not.(size(gk)==anisoJ1_count+anisoJ2_count)) then
          deallocate(gk,gkgrad,gklapli,gklapl)
          allocate(gk(anisoJ1_count+anisoJ2_count), gkgrad(3*ne,anisoJ1_count+anisoJ2_count), &
            gklapli(ne,anisoJ1_count+anisoJ2_count),gklapl(anisoJ1_count+anisoJ2_count))
        endif
      else
        allocate(gk(anisoJ1_count+anisoJ2_count), gkgrad(3*ne,anisoJ1_count+anisoJ2_count), &
          gklapli(ne,anisoJ1_count+anisoJ2_count),gklapl(anisoJ1_count+anisoJ2_count))
      endif
      call JasAniso%eval(gTerm,gDeriv,gLapli,deriveLinParams,nntemp,gk,gkgrad,gklapli,gklapl)
    else
      call JasAniso%eval(gTerm,gDeriv,gLapli,deriveLinParams,nntemp)
    endif
    if(deriveLinParams) then
      uk(unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) = gk(:)
      ukgrad(:,unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) = gkgrad(:,:)
      uklapli(:,unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) = gklapli(:,:)
      uklapl(unum+xnum+fnum+1:unum+xnum+fnum+anisoJ1_count+anisoJ2_count) = gklapl(:)
    endif
  endif

  !---------Calculation of U and its Deriviatives-------------------------
  ju = uTerm + xTerm + fTerm + gTerm

  !!!write(iul,'(a,6g20.10)') 'DBG:jasicall:',ju,uTerm,xTerm,fTerm,uTerm+fTerm,gTerm

  if (calcDerivs == 1) then
    jud = uDeriv(1:3*ne) + xDeriv(1:3*ne) + fDeriv(1:3*ne) + gDeriv(1:3*ne)
    julapli = uLapli + xLapli + fLapli + gLapli
    julapl = sum(julapli(1:ne))
  endif

  !---------Nonlinear parameter derivatives-------------------------------
  if(deriveNumericalParams) then
    if(useAOJasTerms) then
      call numDerivs(x, y, z, rai, rij, optType, curOptMode, nn)
    else
      call numDerivs(x, y, z, rai, rij, optType, curOptMode)
    endif
  endif

  if (present(uuk) .and. allocated(uk)) then
      uuk(:)=uk(:)
      !!!write(iul,'(a,g20.10)') 'DBG:jasicall wpd:u,uk',ju
      !!!write(iul,'(5g20.10)') uuk
  endif
contains

  subroutine initArrays()
  !---------------------!
  ! initialization of arrays
  ju = 0d0
  jud = 0d0
  julapl = 0d0
  julapli = 0d0

  uTerm = 0d0
  xTerm = 0d0
  fTerm = 0d0
  gTerm = 0d0
  uDeriv = 0d0
  xDeriv = 0d0
  fDeriv = 0d0
  gDeriv = 0d0
  uLapli = 0d0
  xLapli = 0d0
  fLapli = 0d0
  gLapli = 0d0

  if(.not. allocatedTerms) then
    call allocateTerms()
  endif

  calcDerivs = 1
  if(present(withDerivs)) calcDerivs = withDerivs

  deriveLinParams = .false.
  deriveNumericalParams = .false.
  if(optType == "jastrow" .or. optType == "jas+ci" .or.  optType == "jas+mo" &
                                             .or. optType == "jas+mo+ci" ) then
    call assert((calcDerivs /= 0), "jasicall: cannot calculate parameter derivatives without " // &
      "jastrow derivatives")

    if(curOptMode == OPT_LIN .or. curOptMode == OPT_ALL) then
      deriveLinParams = .true.
    endif

    if(curOptMode == OPT_NUM_LIN .or. &
       curOptMode == OPT_NONLIN  .or. &
       curOptMode == OPT_ALL     .or. &
       curOptMode == OPT_NUM_ALL) then
      deriveNumericalParams = .true.
    endif

    call assert(allocated(uk), "jastrowParamData not allocated")
    uk = 0d0
    if (calcDerivs /= 2) then
    ukgrad = 0d0
    uklapl = 0d0
    uklapli = 0d0
  endif
  endif

  end subroutine initArrays

  subroutine precalculationOfTerms()

  do i = 1, ne
    do a = 1, nclast
      c = atoms(a)%sa
      scale = scaleEN(c)
      power = powerEN(c)

      ! x/y/z distances for e-n

      if (calcDerivs == 1) then
        xai(a, i) = x(i) - atoms(a)%cx
        yai(a, i) = y(i) - atoms(a)%cy
        zai(a, i) = z(i) - atoms(a)%cz
      endif

      ! powers of e-n distances
      enPowers(a, i, -2) = 0d0
      enPowers(a, i, -1) = 0d0
      enPowers(a, i, 0) = 1d0

      select case(distType)
      case(DIST_SM)
        r = 1 / (1 + scale * rai(a, i))
        enPowers(a, i, 1) = scale * rai(a, i) * r

        if (calcDerivs == 1) then
          v = scale * r * r
          s = v / rai(a, i)
          raiDeriv(1, a, i) = xai(a, i) * s
          raiDeriv(2, a, i) = yai(a, i) * s
          raiDeriv(3, a, i) = zai(a, i) * s

          !raiSquare(a, i) = raiDeriv(a, i, 1) ** 2 + raiDeriv(a, i, 2) ** 2 + raiDeriv(a, i, 3) ** 2
          raiSquare(a, i) = v * v
          raiLapl(a, i) = 2 * s * r
        endif
      case(DIST_DOUBLEEXP)
        r = exp(-scale * rai(a, i))
        enPowers(a, i, 1) = 1 - r

        if (calcDerivs == 1) then
          v = scale * r
          s = v / rai(a, i)
          raiDeriv(1, a, i) = xai(a, i) * s
          raiDeriv(2, a, i) = yai(a, i) * s
          raiDeriv(3, a, i) = zai(a, i) * s

          raiSquare(a, i) = v * v
          raiLapl(a, i) = s * (2 - scale * rai(a, i))
        endif
      case(DIST_NEEDS)
        b = rai(a, i) ** power
        r = 1 / (b + scale)
        enPowers(a, i, 1) = rai(a, i) * r

        if (calcDerivs == 1) then
          v = (scale - (power - 1) * b) * r * r
          s = v / rai(a, i)
          raiDeriv(1, a, i) = xai(a, i) * s
          raiDeriv(2, a, i) = yai(a, i) * s
          raiDeriv(3, a, i) = zai(a, i) * s

          raiSquare(a, i) = v * v
          raiLapl(a, i) = (2 * scale * scale - (power - 1) * (power + 4) * scale * b + &
            (power - 1) * (power - 2) * b * b) * r * r * r/ rai(a, i)
        endif
      case default
        call abortp("Unknown Jastrow distance type in jasicall")
      end select

      do t = 2, max(xmax, fmax)
        enPowers(a, i, t) = enPowers(a, i, t-1) * enPowers(a, i, 1)
      enddo
    enddo
    do j = i + 1, ne
      if (calcDerivs == 1) then
        ! x/y/z distance for e-e
        xij(i, j) = x(i) - x(j)
        yij(i, j) = y(i) - y(j)
        zij(i, j) = z(i) - z(j)
      endif

      ! powers of e-e distances
      eePowers(i, j, -2) = 0d0
      eePowers(i, j, -1) = 0d0
      eePowers(i, j, 0) = 1d0

      select case(distType)
      case(DIST_SM)
        r = 1 / (1 + scaleEE * rij(i, j))
        eePowers(i, j, 1) = scaleEE * rij(i, j) * r

        if (calcDerivs == 1) then
          v = scaleEE * r * r
          s = v / rij(i, j)
          rijDeriv(1, i, j) = xij(i, j) * s
          rijDeriv(2, i, j) = yij(i, j) * s
          rijDeriv(3, i, j) = zij(i, j) * s

          !rijSquare(i, j) = rijDeriv(i, j, 1) ** 2 + rijDeriv(i, j, 2) ** 2 + rijDeriv(i, j, 3) ** 2
          rijSquare(i, j) = v * v
          rijLapl(i, j) = 2 * s * r
        endif
      case(DIST_DOUBLEEXP)
        r = exp(-scaleEE * rij(i, j))
        eePowers(i, j, 1) = 1 - r

        if (calcDerivs == 1) then
          v = scaleEE * r
          s = v / rij(i, j)
          rijDeriv(1, i, j) = xij(i, j) * s
          rijDeriv(2, i, j) = yij(i, j) * s
          rijDeriv(3, i, j) = zij(i, j) * s

          rijSquare(i, j) = v * v
          rijLapl(i, j) = s * (2 - scaleEE * rij(i, j))
        endif
      case(DIST_NEEDS)
        b = rij(i, j) ** powerEE
        r = 1 / (b + scaleEE)
        eePowers(i, j, 1) = rij(i, j) * r

        if (calcDerivs == 1) then
          v = (scaleEE - (powerEE - 1) * b) * r * r
          s = v / rij(i, j)
          rijDeriv(1, i, j) = xij(i, j) * s
          rijDeriv(2, i, j) = yij(i, j) * s
          rijDeriv(3, i, j) = zij(i, j) * s

          rijSquare(i, j) = v * v
          rijLapl(i, j) = (2 * scaleEE * scaleEE - (powerEE - 1) * (powerEE + 4) * scaleEE * b + &
            (powerEE - 1) * (powerEE - 2) * b * b) * r * r * r/ rij(i, j)
        endif
      case default
        call abortp("Unknown Jastrow distance type in jasicall")
      end select

      do t = 2, max(umax, fmax)
        eePowers(i, j, t) = eePowers(i, j, t-1) * eePowers(i, j, 1)
      enddo
    enddo
  enddo

  end subroutine precalculationOfTerms

  ! subroutine precalculationOfTerms()

  ! do i = 1, ne
  !   do a = 1, nclast
  !     c = atoms(a)%sa
  !     scale = scaleEN(c)
  !     power = powerEN(c)

  !     ! x/y/z distances for e-n

  !     ! if (calcDerivs /= 0) then
  !     !   xai(a, i) = x(i) - atoms(a)%cx
  !     !   yai(a, i) = y(i) - atoms(a)%cy
  !     !   zai(a, i) = z(i) - atoms(a)%cz
  !     ! endif

  !     ! powers of e-n distances
  !     enPowers(a, i, -2) = 0d0
  !     enPowers(a, i, -1) = 0d0
  !     enPowers(a, i, 0) = 1d0

  !     select case(distType)
  !     case(DIST_SM)
  !       r = 1 / (1 + scale * rai(a, i))
  !       enPowers(a, i, 1) = scale * rai(a, i) * r

  !       ! if (calcDerivs /= 0) then
  !       !   v = scale * r * r
  !       !   s = v / rai(a, i)
  !       !   raiDeriv(1, a, i) = xai(a, i) * s
  !       !   raiDeriv(2, a, i) = yai(a, i) * s
  !       !   raiDeriv(3, a, i) = zai(a, i) * s

  !       !   !raiSquare(a, i) = raiDeriv(a, i, 1) ** 2 + raiDeriv(a, i, 2) ** 2 + raiDeriv(a, i, 3) ** 2
  !       !   raiSquare(a, i) = v * v
  !       !   raiLapl(a, i) = 2 * s * r
  !       ! endif
  !     case(DIST_DOUBLEEXP)
  !       r = exp(-scale * rai(a, i))
  !       enPowers(a, i, 1) = 1 - r

  !       ! if (calcDerivs /= 0) then
  !       !   v = scale * r
  !       !   s = v / rai(a, i)
  !       !   raiDeriv(1, a, i) = xai(a, i) * s
  !       !   raiDeriv(2, a, i) = yai(a, i) * s
  !       !   raiDeriv(3, a, i) = zai(a, i) * s

  !       !   raiSquare(a, i) = v * v
  !       !   raiLapl(a, i) = s * (2 - scale * rai(a, i))
  !       ! endif
  !     case(DIST_NEEDS)
  !       b = rai(a, i) ** power
  !       r = 1 / (b + scale)
  !       enPowers(a, i, 1) = rai(a, i) * r

  !       ! if (calcDerivs /= 0) then
  !       !   v = (scale - (power - 1) * b) * r * r
  !       !   s = v / rai(a, i)
  !       !   raiDeriv(1, a, i) = xai(a, i) * s
  !       !   raiDeriv(2, a, i) = yai(a, i) * s
  !       !   raiDeriv(3, a, i) = zai(a, i) * s

  !       !   raiSquare(a, i) = v * v
  !       !   raiLapl(a, i) = (2 * scale * scale - (power - 1) * (power + 4) * scale * b + &
  !       !     (power - 1) * (power - 2) * b * b) * r * r * r/ rai(a, i)
  !       ! endif
  !     case default
  !       call abortp("Unknown Jastrow distance type in jasicall")
  !     end select

  !     do t = 2, max(xmax, fmax)
  !       enPowers(a, i, t) = enPowers(a, i, t-1) * enPowers(a, i, 1)
  !     enddo
  !   enddo
  !   do j = i + 1, ne
  !     ! if (calcDerivs /= 0) then
  !     !   ! x/y/z distance for e-e
  !     !   xij(i, j) = x(i) - x(j)
  !     !   yij(i, j) = y(i) - y(j)
  !     !   zij(i, j) = z(i) - z(j)
  !     ! endif

  !     ! powers of e-e distances
  !     eePowers(i, j, -2) = 0d0
  !     eePowers(i, j, -1) = 0d0
  !     eePowers(i, j, 0) = 1d0

  !     select case(distType)
  !     case(DIST_SM)
  !       r = 1 / (1 + scaleEE * rij(i, j))
  !       eePowers(i, j, 1) = scaleEE * rij(i, j) * r

  !       ! if (calcDerivs /= 0) then
  !       !   v = scaleEE * r * r
  !       !   s = v / rij(i, j)
  !       !   rijDeriv(1, i, j) = xij(i, j) * s
  !       !   rijDeriv(2, i, j) = yij(i, j) * s
  !       !   rijDeriv(3, i, j) = zij(i, j) * s

  !       !   !rijSquare(i, j) = rijDeriv(i, j, 1) ** 2 + rijDeriv(i, j, 2) ** 2 + rijDeriv(i, j, 3) ** 2
  !       !   rijSquare(i, j) = v * v
  !       !   rijLapl(i, j) = 2 * s * r
  !       ! endif
  !     case(DIST_DOUBLEEXP)
  !       r = exp(-scaleEE * rij(i, j))
  !       eePowers(i, j, 1) = 1 - r

  !       ! if (calcDerivs /= 0) then
  !       !   v = scaleEE * r
  !       !   s = v / rij(i, j)
  !       !   rijDeriv(1, i, j) = xij(i, j) * s
  !       !   rijDeriv(2, i, j) = yij(i, j) * s
  !       !   rijDeriv(3, i, j) = zij(i, j) * s

  !       !   rijSquare(i, j) = v * v
  !       !   rijLapl(i, j) = s * (2 - scaleEE * rij(i, j))
  !       ! endif
  !     case(DIST_NEEDS)
  !       b = rij(i, j) ** powerEE
  !       r = 1 / (b + scaleEE)
  !       eePowers(i, j, 1) = rij(i, j) * r

  !       ! if (calcDerivs /= 0) then
  !       !   v = (scaleEE - (powerEE - 1) * b) * r * r
  !       !   s = v / rij(i, j)
  !       !   rijDeriv(1, i, j) = xij(i, j) * s
  !       !   rijDeriv(2, i, j) = yij(i, j) * s
  !       !   rijDeriv(3, i, j) = zij(i, j) * s

  !       !   rijSquare(i, j) = v * v
  !       !   rijLapl(i, j) = (2 * scaleEE * scaleEE - (powerEE - 1) * (powerEE + 4) * scaleEE * b + &
  !       !     (powerEE - 1) * (powerEE - 2) * b * b) * r * r * r/ rij(i, j)
  !       ! endif
  !     case default
  !       call abortp("Unknown Jastrow distance type in jasicall")
  !     end select

  !     do t = 2, max(umax, fmax)
  !       eePowers(i, j, t) = eePowers(i, j, t-1) * eePowers(i, j, 1)
  !     enddo
  !   enddo
  ! enddo

  ! end subroutine precalculationOfTerms

  subroutine electronElectronTermsAll()
  do i = 1, ne
    do j = i + 1, ne

      sumr  = eeCusp * eePowers(i, j, 1)
      sumrd = eeCusp
      sumrn = 0d0
      if (diffeecusp .and. ((i<=nalpha) .eqv. (j<=nalpha))) then
          sumr  = 0.5d0*eeCusp * eePowers(i, j, 1)
          sumrd = 0.5d0*eeCusp
          sumrn = 0d0
      endif

      ! start at 2 to satisfy cusp
      do t = 2, umax
        sumr  = sumr  + alpha(t) * eePowers(i, j, t)
      enddo
      uTerm = uTerm + sumr
      do t = 2, umax
        sumrd = sumrd + alpha(t) * t * eePowers(i, j, t-1)
        sumrn = sumrn + alpha(t) * t * (t-1) * eePowers(i, j, t-2)
      enddo

      tmpn(1:3) = sumrd * rijDeriv(1:3, i, j)
      uDeriv(3*i-2:3*i) = uDeriv(3*i-2:3*i) + tmpn(1:3)
      uDeriv(3*j-2:3*j) = uDeriv(3*j-2:3*j) - tmpn(1:3)

      tmpl = sumrd * rijLapl(i, j) + sumrn * rijSquare(i, j)
      uLapli(i) = uLapli(i) + tmpl
      uLapli(j) = uLapli(j) + tmpl

      if(deriveLinParams) then
        do t = 2, umax
          pos = t - 1
          uk(pos) = uk(pos) + eePowers(i, j, t)

          tmpn(1:3) = t * eePowers(i, j, t-1) * rijDeriv(1:3, i, j)
          ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + tmpn(1:3)
          ukgrad(3*j-2:3*j, pos) = ukgrad(3*j-2:3*j, pos) - tmpn(1:3)

          tmpl  = t *         eePowers(i, j, t-1) * rijLapl(i, j) + &
                  t * (t-1) * eePowers(i, j, t-2) * rijSquare(i, j)
          uklapli(i, pos) = uklapli(i, pos) + tmpl
          uklapli(j, pos) = uklapli(j, pos) + tmpl
          uklapl(pos) = uklapl(pos) + 2 * tmpl
        enddo
      endif
    enddo
  enddo

  ! write(iul, *) "e-e corr complete"
  end subroutine electronElectronTermsAll

  subroutine electronElectronTerms()
  do i = 1, ne
    do j = i + 1, ne

      sumr  = eeCusp * eePowers(i, j, 1)
      sumrd = eeCusp
      sumrn = 0d0
      if (diffeecusp .and. ((i<=nalpha) .eqv. (j<=nalpha))) then
          sumr  = 0.5d0*eeCusp * eePowers(i, j, 1)
          sumrd = 0.5d0*eeCusp
          sumrn = 0d0
      endif

      ! start at 2 to satisfy cusp
      do t = 2, umax
        sumr  = sumr  + alpha(t) * eePowers(i, j, t)
      enddo
      uTerm = uTerm + sumr

      if ( calcDerivs == 0) cycle
      if(deriveLinParams) then
        do t = 2, umax
          pos = t - 1
          uk(pos) = uk(pos) + eePowers(i, j, t)
        enddo
      endif
    enddo
  enddo

  ! write(iul, *) "e-e corr complete"
  end subroutine electronElectronTerms

  subroutine electronNucleusTermsAll()
  do a = 1, nclast
    c = atoms(a)%sa
    do i = 1, ne

      sumr  = 0d0
      sumrd = 0d0
      sumrn = 0d0

      if(nucCusp .and. .not. useAOJasTerms) then
        sumr  = atoms(a)%za * enPowers(a, i, 1)
        sumrd = atoms(a)%za
      endif

      ! start at 2 to satisfy cusp
      do t = 2, xmax
        sumr  = sumr  + beta(t, c) * enPowers(a, i, t)
      enddo
      xTerm = xTerm + sumr

      do t = 2, xmax
        sumrd = sumrd + beta(t, c) * t * enPowers(a, i, t-1)
        sumrn = sumrn + beta(t, c) * t * (t-1) * enPowers(a, i, t-2)
      enddo

      xDeriv(3*i-2:3*i) = xDeriv(3*i-2:3*i) + sumrd * raiDeriv(1:3, a, i)
      xLapli(i) = xLapli(i) + sumrd * raiLapl(a, i) + sumrn * raiSquare(a, i)

      if(deriveLinParams) then
        do t = 2, xmax
          pos = unum + (c - 1) * xpnum + t - 1

          uk(pos) = uk(pos) + enPowers(a, i, t)
          ukgrad(3*i-2:3*i, pos) = ukgrad(3*i-2:3*i, pos) + &
                                   t * enPowers(a, i, t-1) * raiDeriv(1:3, a, i)

          tmpl = t * enPowers(a, i, t-1) * raiLapl(a, i) + t * (t-1) * enPowers(a, i, t-2) * raiSquare(a, i)
          uklapli(i, pos) = uklapli(i, pos) + tmpl
          uklapl(pos) = uklapl(pos) + tmpl
        enddo
      endif
    enddo
  enddo
  end subroutine electronNucleusTermsAll

  subroutine electronNucleusTerms()
  do a = 1, nclast
    c = atoms(a)%sa
    do i = 1, ne

      sumr  = 0d0
      sumrd = 0d0
      sumrn = 0d0

      if(nucCusp .and. .not. useAOJasTerms) then
        sumr  = atoms(a)%za * enPowers(a, i, 1)
        sumrd = atoms(a)%za
      endif

      ! start at 2 to satisfy cusp
      do t = 2, xmax
        sumr  = sumr  + beta(t, c) * enPowers(a, i, t)
      enddo
      xTerm = xTerm + sumr

      if ( calcDerivs == 0) cycle

      if(deriveLinParams) then
        do t = 2, xmax
          pos = unum + (c - 1) * xpnum + t - 1
          uk(pos) = uk(pos) + enPowers(a, i, t)
        enddo
      endif
    enddo
  enddo
  end subroutine electronNucleusTerms

end subroutine jasicall

subroutine eenGeneric(fTerm, fDeriv, fLapli, deriveLinParams, calcDerivs)
  ! value of f term + derivs
  real(r8), intent(inout) :: fTerm, fDeriv(:), fLapli(:)
  logical, intent(in) :: deriveLinParams
  integer, intent(in) :: calcDerivs
  ! param derivs of f term
  real(r8) :: fkTerm(fnum+1), fkDeriv(3*ne, fnum+1), fkLapli(ne, fnum+1)
  ! generated non-free e-e-n terms + derivatives
  real(r8) :: nfTerms(fmax), nfDerivs(fmax, 6), nfLapl(fmax, 2)

  real(r8) :: rComb(2)
  real(r8) :: coeff

  real(r8) :: tmp, tmpd, tmpn(3), tmpl
  integer :: a, c, i, j, k, l, m, p, t
  !integer :: x, gt, params
  integer :: g, gstart, gend

  fkTerm = 0d0
  fkDeriv = 0d0
  fkLapli = 0d0

  do a = 1, nclast
    c = atoms(a)%sa

    do i = 1, ne
      do j = i + 1, ne

        rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
        rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

        ! this variable describes the paramter index used by the current term
        !params = 0
        g = 0

        do p = 3, fmax
          !params = (30*(p-1)*(p-1) - 124*(p-1) - 99  + 3*(IAND(p,1)*2 - 1) + 4*((p-1.0)**3))/48 + 5
          !gstart = params;
          gstart = g + 1
          ! free terms without r_ij
          do k = 2, p/2
            g = g + 1
            l = p - k
            !g = params + k - 2
            ! n.b. duplicates terms for k = p/2
            !     e.g. r_i^2*r_j^2+r_i^2*r_j^2
            !     shouldn't matter b/c of coefficient
            tmp = enPowers(a, i, k) * enPowers(a, j, l) + &
                  enPowers(a, i, l) * enPowers(a, j, k)
            fTerm = fTerm + gamma(g, c) * tmp

            if(calcDerivs == 0) cycle

            fkTerm(g) = tmp

            tmpd = k * enPowers(a, i, k - 1) * enPowers(a, j, l) + &
                   l * enPowers(a, i, l - 1) * enPowers(a, j, k)
            tmpn(:) = raiDeriv(:, a, i) * tmpd
            fDeriv(3*i-2:3*i) = fDeriv(3*i-2:3*i) + gamma(g, c) * tmpn(:)
            fkDeriv(3*i-2:3*i, g) = tmpn(:)

            tmpl = tmpd * raiLapl(a, i) + raiSquare(a, i) * &
                   (k * (k - 1) * enPowers(a, i, k - 2) * enPowers(a, j, l) + &
                    l * (l - 1) * enPowers(a, i, l - 2) * enPowers(a, j, k))
            fLapli(i) = fLapli(i) + gamma(g, c) * tmpl
            fkLapli(i, g) = tmpl

            tmpd = l * enPowers(a, i, k) * enPowers(a, j, l - 1) + &
                   k * enPowers(a, i, l) * enPowers(a, j, k - 1)
            tmpn(:) = raiDeriv(:, a, j) * tmpd
            fDeriv(3*j-2:3*j) = fDeriv(3*j-2:3*j) + gamma(g, c) * tmpn(:)
            fkDeriv(3*j-2:3*j, g) = tmpn(:)

            tmpl = tmpd * raiLapl(a, j) + raiSquare(a, j) * &
                   (k * (k - 1) * enPowers(a, i, l) * enPowers(a, j, k - 2) + &
                    l * (l - 1) * enPowers(a, i, k) * enPowers(a, j, l - 2))
            fLapli(j) = fLapli(j) + gamma(g, c) * tmpl
            fkLapli(j, g) = tmpl
          enddo

          !params = params + p / 2 - 1

          ! free terms with linear r_ij
          do k = 2, (p - 1)/2
            g = g + 1
            l = p - k
            !g = params + k - 2
            tmp = enPowers(a, i, p - 1) + enPowers(a, j, p - 1) - &
                  enPowers(a, i, k)     * enPowers(a, j, l - 1) - &
                  enPowers(a, i, l - 1) * enPowers(a, j, k)
            fTerm = fTerm + gamma(g, c) * eePowers(i, j, 1) * tmp

            if( calcDerivs == 0) cycle

            fkTerm(g) = eePowers(i, j, 1) * tmp

            tmpd = (p - 1) * enPowers(a, i, p - 2) - &
                   k       * enPowers(a, i, k - 1) * enPowers(a, j, l - 1) - &
                   (l - 1) * enPowers(a, i, l - 2) * enPowers(a, j, k)
            tmpn(:) = rijDeriv(:, i, j) * tmp + &
                   eePowers(i, j, 1) * raiDeriv(:, a, i) * tmpd
            fDeriv(3*i-2:3*i) = fDeriv(3*i-2:3*i) + gamma(g, c) * tmpn(:)
            fkDeriv(3*i-2:3*i, g) = tmpn(:)

            tmpl = rijLapl(i, j) * tmp + &
                   2 * rComb(1) * tmpd + &
                   eePowers(i, j, 1) * &
                   (raiLapl(a, i) * tmpd + raiSquare(a, i) * &
                    ((p - 1) * (p - 2) * enPowers(a, i, p - 3) - &
                     k       * (k - 1) * enPowers(a, i, k - 2) * enPowers(a, j, l - 1) - &
                     (l - 1) * (l - 2) * enPowers(a, i, l - 3) * enPowers(a, j, k)))
            fLapli(i) = fLapli(i) + gamma(g, c) * tmpl
            fkLapli(i, g) = tmpl

            tmpd = (p - 1) * enPowers(a, j, p - 2) - &
                   (l - 1) * enPowers(a, i, k) * enPowers(a, j, l - 2) - &
                   k       * enPowers(a, i, l - 1) * enPowers(a, j, k - 1)
            tmpn(:) = -tmp * rijDeriv(:, i, j) + &
                    eePowers(i, j, 1) * raiDeriv(:, a, j) * tmpd
            fDeriv(3*j-2:3*j) = fDeriv(3*j-2:3*j) + gamma(g, c) * tmpn(:)
            fkDeriv(3*j-2:3*j, g) = tmpn(:)

            tmpl = rijLapl(i, j) * tmp + &
                   2 * rComb(2) * tmpd + &
                   eePowers(i, j, 1) * &
                   (raiLapl(a, j) * tmpd + raiSquare(a, j) * &
                    ((p - 1) * (p - 2)                         * enPowers(a, j, p - 3) - &
                     k       * (k - 1) * enPowers(a, i, l - 1) * enPowers(a, j, k - 2) - &
                     (l - 1) * (l - 2) * enPowers(a, i, k)     * enPowers(a, j, l - 3)))
            fLapli(j) = fLapli(j) + gamma(g, c) * tmpl
            fkLapli(j, g) = tmpl
          enddo

          !params = params + (p-1) / 2 - 1

          ! free terms with r_ij^k, k > 1
          do k = 2, p - 2
            g = g + 1
            l = p - k
            !   g = params + k - 2

            !   if (p >= 6 .and. k > 2) then
            !     do x = 2, k - 1
            !       g = g + (p - x) / 2 - 1;
            !     enddo
            !   endif

            tmp = enPowers(a, i, k) + enPowers(a, j, k)
            fTerm = fTerm + gamma(g, c) * eePowers(i, j, l) * tmp

            if(calcDerivs /= 0) then
              fkTerm(g) = eePowers(i, j, l) * tmp

              tmpd = k * enPowers(a, i, k - 1)
              tmpn(:) = l * eePowers(i, j, l - 1) * rijDeriv(:, i, j) * tmp + &
                         eePowers(i, j, l)     * raiDeriv(:, a, i) * tmpd
              fDeriv(3*i-2:3*i) = fDeriv(3*i-2:3*i) + gamma(g, c) * tmpn(:)
              fkDeriv(3*i-2:3*i, g) = tmpn(:)

              tmpl = l * (l - 1) * eePowers(i, j, l - 2) * rijSquare(i, j) * tmp  + &
                     l           * eePowers(i, j, l - 1) * rijLapl(i, j)   * tmp  + &
                     2 * l       * eePowers(i, j, l - 1) * rComb(1)        * tmpd + &
                                   eePowers(i, j, l)     * raiLapl(a, i)   * tmpd + &
                                   eePowers(i, j, l)     * raiSquare(a, i) * k * (k - 1) * enPowers(a, i, k - 2)
              fLapli(i) = fLapli(i) + gamma(g, c) * tmpl
              fkLapli(i, g) = tmpl

              tmpd = k * enPowers(a, j, k - 1)
              tmpn(:) = -l * eePowers(i, j, l - 1) * rijDeriv(:, i, j) * tmp + &
                          eePowers(i, j, l)     * raiDeriv(:, a, j) * tmpd
              fDeriv(3*j-2:3*j) = fDeriv(3*j-2:3*j) + gamma(g, c) * tmpn(:)
              fkDeriv(3*j-2:3*j, g) = tmpn(:)

              tmpl = l * (l - 1) * eePowers(i, j, l - 2) * rijSquare(i, j) * tmp  + &
                     l           * eePowers(i, j, l - 1) * rijLapl(i, j)   * tmp  + &
                     2 * l       * eePowers(i, j, l - 1) * rComb(2)        * tmpd + &
                                   eePowers(i, j, l)     * raiLapl(a, j)   * tmpd + &
                                   eePowers(i, j, l)     * raiSquare(a, j) * k * (k - 1) * enPowers(a, j, k - 2)
              fLapli(j) = fLapli(j) + gamma(g, c) * tmpl
              fkLapli(j, g) = tmpl
            endif

            do m = 2, l/2
              !gt = g + m - 1
              g = g + 1
              ! duplicates for m = k/2, see above
              tmp = enPowers(a, i, m) * enPowers(a, j, l - m) + &
                    enPowers(a, i, l - m) * enPowers(a, j, m)
              fTerm = fTerm + gamma(g, c) * eePowers(i, j, k) * tmp

              if(calcDerivs == 0) cycle

              fkTerm(g) = eePowers(i, j, k) * tmp

              tmpd = (l-m) * enPowers(a, i, l - m - 1) * enPowers(a, j, m) + &
                     m     * enPowers(a, i, m - 1)     * enPowers(a, j, l - m)
              tmpn(:) = k * eePowers(i, j, k - 1) * rijDeriv(:, i, j) * tmp + &
                         eePowers(i, j, k)     * raiDeriv(:, a, i) * tmpd
              fDeriv(3*i-2:3*i) = fDeriv(3*i-2:3*i) + gamma(g, c) * tmpn(:)
              fkDeriv(3*i-2:3*i, g) = tmpn(:)

              tmpl = k * (k - 1) * eePowers(i, j, k - 2) * rijSquare(i, j) * tmp  + &
                     k           * eePowers(i, j, k - 1) * rijLapl(i, j)   * tmp  + &
                     2 * k       * eePowers(i, j, k - 1) * rComb(1)        * tmpd + &
                                   eePowers(i, j, k)     * raiLapl(a, i)   * tmpd + &
                                   eePowers(i, j, k)     * raiSquare(a, i) * &
                     (m * (m - 1)           * enPowers(a, i, m - 2)     * enPowers(a, j, l - m) + &
                      (l - m) * (l - m - 1) * enPowers(a, i, l - m - 2) * enPowers(a, j, m))
              fLapli(i) = fLapli(i) + gamma(g, c) * tmpl
              fkLapli(i, g) = tmpl


              tmpd = (l-m) * enPowers(a, i, m) * enPowers(a, j, l - m - 1) + &
                     m     * enPowers(a, i, l - m) * enPowers(a, j, m - 1)
              tmpn(:) = -k * eePowers(i, j, k - 1) * rijDeriv(:, i, j) * tmp + &
                          eePowers(i, j, k)     * raiDeriv(:, a, j) * tmpd
              fDeriv(3*j-2:3*j) = fDeriv(3*j-2:3*j) + gamma(g, c) * tmpn(:)
              fkDeriv(3*j-2:3*j, g) = tmpn(:)

              tmpl = k * (k - 1) * eePowers(i, j, k - 2) * rijSquare(i, j) * tmp  + &
                     k           * eePowers(i, j, k - 1) * rijLapl(i, j)   * tmp  + &
                     2 * k       * eePowers(i, j, k - 1) * rComb(2)        * tmpd + &
                                   eePowers(i, j, k)     * raiLapl(a, j)   * tmpd + &
                                   eePowers(i, j, k)     * raiSquare(a, j) * &
                     (m * (m - 1)           * enPowers(a, i, l - m) * enPowers(a, j, m - 2) + &
                      (l - m) * (l - m - 1) * enPowers(a, i, m)     * enPowers(a, j, l - m - 2))
              fLapli(j) = fLapli(j) + gamma(g, c) * tmpl
              fkLapli(j, g) = tmpl
            enddo
          enddo

          !params = params + ((p - 2)**2) / 4
          l = p - 1
          ! in the following, non-free means terms with linear r_i
          ! non-free because they don't satisfy the cusp by themselves,
          ! but linear combinations do so
          ! non-free term without r_ij
          tmp = enPowers(a, i, 1) * enPowers(a, j, l) + &
                enPowers(a, i, l) * enPowers(a, j, 1)
          nfTerms(1) = tmp

          if (calcDerivs /= 0) then
            tmpd = enPowers(a, j, l) + &
                   l * enPowers(a, i, l - 1) * enPowers(a, j, 1)
            nfDerivs(1, 1:3) = tmpd * raiDeriv(:, a, i)

            nfLapl(1, 1) = raiLapl(a, i)   * tmpd + &
                           raiSquare(a, i) * l * (l - 1) * enPowers(a, i, l - 2) * enPowers(a, j, 1)

            tmpd = enPowers(a, i, l) + &
                   l * enPowers(a, i, 1) * enPowers(a, j, l - 1)
            nfDerivs(1, 4:6) = tmpd * raiDeriv(:, a, j)

            nfLapl(1, 2) = raiLapl(a, j)   * tmpd + &
                           raiSquare(a, j) * l * (l - 1) * enPowers(a, i, 1) * enPowers(a, j, l - 2)
          endif

          ! non-free term with linear r_ij
          if(l == 2) then
            ! whole term is duplicated, so divide by 2 to satisfy cusp exactly
            coeff = -0.5d0
          else
            coeff = -1d0
          endif
          tmp = enPowers(a, i, l) + enPowers(a, j, l) - &
                enPowers(a, i, 1) * enPowers(a, j, l - 1) - &
                enPowers(a, i, l - 1) * enPowers(a, j, 1)
          nfTerms(2) = coeff * eePowers(i, j, 1) * tmp

          if (calcDerivs /= 0) then
            tmpd = l       * enPowers(a, i, l - 1) - &
                             enPowers(a, j, l - 1) - &
                   (l - 1) * enPowers(a, i, l - 2) * enPowers(a, j, 1)
            nfDerivs(2, 1:3) = coeff * (rijDeriv(:, i, j) * tmp + &
                                        eePowers(i, j, 1) * raiDeriv(:, a, i) * tmpd)

            nfLapl(2, 1) = coeff * (&
                            rijLapl(i, j) * tmp + &
                            2 * rComb(1) * tmpd + &
                            eePowers(i, j, 1) * &
                            (raiLapl(a, i) * tmpd + raiSquare(a, i) * &
                             (l * (l - 1)       * enPowers(a, i, l - 2) - &
                              (l - 1) * (l - 2) * enPowers(a, i, l - 3) * enPowers(a, j, 1))))

            tmpd = l       * enPowers(a, j, l - 1) - &
                             enPowers(a, i, l - 1) - &
                   (l - 1) * enPowers(a, j, l - 2) * enPowers(a, i, 1)
            nfDerivs(2, 4:6) = coeff * (-rijDeriv(:, i, j) * tmp + &
                                        eePowers(i, j, 1) * raiDeriv(:, a, j) * tmpd)

            nfLapl(2, 2) = coeff * ( &
                            rijLapl(i, j) * tmp + &
                            2 * rComb(2) * tmpd + &
                            eePowers(i, j, 1) * &
                            (raiLapl(a, j) * tmpd + raiSquare(a, j) * &
                             (l * (l - 1)                           * enPowers(a, j, l - 2) - &
                              (l - 1) * (l - 2) * enPowers(a, i, 1) * enPowers(a, j, l - 3))))
          endif

          ! non-free term with r_ij^k, k > 1
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfTerms(3) = eePowers(i, j, l) * tmp

          if (calcDerivs /= 0) then
            ! tmpd = 1
            nfDerivs(3, 1:3) = l * eePowers(i, j, l - 1) * rijDeriv(:, i, j) * tmp + &
                                   eePowers(i, j, l)     * raiDeriv(:, a, i)

            nfLapl(3, 1) = l * (l - 1) * eePowers(i, j, l - 2) * rijSquare(i, j) * tmp + &
                           l           * eePowers(i, j, l - 1) * rijLapl(i, j)   * tmp + &
                           2 * l       * eePowers(i, j, l - 1) * rComb(1)              + &
                                         eePowers(i, j, l)     * raiLapl(a, i)

            nfDerivs(3, 4:6) = -l * eePowers(i, j, l - 1) * rijDeriv(:, i, j) * tmp + &
                                    eePowers(i, j, l)     * raiDeriv(:, a, j)

            nfLapl(3, 2) = l * (l - 1) * eePowers(i, j, l - 2) * rijSquare(i, j) * tmp + &
                           l           * eePowers(i, j, l - 1) * rijLapl(i, j)   * tmp + &
                           2 * l       * eePowers(i, j, l - 1) * rComb(2)              + &
                                         eePowers(i, j, l)     * raiLapl(a, j)
          endif

          ! non-free terms with products of r_i and r_j
          do m = 2, p - 2
            k = l - m

            if(k == 1) then
              coeff = 0.5d0
            else
              coeff = 1d0
            endif

            tmp = enPowers(a, i, 1) * enPowers(a, j, k) + &
                  enPowers(a, i, k) * enPowers(a, j, 1)
            nfTerms(2 + m) = coeff * eePowers(i, j, m) * tmp

            if (calcDerivs == 0) cycle

            tmpd = k * enPowers(a, i, k - 1) * enPowers(a, j, 1) + &
                                               enPowers(a, j, k)
            nfDerivs(2 + m, 1:3) = coeff * (m * eePowers(i, j, m - 1) * rijDeriv(:, i, j) * tmp + &
                                                eePowers(i, j, m)     * raiDeriv(:, a, i) * tmpd)

            nfLapl(2 + m, 1) = coeff * (m * (m - 1) * eePowers(i, j, m - 2) * rijSquare(i, j) * tmp  + &
                                        m           * eePowers(i, j, m - 1) * rijLapl(i, j)   * tmp  + &
                                        2 * m       * eePowers(i, j, m - 1) * rComb(1)        * tmpd + &
                                                      eePowers(i, j, m)     * raiLapl(a, i)   * tmpd + &
                                                      eePowers(i, j, m)     * raiSquare(a, i) * &
                                        (k * (k - 1) * enPowers(a, i, k - 2) * enPowers(a, j, 1)))

            tmpd = (p - m - 1) * enPowers(a, i, 1) * enPowers(a, j, k - 1) + &
                                 enPowers(a, i, k)
            nfDerivs(2 + m, 4:6) = coeff * (-m * eePowers(i, j, m - 1) * rijDeriv(:, i, j) * tmp + &
                                                 eePowers(i, j, m)     * raiDeriv(:, a, j) * tmpd)

            nfLapl(2 + m, 2) = coeff * (m * (m - 1) * eePowers(i, j, m - 2) * rijSquare(i, j) * tmp  + &
                                        m           * eePowers(i, j, m - 1) * rijLapl(i, j)   * tmp  + &
                                        2 * m       * eePowers(i, j, m - 1) * rComb(2)        * tmpd + &
                                                      eePowers(i, j, m)     * raiLapl(a, j)   * tmpd + &
                                                      eePowers(i, j, m)     * raiSquare(a, j) * &
                                        (k * (k - 1) * enPowers(a, i, 1) * enPowers(a, j, k - 2)))
          enddo

          !gend = params
          gend = g

          ! generate the linear combinations of non-free terms
          ! to do so, just substract one of the terms from all others
          do k = 1, p - 1
            !g = params + k - 1
            g = g + 1
            fTerm = fTerm + gamma(g, c) * (nfTerms(k) - nfTerms(p))

            if(calcDerivs == 0) cycle

            fDeriv(3*i-2:3*i) = fDeriv(3*i-2:3*i) + gamma(g, c) * (nfDerivs(k, 1:3) - nfDerivs(p, 1:3))
            fDeriv(3*j-2:3*j) = fDeriv(3*j-2:3*j) + gamma(g, c) * (nfDerivs(k, 4:6) - nfDerivs(p, 4:6))
            fLapli(i) = fLapli(i) + gamma(g, c) * (nfLapl(k, 1) - nfLapl(p, 1))
            fLapli(j) = fLapli(j) + gamma(g, c) * (nfLapl(k, 2) - nfLapl(p, 2))
          enddo

          !params = params + p - 1

          if(deriveLinParams) then
            ! g has already been incremented, so only loop until gend - 1
            do t = gstart, gend
              m = unum + xnum + (c - 1) * fpnum + t
              uk(m) = uk(m) + fkTerm(t)
              ukgrad(3*i-2:3*i, m) = ukgrad(3*i-2:3*i, m) + fkDeriv(3*i-2:3*i, t)
              ukgrad(3*j-2:3*j, m) = ukgrad(3*j-2:3*j, m) + fkDeriv(3*j-2:3*j, t)
              uklapli(i, m) = uklapli(i, m) + fkLapli(i, t)
              uklapli(j, m) = uklapli(j, m) + fkLapli(j, t)
              uklapl(m) = uklapl(m) + fkLapli(i, t) + fkLapli(j, t)
            enddo

            do k = 1, p - 1
              t = gend + k
              m = unum + xnum + (c - 1) * fpnum + t
              uk(m) = uk(m) + nfTerms(k) - nfTerms(p)
              ukgrad(3*i-2:3*i, m) = ukgrad(3*i-2:3*i, m) + nfDerivs(k, 1:3) - nfDerivs(p, 1:3)
              ukgrad(3*j-2:3*j, m) = ukgrad(3*j-2:3*j, m) + nfDerivs(k, 4:6) - nfDerivs(p, 4:6)
              uklapli(i, m) = uklapli(i, m) + nfLapl(k, 1) - nfLapl(p, 1)
              uklapli(j, m) = uklapli(j, m) + nfLapl(k, 2) - nfLapl(p, 2)
              uklapl(m) = uklapl(m) + nfLapl(k, 1) - nfLapl(p, 1) + &
                                      nfLapl(k, 2) - nfLapl(p, 2)
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
end subroutine eenGeneric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  _____ ______  __     ______  _    _
!!!! |_   _|  ____| \ \   / / __ \| |  | |
!!!!   | | | |__     \ \_/ / |  | | |  | |
!!!!   | | |  __|     \   /| |  | | |  | |
!!!!  _| |_| |         | | | |__| | |__| |
!!!! |_____|_|    _    |_|  \____/ \____/ ______
!!!!  / ____| |  | |   /\   | \ | |/ ____|  ____|
!!!! | |    | |__| |  /  \  |  \| | |  __| |__
!!!! | |    |  __  | / /\ \ | . ` | | |_ |  __|
!!!! | |____| |  | |/ ____ \| |\  | |__| | |____
!!!!  \_____|_|_ |_/_/    \_\_|_\_|\_____|______| _   _  _____
!!!!     /\   | \ | \ \   / /__   __| |  | |_   _| \ | |/ ____|
!!!!    /  \  |  \| |\ \_/ /   | |  | |__| | | | |  \| | |  __
!!!!   / /\ \ | . ` | \   /    | |  |  __  | | | | . ` | | |_ |
!!!!  / ____ \| |\  |  | |     | |  | |  | |_| |_| |\  | |__| |
!!!! /_/   _\_\_|_\_|__|_| ____|_|  |_|  |_|_____|_| \_|\_____|
!!!! | |  | |  ____|  __ \|  ____|
!!!! | |__| | |__  | |__) | |__
!!!! |  __  |  __| |  _  /|  __|
!!!! | |  | | |____| | \ \| |____ _
!!!! |_|  |_|______|_|  \_\______( )
!!!!  _____  ______ __  __ ______|/_  __ ____  ______ _____
!!!! |  __ \|  ____|  \/  |  ____|  \/  |  _ \|  ____|  __ \
!!!! | |__) | |__  | \  / | |__  | \  / | |_) | |__  | |__) |
!!!! |  _  /|  __| | |\/| |  __| | |\/| |  _ <|  __| |  _  /
!!!! | | \ \| |____| |  | | |____| |  | | |_) | |____| | \ \
!!!! |_|__\_\______|_|  |_|______|_|  |_|____/|______|_|  \_\
!!!! |__   __/ __ \
!!!!    | | | |  | |
!!!!    | | | |  | |
!!!!    | | | |__| |
!!!!   _|_|_ \____/          _   _  _____ ______
!!!!  / ____| |  | |   /\   | \ | |/ ____|  ____|
!!!! | |    | |__| |  /  \  |  \| | |  __| |__
!!!! | |    |  __  | / /\ \ | . ` | | |_ |  __|
!!!! | |____| |  | |/ ____ \| |\  | |__| | |____
!!!!  \_____|_|__|_/_/  _ \_\_|_\_|\_____|______|_____   _____ ____  _____  ______ _____
!!!! |  ____|  ____| \ | | |  | |   /\   |  __ \|  __ \ / ____/ __ \|  __ \|  ____|  __ \
!!!! | |__  | |__  |  \| | |__| |  /  \  | |__) | |  | | |   | |  | | |  | | |__  | |  | |______
!!!! |  __| |  __| | . ` |  __  | / /\ \ |  _  /| |  | | |   | |  | | |  | |  __| | |  | |______|
!!!! | |____| |____| |\  | |  | |/ ____ \| | \ \| |__| | |___| |__| | |__| | |____| |__| |
!!!! |______|______|_|_\_|_|__|_/_/___ \_\_|__\_\_____/_\_____\____/|_____/|______|_____/
!!!! | \ | |/ __ \|  __ \|  ____|  __ \|_   _\ \    / / ____|
!!!! |  \| | |  | | |  | | |__  | |__) | | |  \ \  / / (___
!!!! | . ` | |  | | |  | |  __| |  _  /  | |   \ \/ / \___ \
!!!! | |\  | |__| | |__| | |____| | \ \ _| |_   \  /  ____) |
!!!! |_| \_|\____/|_____/|______|_|  \_\_____|   \/  |_____/
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eenHardcodedAll(fTerm, fDeriv, fLapli, deriveLinParams)
  ! value of f term + derivs
  real(r8), intent(inout) :: fTerm, fDeriv(:), fLapli(:)
  logical, intent(in) :: deriveLinParams

  real(r8) :: rComb(2)

  real(r8) :: tmp, tmpd, tmpn(6), tmpl
  real(r8) :: terms(fpnum), termsn(6, fpnum), termsl(2, fpnum)
  real(r8) :: nfterms(fmax), nftermsd(2, fmax), nftermsn(6, fmax), nftermsl(2, fmax)
  integer :: a, c, g, i, j, m, t, ii, jj

  do a = 1, nclast
    c = atoms(a)%sa

    do i = 1, ne
      do j = i + 1, ne
        !DIR$ vector
        rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
        !DIR$ vector
        rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

        if(fmax > 2) then
          ! r_i^2 r_j + r_i r_j^2
          nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                       enPowers(a, i, 1) * enPowers(a, j, 2)

          nftermsd(1, 1) = 2 * enPowers(a, i, 1) * enPowers(a, j, 1) + &
                                                   enPowers(a, j, 2)
          nftermsd(2, 1) =     enPowers(a, i, 2) + &
                           2 * enPowers(a, i, 1) * enPowers(a, j, 1)
          nftermsn(1:3, 1) = nftermsd(1, 1) * raiDeriv(1:3, a, i)
          nftermsn(4:6, 1) = nftermsd(2, 1) * raiDeriv(1:3, a, j)
          nftermsl(1, 1) = 2 * enPowers(a, j, 1) * raiSquare(a, i) + nftermsd(1, 1) * raiLapl(a, i)
          nftermsl(2, 1) = 2 * enPowers(a, i, 1) * raiSquare(a, j) + nftermsd(2, 1) * raiLapl(a, j)

          ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                2 * enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(2) = -0.5d0 * eePowers(i, j, 1) * tmp

          nftermsd(1, 2) = 2 * enPowers(a, i, 1) - 2 * enPowers(a, j, 1)
          nftermsd(2, 2) = -nftermsd(1, 2)
          tmpn(1:3) = rijDeriv(1:3, i, j) * tmp
          nftermsn(1:3, 2) = -0.5d0 * (tmpn(1:3) + eePowers(i, j, 1) * raiDeriv(1:3, a, i) * nftermsd(1, 2))
          nftermsn(4:6, 2) =  0.5d0 * (tmpn(1:3) - eePowers(i, j, 1) * raiDeriv(1:3, a, j) * nftermsd(2, 2))
          tmpl = rijLapl(i, j) * tmp
          nftermsl(1, 2) = -0.5d0 * (tmpl + 2 * rComb(1) * nftermsd(1, 2) + &
                                     eePowers(i, j, 1) * (raiLapl(a, i) * nftermsd(1, 2) + 2 * raiSquare(a, i)))
          nftermsl(2, 2) = -0.5d0 * (tmpl + 2 * rComb(2) * nftermsd(2, 2) + &
                                     eePowers(i, j, 1) * (raiLapl(a, j) * nftermsd(2, 2) + 2 * raiSquare(a, j)))

          ! r_ij^2 (r_i + r_j)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 2) * tmp

          tmpn(1:3) = 2 * eePowers(i, j, 1) * rijDeriv(1:3, i, j) * tmp
          nftermsn(1:3, 3) =  tmpn(1:3) + eePowers(i, j, 2) * raiDeriv(1:3, a, i)
          nftermsn(4:6, 3) = -tmpn(1:3) + eePowers(i, j, 2) * raiDeriv(1:3, a, j)
          tmpd = 2 * rijSquare(i, j) * tmp + &
                 2 * eePowers(i, j, 1) * rijLapl(i, j) * tmp
          nftermsl(1, 3) = tmpd + &
                           4 * eePowers(i, j, 1) * rComb(1) + &
                               eePowers(i, j, 2) * raiLapl(a, i)
          nftermsl(2, 3) = tmpd + &
                           4 * eePowers(i, j, 1) * rComb(2) + &
                               eePowers(i, j, 2) * raiLapl(a, j)

          do t = 1, 2
            terms(t) = nfterms(t) - nfterms(3)

            termsn(:, t) = nftermsn(:, t) - nftermsn(:, 3)
            termsl(:, t) = nftermsl(:, t) - nftermsl(:, 3)
          enddo
        endif

        if(fmax > 3) then
          ! 2 ri^2 rj^2
          terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)

          tmpd = 4 * enPowers(a, i, 1) * enPowers(a, j, 2)
          termsn(1:3, 3) = tmpd * raiDeriv(:, a, i)
          termsl(1, 3) = 4 * enPowers(a, j, 2) * raiSquare(a, i) + tmpd * raiLapl(a, i)

          tmpd = 4 * enPowers(a, i, 2) * enPowers(a, j, 1)
          termsn(4:6, 3) = tmpd * raiDeriv(:, a, j)
          termsl(2, 3) = 4 * enPowers(a, i, 2) * raiSquare(a, j) + tmpd * raiLapl(a, j)

          ! rij^2 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(4) = eePowers(i, j, 2) * tmp

          tmpn(1:3) = eePowers(i, j, 1) * rijDeriv(:, i, j) * tmp
          termsn(1:3, 4) = 2.0D0 * ( tmpn(1:3) + eePowers(i, j, 2) * enPowers(a, i, 1) * raiDeriv(:, a, i))
          termsn(4:6, 4) = 2.0D0 * (-tmpn(1:3) + eePowers(i, j, 2) * enPowers(a, j, 1) * raiDeriv(:, a, j))
          tmpd = (eePowers(i, j, 1) * rijLapl(i, j) + rijSquare(i, j)) * tmp
          termsl(1, 4) = 2 * (4 * eePowers(i, j, 1) * rComb(1) * enPowers(a, i, 1) + &
                                  eePowers(i, j, 2) * raiLapl(a, i) * enPowers(a, i, 1) + &
                                  eePowers(i, j, 2) * raiSquare(a, i) + &
                                  tmpd)
          termsl(2, 4) = 2 * (4 * eePowers(i, j, 1) * rComb(2) * enPowers(a, j, 1) + &
                                  eePowers(i, j, 2) * raiLapl(a, j) * enPowers(a, j, 1) + &
                                  eePowers(i, j, 2) * raiSquare(a, j) + &
                                  tmpd)

          ! ri^3 rj + ri rj^3
          nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 3)

          nftermsd(1, 1) = 3 * enPowers(a, i, 2) * enPowers(a, j, 1) + enPowers(a, j, 3)
          nftermsd(2, 1) = 3 * enPowers(a, j, 2) * enPowers(a, i, 1) + enPowers(a, i, 3)
          nftermsn(1:3, 1) = nftermsd(1, 1) * raiDeriv(:, a, i)
          nftermsn(4:6, 1) = nftermsd(2, 1) * raiDeriv(:, a, j)

          tmp = 6 * enPowers(a, i, 1) * enPowers(a, j, 1)
          nftermsl(1, 1) = nftermsd(1, 1) * raiLapl(a, i) + tmp * raiSquare(a, i)
          nftermsl(2, 1) = nftermsd(2, 1) * raiLapl(a, j) + tmp * raiSquare(a, j)

          ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                enPowers(a, i, 2) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(2) = - eePowers(i, j, 1) * tmp

          nftermsd(1, 2) = 3 * enPowers(a, i, 2) - &
                           2 * enPowers(a, i, 1) * enPowers(a, j, 1) - &
                                                   enPowers(a, j, 2)
          nftermsd(2, 2) = 3 * enPowers(a, j, 2) - &
                           2 * enPowers(a, j, 1) * enPowers(a, i, 1) - &
                                                   enPowers(a, i, 2)
          tmpn(1:3) = rijDeriv(1:3, i, j) * tmp
          nftermsn(1:3, 2) = -tmpn(1:3) - eePowers(i, j, 1) * nftermsd(1, 2) * raiDeriv(:, a, i)
          nftermsn(4:6, 2) =  tmpn(1:3) - eePowers(i, j, 1) * nftermsd(2, 2) * raiDeriv(:, a, j)

          tmpd = -rijLapl(i, j) * tmp
          nftermsl(1, 2) = tmpd - 2 * rComb(1) * nftermsd(1, 2) - &
                           eePowers(i, j, 1) * (raiLapl(a, i) * nftermsd(1, 2) + &
                            raiSquare(a, i) * (6 * enPowers(a, i, 1) - 2 * enPowers(a, j, 1)))
          nftermsl(2, 2) = tmpd - 2 * rComb(2) * nftermsd(2, 2) - &
                           eePowers(i, j, 1) * (raiLapl(a, j) * nftermsd(2, 2) + &
                            raiSquare(a, j) * (6 * enPowers(a, j, 1) - 2 * enPowers(a, i, 1)))

          ! rij^3 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 3) * tmp

          tmpn(1:3) = 3 * eePowers(i, j, 2) * rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 3) =  tmpn(1:3) + eePowers(i, j, 3) * raiDeriv(:, a, i)
          nftermsn(4:6, 3) = -tmpn(1:3) + eePowers(i, j, 3) * raiDeriv(:, a, j)
          tmpd = 3 * tmp * (2 * eePowers(i, j, 1) * rijSquare(i, j) + &
                                eePowers(i, j, 2) * rijLapl(i, j))
          nftermsl(1, 3) = tmpd + 6 * eePowers(i, j, 2) * rComb(1) + eePowers(i, j, 3) * raiLapl(a, i)

          nftermsl(2, 3) = tmpd + 6 * eePowers(i, j, 2) * rComb(2) + eePowers(i, j, 3) * raiLapl(a, j)

          ! rij^2 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(4) = eePowers(i, j, 2) * tmp

          tmpn(1:3) = 2 * eePowers(i, j, 1) * rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 4) =  tmpn(1:3) + eePowers(i, j, 2) * enPowers(a, j, 1) * raiDeriv(:, a, i)
          nftermsn(4:6, 4) = -tmpn(1:3) + eePowers(i, j, 2) * enPowers(a, i, 1) * raiDeriv(:, a, j)
          tmpd = 2 * (tmp * rijSquare(i, j) + eePowers(i, j, 1) * tmp * rijLapl(i, j))
          nftermsl(1, 4) = tmpd + eePowers(i, j, 2) * enPowers(a, j, 1) * raiLapl(a, i) + &
                              4 * eePowers(i, j, 1) * enPowers(a, j, 1) * rComb(1)
          nftermsl(2, 4) = tmpd + eePowers(i, j, 2) * enPowers(a, i, 1) * raiLapl(a, j) + &
                              4 * eePowers(i, j, 1) * enPowers(a, i, 1) * rComb(2)

          do t = 1, 3
            terms(4+t) = nfterms(t) - nfterms(4)

            termsn(:, 4+t) = nftermsn(:, t) - nftermsn(:, 4)
            termsl(:, 4+t) = nftermsl(:, t) - nftermsl(:, 4)
          enddo
        endif

        if(fmax > 4) then
          ! ri^3 rj^2 + ri^2 rj^3
          terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                     enPowers(a, i, 2) * enPowers(a, j, 3)

          termsn(1:3, 8) = (3 * enPowers(a, i, 2) * enPowers(a, j, 2) + &
                            2 * enPowers(a, i, 1) * enPowers(a, j, 3)) * raiDeriv(:, a, i)
          termsn(4:6, 8) = (3 * enPowers(a, j, 2) * enPowers(a, i, 2) + &
                            2 * enPowers(a, j, 1) * enPowers(a, i, 3)) * raiDeriv(:, a, j)
          termsl(1, 8) = 3 * enPowers(a, j, 2) * (2 * enPowers(a, i, 1) * raiSquare(a, i) + &
                                                      enPowers(a, i, 2) * raiLapl(a, i)) + &
                         2 * enPowers(a, j, 3) * (raiSquare(a, i) + enPowers(a, i, 1) * raiLapl(a, i))
          termsl(2, 8) = 3 * enPowers(a, i, 2) * (2 * enPowers(a, j, 1) * raiSquare(a, j) + &
                                                      enPowers(a, j, 2) * raiLapl(a, j)) + &
                         2 * enPowers(a, i, 3) * (raiSquare(a, j) + enPowers(a, j, 1) * raiLapl(a, j))

          ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
          tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                2 * enPowers(a, i, 2) * enPowers(a, j, 2)
          terms(9)  = eePowers(i, j, 1) * tmp

          tmpd = 4 * (enPowers(a, i, 3) - enPowers(a, i, 1) * enPowers(a, j, 2))
          termsn(1:3, 9) = rijDeriv(:, i, j) * tmp + eePowers(i, j, 1) * tmpd * raiDeriv(1:3, a, i)
          tmpl = rijLapl(i, j) * tmp
          termsl(1, 9) = tmpl + 2 * rComb(1) * tmpd + eePowers(i, j, 1) * &
                          (tmpd * raiLapl(a, i) + raiSquare(a, i) * (12 * enPowers(a, i, 2) - 4 * enPowers(a, j, 2)))

          tmpd = 4 * (enPowers(a, j, 3) - enPowers(a, j, 1) * enPowers(a, i, 2))
          termsn(4:6, 9) = -rijDeriv(:, i, j) * tmp + eePowers(i, j, 1) * tmpd * raiDeriv(1:3, a, j)
          termsl(2, 9) = tmpl + 2 * rComb(2) * tmpd + eePowers(i, j, 1) * &
                          (tmpd * raiLapl(a, j) + raiSquare(a, j) * (12 * enPowers(a, j, 2) - 4 * enPowers(a, i, 2)))

          ! rij^3 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(10) = eePowers(i, j, 3) * tmp

          tmpd = 2 * enPowers(a, i, 1)
          tmpn(1:3) = 3 * eePowers(i, j, 2) * rijDeriv(:, i, j) * tmp
          termsn(1:3, 10) = tmpn(1:3) + eePowers(i, j, 3) * raiDeriv(:, a, i) * tmpd
          tmpl = 6 * eePowers(i, j, 1) * rijSquare(i, j) * tmp + &
                 3 * eePowers(i, j, 2) * rijLapl(i, j) * tmp
          termsl(1, 10) = tmpl + 6 * eePowers(i, j, 2) * rComb(1) * tmpd + &
                          eePowers(i, j, 3) * (raiLapl(a, i) * tmpd + 2 * raiSquare(a, i))

          tmpd = 2 * enPowers(a, j, 1)
          termsn(4:6, 10) = -tmpn(1:3) + eePowers(i, j, 3) * raiDeriv(:, a, j) * tmpd
          termsl(2, 10) = tmpl + 6 * eePowers(i, j, 2) * rComb(2) * tmpd + &
                          eePowers(i, j, 3) * (raiLapl(a, j) * tmpd + 2 * raiSquare(a, j))

          ! rij^2 (ri^3 + rj^3)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
          terms(11) = eePowers(i, j, 2) * tmp

          tmpd = 3 * enPowers(a, i, 2)
          tmpn(1:3) = 2 * eePowers(i, j, 1) * rijDeriv(1:3, i, j) * tmp
          termsn(1:3, 11) = tmpn(1:3) + eePowers(i, j, 2) * raiDeriv(:, a, i) * tmpd
          tmpl = 2 * (rijSquare(i, j) + eePowers(i, j, 1) * rijLapl(i, j)) * tmp
          termsl(1, 11) = tmpl + 4 * eePowers(i, j, 1) * rComb(1) * tmpd + &
                          eePowers(i, j, 2) * (tmpd * raiLapl(a, i) + 6 * raiSquare(a, i) * enPowers(a, i, 1))

          tmpd = 3 * enPowers(a, j, 2)
          termsn(4:6, 11) = -tmpn(1:3) + eePowers(i, j, 2) * raiDeriv(:, a, j) * tmpd
          termsl(2, 11) = tmpl + 4 * eePowers(i, j, 1) * rComb(2) * tmpd + &
                          eePowers(i, j, 2) * (tmpd * raiLapl(a, j) + 6 * raiSquare(a, j) * enPowers(a, j, 1))

          ! ri^4 rj + ri rj^4
          nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 4)

          nftermsd(1, 1) = 4 * enPowers(a, i, 3) * enPowers(a, j, 1) + enPowers(a, j, 4)
          nftermsd(2, 1) = 4 * enPowers(a, j, 3) * enPowers(a, i, 1) + enPowers(a, i, 4)
          nftermsn(1:3, 1) = nftermsd(1, 1) * raiDeriv(:, a, i)
          nftermsn(4:6, 1) = nftermsd(2, 1) * raiDeriv(:, a, j)
          nftermsl(1, 1) = 12 * enPowers(a, i, 2) * enPowers(a, j, 1) * raiSquare(a, i) + nftermsd(1, 1) * raiLapl(a, i)
          nftermsl(2, 1) = 12 * enPowers(a, j, 2) * enPowers(a, i, 1) * raiSquare(a, j) + nftermsd(2, 1) * raiLapl(a, j)

          ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
          tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                enPowers(a, i, 3) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 3)
          nfterms(2) = -eePowers(i, j, 1) * tmp

          nftermsd(1, 2) = 4 * enPowers(a, i, 3) - &
                           3 * enPowers(a, i, 2) * enPowers(a, j, 1) - &
                                                   enPowers(a, j, 3)
          nftermsd(2, 2) = 4 * enPowers(a, j, 3) - &
                           3 * enPowers(a, j, 2) * enPowers(a, i, 1) - &
                                                   enPowers(a, i, 3)
          tmpn(1:3) = rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 2) = -tmpn(1:3) - eePowers(i, j, 1) * nftermsd(1, 2) * raiDeriv(:, a, i)
          nftermsn(4:6, 2) =  tmpn(1:3) - eePowers(i, j, 1) * nftermsd(2, 2) * raiDeriv(:, a, j)
          tmpd = -rijLapl(i, j) * tmp
          nftermsl(1, 2) = tmpd - 2 * rComb(1) * nftermsd(1, 2) - eePowers(i, j, 1) * &
                           (nftermsd(1, 2) * raiLapl(a, i) + 6 * raiSquare(a, i) * ( &
                            2 * enPowers(a, i, 2) - enPowers(a, i, 1) * enPowers(a, j, 1)))
          nftermsl(2, 2) = tmpd - 2 * rComb(2) * nftermsd(2, 2) - eePowers(i, j, 1) * &
                           (nftermsd(2, 2) * raiLapl(a, j) + 6 * raiSquare(a, j) * ( &
                            2 * enPowers(a, j, 2) - enPowers(a, j, 1) * enPowers(a, i, 1)))

          ! rij^4 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 4) * tmp

          tmpn(1:3) = 4 * eePowers(i, j, 3) * rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 3) =  tmpn(1:3) + eePowers(i, j, 4) * raiDeriv(:, a, i)
          nftermsn(4:6, 3) = -tmpn(1:3) + eePowers(i, j, 4) * raiDeriv(:, a, j)
          tmpd = (12 * eePowers(i, j, 2) * rijSquare(i, j) + &
                   4 * eePowers(i, j, 3) * rijLapl(i, j)) * tmp
          nftermsl(1, 3) = tmpd + 8 * eePowers(i, j, 3) * rComb(1) + &
                                      eePowers(i, j, 4) * raiLapl(a, i)
          nftermsl(2, 3) = tmpd + 8 * eePowers(i, j, 3) * rComb(2) + &
                                      eePowers(i, j, 4) * raiLapl(a, j)

          ! rij^2 (ri^2 rj + ri rj^2)
          tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(4) = eePowers(i, j, 2) * tmp

          nftermsd(1, 4) = 2 * enPowers(a, i, 1) * enPowers(a, j, 1) + enPowers(a, j, 2)
          nftermsd(2, 4) = 2 * enPowers(a, j, 1) * enPowers(a, i, 1) + enPowers(a, i, 2)
          tmpn(1:3) = 2 * eePowers(i, j, 1) * rijDeriv(1:3, i, j) * tmp
          nftermsn(1:3, 4) =  tmpn(1:3) + eePowers(i, j, 2) * nftermsd(1, 4) * raiDeriv(:, a, i)
          nftermsn(4:6, 4) = -tmpn(1:3) + eePowers(i, j, 2) * nftermsd(2, 4) * raiDeriv(:, a, j)
          tmpd = (2 * rijSquare(i, j) + 2 * eePowers(i, j, 1) * rijLapl(i, j)) * tmp
          nftermsl(1, 4) = tmpd + 4 * eePowers(i, j, 1) * rComb(1) * nftermsd(1, 4) + &
                           eePowers(i, j, 2) * (2 * enPowers(a, j, 1) * raiSquare(a, i) + nftermsd(1, 4) * raiLapl(a, i))
          nftermsl(2, 4) = tmpd + 4 * eePowers(i, j, 1) * rComb(2) * nftermsd(2, 4) + &
                           eePowers(i, j, 2) * (2 * enPowers(a, i, 1) * raiSquare(a, j) + nftermsd(2, 4) * raiLapl(a, j))

          ! rij^3 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(5) = eePowers(i, j, 3) * tmp

          tmpn(1:3) = 3 * eePowers(i, j, 2) * rijDeriv(1:3, i, j) * tmp
          nftermsn(1:3, 5) =  tmpn(1:3) + eePowers(i, j, 3) * enPowers(a, j, 1) * raiDeriv(:, a, i)
          nftermsn(4:6, 5) = -tmpn(1:3) + eePowers(i, j, 3) * enPowers(a, i, 1) * raiDeriv(:, a, j)
          tmpd = (6 * eePowers(i, j, 1) * rijSquare(i, j) + &
                  3 * eePowers(i, j, 2) * rijLapl(i, j)) * tmp
          nftermsl(1, 5) = tmpd + enPowers(a, j, 1) * &
                           (6 * eePowers(i, j, 2) * rComb(1) + &
                                eePowers(i, j, 3) * raiLapl(a, i))
          nftermsl(2, 5) = tmpd + enPowers(a, i, 1) * &
                           (6 * eePowers(i, j, 2) * rComb(2) + &
                                eePowers(i, j, 3) * raiLapl(a, j))

          do t = 1, 4
            terms(11+t) = nfterms(t) - nfterms(5)

            termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
            termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
          enddo
        endif

        if(fmax > 5) then
          ! ri^2 rj^4 + ri^4 rj^2
          terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                      enPowers(a, i, 2) * enPowers(a, j, 4)

          tmpd = 4 * enPowers(a, i, 3) * enPowers(a, j, 2) + &
                 2 * enPowers(a, i, 1) * enPowers(a, j, 4)
          termsn(1:3, 16) = tmpd * raiDeriv(:, a, i)
          tmp = 12 * enPowers(a, i, 2) * enPowers(a, j, 2)
          termsl(1, 16) = tmpd * raiLapl(a, i) + &
                          (tmp + 2 * enPowers(a, j, 4)) * raiSquare(a, i)

          tmpd = 4 * enPowers(a, j, 3) * enPowers(a, i, 2) + &
                 2 * enPowers(a, j, 1) * enPowers(a, i, 4)
          termsn(4:6, 16) = tmpd * raiDeriv(:, a, j)
          termsl(2, 16) = tmpd * raiLapl(a, j) + &
                          (tmp + 2 * enPowers(a, i, 4)) * raiSquare(a, j)

          ! 2 ri^3 rj^3
          terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

          tmpd = 6 * enPowers(a, i, 2) * enPowers(a, j, 3)
          termsn(1:3, 17) = tmpd * raiDeriv(:, a, i)
          termsl(1, 17) = 12 * enPowers(a, i, 1) * enPowers(a, j, 3) * raiSquare(a, i) + &
                          tmpd * raiLapl(a, i)

          tmpd = 6 * enPowers(a, j, 2) * enPowers(a, i, 3)
          termsn(4:6, 17) = tmpd * raiDeriv(:, a, j)
          termsl(2, 17) = 12 * enPowers(a, j, 1) * enPowers(a, i, 3) * raiSquare(a, j) + &
                          tmpd * raiLapl(a, j)

          ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
          tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 3) * enPowers(a, j, 2) - &
                enPowers(a, i, 2) * enPowers(a, j, 3)
          terms(18) = eePowers(i, j, 1) * tmp

          tmpd = 5 * enPowers(a, i, 4) - 3 * enPowers(a, i, 2) * enPowers(a, j, 2) - &
                                         2 * enPowers(a, i, 1) * enPowers(a, j, 3)
          tmpn(1:3) = rijDeriv(:, i, j) * tmp
          termsn(1:3, 18) = tmpn(1:3) + eePowers(i, j, 1) * tmpd * raiDeriv(:, a, i)
          tmpl = rijLapl(i, j) * tmp
          termsl(1, 18) = tmpl + 2 * rComb(1) * tmpd + &
                          eePowers(i, j, 1) * (tmpd * raiLapl(a, i) + &
                           (20 * enPowers(a, i, 3) - 6 * enPowers(a, i, 1) * enPowers(a, j, 2) -&
                             2 * enPowers(a, j, 3)) * raiSquare(a, i))

          tmpd = 5 * enPowers(a, j, 4) - 3 * enPowers(a, j, 2) * enPowers(a, i, 2) - &
                                         2 * enPowers(a, j, 1) * enPowers(a, i, 3)
          termsn(4:6, 18) = -tmpn(1:3) + eePowers(i, j, 1) * tmpd * raiDeriv(:, a, j)
          termsl(2, 18) = tmpl + 2 * rComb(2) * tmpd + &
                          eePowers(i, j, 1) * (tmpd * raiLapl(a, j) + &
                           (20 * enPowers(a, j, 3) - 6 * enPowers(a, j, 1) * enPowers(a, i, 2) -&
                             2 * enPowers(a, i, 3)) * raiSquare(a, j))

          ! rij^4 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(19) = eePowers(i, j, 4) * tmp

          tmpd = 2 * enPowers(a, i, 1)
          tmpn(1:3) = 4 * eePowers(i, j, 3) * rijDeriv(:, i, j) * tmp
          termsn(1:3, 19) = tmpn(1:3) + eePowers(i, j, 4) * tmpd * raiDeriv(:, a, i)
          tmpl = 4 * tmp * (3 * eePowers(i, j, 2) * rijSquare(i, j) + &
                                eePowers(i, j, 3) * rijLapl(i, j))
          termsl(1, 19) = tmpl + eePowers(i, j, 4) * 2 * raiSquare(a, i) + &
                                 (8 * eePowers(i, j, 3) * rComb(1) + &
                                      eePowers(i, j, 4) * raiLapl(a, i)) * tmpd
          tmpd = 2 * enPowers(a, j, 1)
          termsn(4:6, 19) = -tmpn(1:3) + eePowers(i, j, 4) * tmpd * raiDeriv(:, a, j)
          termsl(2, 19) = tmpl + eePowers(i, j, 4) * 2 * raiSquare(a, j) + &
                                 (8 * eePowers(i, j, 3) * rComb(2) + &
                                      eePowers(i, j, 4) * raiLapl(a, j)) * tmpd

          ! 2 rij^2 ri^2 rj^2
          tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
          terms(20) = 2 * eePowers(i, j, 2) * tmp

          tmpd = 2 * enPowers(a, i, 1) * enPowers(a, j, 2)
          tmpn(1:3) = 4 * eePowers(i, j, 1) * rijDeriv(:, i, j) * tmp
          termsn(1:3, 20) = tmpn(1:3) + 2 * eePowers(i, j, 2) * tmpd * raiDeriv(:, a, i)
          tmpl = 4 * tmp * (rijSquare(i, j) + eePowers(i, j, 1) * rijLapl(i, j))
          termsl(1, 20) = tmpl + 4 * eePowers(i, j, 2) * enPowers(a, j, 2) * raiSquare(a, i) + &
                                 (8 * eePowers(i, j, 1) * rComb(1) + &
                                  2 * eePowers(i, j, 2) * raiLapl(a, i)) * tmpd
          tmpd = 2 * enPowers(a, j, 1) * enPowers(a, i, 2)
          termsn(4:6, 20) = -tmpn(1:3) + 2 * eePowers(i, j, 2) * tmpd * raiDeriv(:, a, j)
          termsl(2, 20) = tmpl + 4 * eePowers(i, j, 2) * enPowers(a, i, 2) * raiSquare(a, j) + &
                                 (8 * eePowers(i, j, 1) * rComb(2) + &
                                  2 * eePowers(i, j, 2) * raiLapl(a, j)) * tmpd

          ! rij^3 (ri^3 + rj^3)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
          terms(21) = eePowers(i, j, 3) * tmp

          tmpd = 3 * enPowers(a, i, 2)
          tmpn(1:3) = 3 * eePowers(i, j, 2) * rijDeriv(:, i, j) * tmp
          termsn(1:3, 21) = tmpn(1:3) + eePowers(i, j, 3) * tmpd * raiDeriv(:, a, i)
          tmpl = 3 * tmp * (2 * eePowers(i, j, 1) * rijSquare(i, j) + eePowers(i, j, 2) * rijLapl(i, j))
          termsl(1, 21) = tmpl +  6 * eePowers(i, j, 3) * enPowers(a, i, 1) * raiSquare(a, i) + &
                                 (6 * eePowers(i, j, 2) * rComb(1) + &
                                      eePowers(i, j, 3) * raiLapl(a, i)) * tmpd
          tmpd = 3 * enPowers(a, j, 2)
          termsn(4:6, 21) = -tmpn(1:3) + eePowers(i, j, 3) * tmpd * raiDeriv(:, a, j)
          termsl(2, 21) = tmpl +  6 * eePowers(i, j, 3) * enPowers(a, j, 1) * raiSquare(a, j) + &
                                 (6 * eePowers(i, j, 2) * rComb(2) + &
                                      eePowers(i, j, 3) * raiLapl(a, j)) * tmpd

          ! rij^2 (ri^4 + rj^4)
          tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
          terms(22) = eePowers(i, j, 2) * tmp

          tmpd = 4 * enPowers(a, i, 3)
          tmpn(1:3) = 2 * eePowers(i, j, 1) * rijDeriv(:, i, j) * tmp
          termsn(1:3, 22) = tmpn(1:3) + eePowers(i, j, 2) * tmpd * raiDeriv(:, a, i)
          tmpl = 2 * tmp * (rijSquare(i, j) + eePowers(i, j, 1) * rijLapl(i, j))
          termsl(1, 22) = tmpl + (4 * eePowers(i, j, 1) * rComb(1) + &
                                      eePowers(i, j, 2) * raiLapl(a, i)) * tmpd + &
                                 12 * eePowers(i, j, 2) * raiSquare(a, i) * enPowers(a, i, 2)

          tmpd = 4 * enPowers(a, j, 3)
          termsn(4:6, 22) = -tmpn(1:3) + eePowers(i, j, 2) * tmpd * raiDeriv(:, a, j)
          termsl(2, 22) = tmpl + (4 * eePowers(i, j, 1) * rComb(2) + &
                                      eePowers(i, j, 2) * raiLapl(a, j)) * tmpd + &
                                 12 * eePowers(i, j, 2) * raiSquare(a, j) * enPowers(a, j, 2)

          ! ri^5 rj + ri rj^5
          nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 5)

          nftermsd(1, 1) = 5 * enPowers(a, i, 4) * enPowers(a, j, 1) + &
                                                   enPowers(a, j, 5)
          nftermsd(2, 1) = 5 * enPowers(a, j, 4) * enPowers(a, i, 1) + &
                                                   enPowers(a, i, 5)
          nftermsn(1:3, 1) = nftermsd(1, 1) * raiDeriv(:, a, i)
          nftermsn(4:6, 1) = nftermsd(2, 1) * raiDeriv(:, a, j)
          nftermsl(1, 1) = nftermsd(1, 1) * raiLapl(a, i) + &
                           (20 * enPowers(a, i, 3) * enPowers(a, j, 1)) * raiSquare(a, i)
          nftermsl(2, 1) = nftermsd(2, 1) * raiLapl(a, j) + &
                           (20 * enPowers(a, j, 3) * enPowers(a, i, 1)) * raiSquare(a, j)

          ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
          tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 4) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 4)
          nfterms(2) = - eePowers(i, j, 1) * tmp

          nftermsd(1, 2) = 5 * enPowers(a, i, 4) - enPowers(a, j, 4) - &
                           4 * enPowers(a, i, 3) * enPowers(a, j, 1)
          nftermsd(2, 2) = 5 * enPowers(a, j, 4) - enPowers(a, i, 4) - &
                           4 * enPowers(a, j, 3) * enPowers(a, i, 1)
          tmpn(1:3) = rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 2) = -tmpn(1:3) - eePowers(i, j, 1) * nftermsd(1, 2) * raiDeriv(:, a, i)
          nftermsn(4:6, 2) =  tmpn(1:3) - eePowers(i, j, 1) * nftermsd(2, 2) * raiDeriv(:, a, j)
          tmpd = -rijLapl(i, j) * tmp
          nftermsl(1, 2) = tmpd - (2 * rComb(1) + eePowers(i, j, 1) * raiLapl(a, i)) * nftermsd(1, 2) - &
                                  (20 * enPowers(a, i, 3) - 12 * enPowers(a, i, 2) * enPowers(a, j, 1)) * &
                                   eePowers(i, j, 1) * raiSquare(a, i)
          nftermsl(2, 2) = tmpd - (2 * rComb(2) + eePowers(i, j, 1) * raiLapl(a, j)) * nftermsd(2, 2) - &
                                  (20 * enPowers(a, j, 3) - 12 * enPowers(a, j, 2) * enPowers(a, i, 1)) * &
                                   eePowers(i, j, 1) * raiSquare(a, j)

          ! rij^5 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 5) * tmp

          tmpn(1:3) = 5 * eePowers(i, j, 4) * rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 3) =  tmpn(1:3) + eePowers(i, j, 5) * raiDeriv(:, a, i)
          nftermsn(4:6, 3) = -tmpn(1:3) + eePowers(i, j, 5) * raiDeriv(:, a, j)
          tmpl = 5 * tmp * (4 * eePowers(i, j, 3) * rijSquare(i, j) + &
                                eePowers(i, j, 4) * rijLapl(i, j))
          nftermsl(1, 3) = tmpl + 10 * eePowers(i, j, 4) * rComb(1) + &
                                       eePowers(i, j, 5) * raiLapl(a, i)
          nftermsl(2, 3) = tmpl + 10 * eePowers(i, j, 4) * rComb(2) + &
                                       eePowers(i, j, 5) * raiLapl(a, j)

          ! rij^2 (ri rj^3 + ri^3 rj)
          tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 3)
          nfterms(4) = eePowers(i, j, 2) * tmp

          nftermsd(1, 4) = 3 * enPowers(a, i, 2) * enPowers(a, j, 1) + enPowers(a, j, 3)
          nftermsd(2, 4) = 3 * enPowers(a, j, 2) * enPowers(a, i, 1) + enPowers(a, i, 3)
          tmpn(1:3) = 2 * eePowers(i, j, 1) * rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 4) =  tmpn(1:3) + eePowers(i, j, 2) * nftermsd(1, 4) * raiDeriv(:, a, i)
          nftermsn(4:6, 4) = -tmpn(1:3) + eePowers(i, j, 2) * nftermsd(2, 4) * raiDeriv(:, a, j)
          tmpl = 2 * tmp * (rijSquare(i, j) + eePowers(i, j, 1) * rijLapl(i, j))
          tmp = 6 * eePowers(i, j, 2) * enPowers(a, i, 1) * enPowers(a, j, 1)
          nftermsl(1, 4) = tmpl + tmp * raiSquare(a, i) + &
                                  (4 * eePowers(i, j, 1) * rComb(1) + &
                                       eePowers(i, j, 2) * raiLapl(a, i)) * nftermsd(1, 4)
          nftermsl(2, 4) = tmpl + tmp * raiSquare(a, j) + &
                                  (4 * eePowers(i, j, 1) * rComb(2) + &
                                       eePowers(i, j, 2) * raiLapl(a, j)) * nftermsd(2, 4)

          ! rij^3 (ri rj^2 + ri^2 rj)
          tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(5) = eePowers(i, j, 3) * tmp

          nftermsd(1, 5) = 2 * enPowers(a, i, 1) * enPowers(a, j, 1) + enPowers(a, j, 2)
          nftermsd(2, 5) = 2 * enPowers(a, i, 1) * enPowers(a, j, 1) + enPowers(a, i, 2)
          tmpn(1:3) = 3 * eePowers(i, j, 2) * rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 5) =  tmpn(1:3) + eePowers(i, j, 3) * nftermsd(1, 5) * raiDeriv(:, a, i)
          nftermsn(4:6, 5) = -tmpn(1:3) + eePowers(i, j, 3) * nftermsd(2, 5) * raiDeriv(:, a, j)
          tmpl = 3 * tmp * (2 * eePowers(i, j, 1) * rijSquare(i, j) + eePowers(i, j, 2) * rijLapl(i, j))
          nftermsl(1, 5) = tmpl + (6 * eePowers(i, j, 2) * rComb(1) + &
                                       eePowers(i, j, 3) * raiLapl(a, i)) * nftermsd(1, 5) + &
                                  2 * eePowers(i, j, 3) * raiSquare(a, i) * enPowers(a, j, 1)
          nftermsl(2, 5) = tmpl + (6 * eePowers(i, j, 2) * rComb(2) + &
                                       eePowers(i, j, 3) * raiLapl(a, j)) * nftermsd(2, 5) + &
                                  2 * eePowers(i, j, 3) * raiSquare(a, j) * enPowers(a, i, 1)

          ! rij^4 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(6) = eePowers(i, j, 4) * tmp

          tmpn(1:3) = 4 * eePowers(i, j, 3) * rijDeriv(:, i, j) * tmp
          nftermsn(1:3, 6) =  tmpn(1:3) + eePowers(i, j, 4) * enPowers(a, j, 1) * raiDeriv(:, a, i)
          nftermsn(4:6, 6) = -tmpn(1:3) + eePowers(i, j, 4) * enPowers(a, i, 1) * raiDeriv(:, a, j)
          tmpl = 4 * tmp * (3 * eePowers(i, j, 2) * rijSquare(i, j) + &
                                eePowers(i, j, 3) * rijLapl(i, j))
          nftermsl(1, 6) = tmpl + (8 * eePowers(i, j, 3) * rComb(1) + &
                                       eePowers(i, j, 4) * raiLapl(a, i)) * enPowers(a, j, 1)
          nftermsl(2, 6) = tmpl + (8 * eePowers(i, j, 3) * rComb(2) + &
                                       eePowers(i, j, 4) * raiLapl(a, j)) * enPowers(a, i, 1)

          do t = 1, 5
            terms(22+t) = nfterms(t) - nfterms(6)

            termsn(:, 22+t) = nftermsn(:, t) - nftermsn(:, 6)
            termsl(:, 22+t) = nftermsl(:, t) - nftermsl(:, 6)
          enddo
        endif

        ii = 3 * i - 2
        jj = 3 * j - 2

        !DIR$ VECTOR
        do g = 1, fpnum
          fTerm = fTerm + gamma(g, c) * terms(g)

          fDeriv(ii) = fDeriv(ii) + gamma(g, c) * termsn(1, g)
          fDeriv(ii+1) = fDeriv(ii+1) + gamma(g, c) * termsn(2, g)
          fDeriv(ii+2) = fDeriv(ii+2) + gamma(g, c) * termsn(3, g)

          fDeriv(jj) = fDeriv(jj) + gamma(g, c) * termsn(4, g)
          fDeriv(jj+1) = fDeriv(jj+1) + gamma(g, c) * termsn(5, g)
          fDeriv(jj+2) = fDeriv(jj+2) + gamma(g, c) * termsn(6, g)

          fLapli(i) = fLapli(i) + gamma(g, c) * termsl(1, g)
          fLapli(j) = fLapli(j) + gamma(g, c) * termsl(2, g)
        enddo

        if(deriveLinParams) then
          do g = 1, fpnum
            m = unum + xnum + (c - 1) * fpnum + g
            uk(m) = uk(m) + terms(g)
            ukgrad(3*i-2:3*i, m) = ukgrad(3*i-2:3*i, m) + termsn(1:3, g)
            ukgrad(3*j-2:3*j, m) = ukgrad(3*j-2:3*j, m) + termsn(4:6, g)
            uklapli(i, m) = uklapli(i, m) + termsl(1, g)
            uklapli(j, m) = uklapli(j, m) + termsl(2, g)
            uklapl(m) = uklapl(m) + termsl(1, g) + termsl(2, g)
          enddo
        endif
      enddo
    enddo
  enddo
end subroutine eenHardcodedAll

subroutine eenHardcodedNoDerivs(fTerm)
  ! value of f term + derivs
  real(r8), intent(inout) :: fTerm

  real(r8) :: rComb(2)

  real(r8) :: tmp
  real(r8) :: terms(fpnum), termsn(6, fpnum), termsl(2, fpnum)
  real(r8) :: nfterms(fmax), nftermsn(6, fmax), nftermsl(2, fmax)
  integer :: a, c, g, i, j, t

  do a = 1, nclast
    c = atoms(a)%sa

    do i = 1, ne
      do j = i + 1, ne
        ! rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
        ! rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

        if(fmax > 2) then
          ! r_i^2 r_j + r_i r_j^2
          nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                       enPowers(a, i, 1) * enPowers(a, j, 2)

          ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                2 * enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(2) = -0.5d0 * eePowers(i, j, 1) * tmp

          ! r_ij^2 (r_i + r_j)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 2) * tmp

          do t = 1, 2
            terms(t) = nfterms(t) - nfterms(3)
          enddo
        endif

        if(fmax > 3) then
          ! 2 ri^2 rj^2
          terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)

          ! rij^2 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(4) = eePowers(i, j, 2) * tmp

          ! ri^3 rj + ri rj^3
          nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 3)

          ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                enPowers(a, i, 2) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(2) = - eePowers(i, j, 1) * tmp

          ! rij^3 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 3) * tmp

          ! rij^2 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(4) = eePowers(i, j, 2) * tmp

          do t = 1, 3
            terms(4+t) = nfterms(t) - nfterms(4)
          enddo
        endif

        if(fmax > 4) then
          ! ri^3 rj^2 + ri^2 rj^3
          terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                     enPowers(a, i, 2) * enPowers(a, j, 3)

          ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
          tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                2 * enPowers(a, i, 2) * enPowers(a, j, 2)
          terms(9)  = eePowers(i, j, 1) * tmp

          ! rij^3 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(10) = eePowers(i, j, 3) * tmp

          ! rij^2 (ri^3 + rj^3)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
          terms(11) = eePowers(i, j, 2) * tmp

          ! ri^4 rj + ri rj^4
          nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 4)
          ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
          tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                enPowers(a, i, 3) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 3)
          nfterms(2) = -eePowers(i, j, 1) * tmp

          ! rij^4 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 4) * tmp

          ! rij^2 (ri^2 rj + ri rj^2)
          tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(4) = eePowers(i, j, 2) * tmp

          ! rij^3 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(5) = eePowers(i, j, 3) * tmp

          do t = 1, 4
            terms(11+t) = nfterms(t) - nfterms(5)
            termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
            termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
          enddo
        endif

        if(fmax > 5) then
          ! ri^2 rj^4 + ri^4 rj^2
          terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                      enPowers(a, i, 2) * enPowers(a, j, 4)

          ! 2 ri^3 rj^3
          terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

          ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
          tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 3) * enPowers(a, j, 2) - &
                enPowers(a, i, 2) * enPowers(a, j, 3)
          terms(18) = eePowers(i, j, 1) * tmp

          ! rij^4 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(19) = eePowers(i, j, 4) * tmp

          ! 2 rij^2 ri^2 rj^2
          tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
          terms(20) = 2 * eePowers(i, j, 2) * tmp

          ! rij^3 (ri^3 + rj^3)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
          terms(21) = eePowers(i, j, 3) * tmp

          ! rij^2 (ri^4 + rj^4)
          tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
          terms(22) = eePowers(i, j, 2) * tmp

          ! ri^5 rj + ri rj^5
          nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 5)

          ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
          tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 4) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 4)
          nfterms(2) = - eePowers(i, j, 1) * tmp

          ! rij^5 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 5) * tmp

          ! rij^2 (ri rj^3 + ri^3 rj)
          tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 3)
          nfterms(4) = eePowers(i, j, 2) * tmp

          ! rij^3 (ri rj^2 + ri^2 rj)
          tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(5) = eePowers(i, j, 3) * tmp

          ! rij^4 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(6) = eePowers(i, j, 4) * tmp

          do t = 1, 5
            terms(22+t) = nfterms(t) - nfterms(6)
          enddo
        endif

        do g = 1, fpnum
          fTerm = fTerm + gamma(g, c) * terms(g)
        enddo
      enddo
    enddo
  enddo
end subroutine eenHardcodedNoDerivs
subroutine eenHardcodedWithUk(fTerm, deriveLinParams)
  ! value of f term + derivs
  real(r8), intent(inout) :: fTerm
  logical, intent(in) :: deriveLinParams

  real(r8) :: rComb(2)

  real(r8) :: tmp
  real(r8) :: terms(fpnum), termsn(6, fpnum), termsl(2, fpnum)
  real(r8) :: nfterms(fmax), nftermsn(6, fmax), nftermsl(2, fmax)
  integer :: a, c, g, i, j, m, t

  do a = 1, nclast
    c = atoms(a)%sa

    do i = 1, ne
      do j = i + 1, ne
        ! !DIR$ vector
        ! rComb(1) =  sum(raiDeriv(:, a, i) * rijDeriv(:, i, j))
        ! !DIR$ vector
        ! rComb(2) = -sum(raiDeriv(:, a, j) * rijDeriv(:, i, j))

        if(fmax > 2) then
          ! r_i^2 r_j + r_i r_j^2
          nfterms(1) = enPowers(a, i, 2) * enPowers(a, j, 1)  + &
                       enPowers(a, i, 1) * enPowers(a, j, 2)

          ! -0.5 * r_ij (r_i^2 + r_j^2 - 2 r_i r_j)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2) - &
                2 * enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(2) = -0.5d0 * eePowers(i, j, 1) * tmp
          ! r_ij^2 (r_i + r_j)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 2) * tmp
          do t = 1, 2
            terms(t) = nfterms(t) - nfterms(3)
          enddo
        endif

        if(fmax > 3) then
          ! 2 ri^2 rj^2
          terms(3) = 2 * enPowers(a, i, 2) * enPowers(a, j, 2)
          ! rij^2 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(4) = eePowers(i, j, 2) * tmp

          ! ri^3 rj + ri rj^3
          nfterms(1) = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 3)
          ! -rij (ri^3 + rj^3 - ri^2 rj - ri rj^2)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3) - &
                enPowers(a, i, 2) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(2) = - eePowers(i, j, 1) * tmp

          ! rij^3 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 3) * tmp
          ! rij^2 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(4) = eePowers(i, j, 2) * tmp

          do t = 1, 3
            terms(4+t) = nfterms(t) - nfterms(4)
          enddo
        endif

        if(fmax > 4) then
          ! ri^3 rj^2 + ri^2 rj^3
          terms(8) = enPowers(a, i, 3) * enPowers(a, j, 2) + &
                     enPowers(a, i, 2) * enPowers(a, j, 3)

          ! rij (ri^4 + rj^4 - 2 ri^2 rj^2)
          tmp =     enPowers(a, i, 4) + enPowers(a, j, 4) - &
                2 * enPowers(a, i, 2) * enPowers(a, j, 2)
          terms(9)  = eePowers(i, j, 1) * tmp
          ! rij^3 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(10) = eePowers(i, j, 3) * tmp
          ! rij^2 (ri^3 + rj^3)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
          terms(11) = eePowers(i, j, 2) * tmp
          ! ri^4 rj + ri rj^4
          nfterms(1) = enPowers(a, i, 4) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 4)
          ! -rij (ri^4 + rj^4 - ri^3 rj - ri rj^3)
          tmp = enPowers(a, i, 4) + enPowers(a, j, 4) - &
                enPowers(a, i, 3) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 3)
          nfterms(2) = -eePowers(i, j, 1) * tmp
          ! rij^4 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 4) * tmp
          ! rij^2 (ri^2 rj + ri rj^2)
          tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(4) = eePowers(i, j, 2) * tmp
          ! rij^3 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(5) = eePowers(i, j, 3) * tmp

          do t = 1, 4
            terms(11+t) = nfterms(t) - nfterms(5)

            termsn(:, 11+t) = nftermsn(:, t) - nftermsn(:, 5)
            termsl(:, 11+t) = nftermsl(:, t) - nftermsl(:, 5)
          enddo
        endif

        if(fmax > 5) then
          ! ri^2 rj^4 + ri^4 rj^2
          terms(16) = enPowers(a, i, 4) * enPowers(a, j, 2) + &
                      enPowers(a, i, 2) * enPowers(a, j, 4)
          ! 2 ri^3 rj^3
          terms(17) = 2 * enPowers(a, i, 3) * enPowers(a, j, 3)

          ! rij (ri^5 + rj^5 - ri^3 rj^2 - ri^2 rj^3)
          tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 3) * enPowers(a, j, 2) - &
                enPowers(a, i, 2) * enPowers(a, j, 3)
          terms(18) = eePowers(i, j, 1) * tmp

          ! rij^4 (ri^2 + rj^2)
          tmp = enPowers(a, i, 2) + enPowers(a, j, 2)
          terms(19) = eePowers(i, j, 4) * tmp
          ! 2 rij^2 ri^2 rj^2
          tmp = enPowers(a, i, 2) * enPowers(a, j, 2)
          terms(20) = 2 * eePowers(i, j, 2) * tmp

          ! rij^3 (ri^3 + rj^3)
          tmp = enPowers(a, i, 3) + enPowers(a, j, 3)
          terms(21) = eePowers(i, j, 3) * tmp
          ! rij^2 (ri^4 + rj^4)
          tmp = enPowers(a, i, 4) + enPowers(a, j, 4)
          terms(22) = eePowers(i, j, 2) * tmp

          ! ri^5 rj + ri rj^5
          nfterms(1) = enPowers(a, i, 5) * enPowers(a, j, 1) + &
                       enPowers(a, i, 1) * enPowers(a, j, 5)
          ! -rij (ri^5 + rj^5 - ri rj^4 - ri^4 rj)
          tmp = enPowers(a, i, 5) + enPowers(a, j, 5) - &
                enPowers(a, i, 4) * enPowers(a, j, 1) - &
                enPowers(a, i, 1) * enPowers(a, j, 4)
          nfterms(2) = - eePowers(i, j, 1) * tmp
          ! rij^5 (ri + rj)
          tmp = enPowers(a, i, 1) + enPowers(a, j, 1)
          nfterms(3) = eePowers(i, j, 5) * tmp
          ! rij^2 (ri rj^3 + ri^3 rj)
          tmp = enPowers(a, i, 3) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 3)
          nfterms(4) = eePowers(i, j, 2) * tmp
          ! rij^3 (ri rj^2 + ri^2 rj)
          tmp = enPowers(a, i, 2) * enPowers(a, j, 1) + &
                enPowers(a, i, 1) * enPowers(a, j, 2)
          nfterms(5) = eePowers(i, j, 3) * tmp
          ! rij^4 ri rj
          tmp = enPowers(a, i, 1) * enPowers(a, j, 1)
          nfterms(6) = eePowers(i, j, 4) * tmp

          do t = 1, 5
            terms(22+t) = nfterms(t) - nfterms(6)
          enddo
        endif

        !DIR$ VECTOR
        do g = 1, fpnum
          fTerm = fTerm + gamma(g, c) * terms(g)
        enddo

        if(deriveLinParams) then
          do g = 1, fpnum
            m = unum + xnum + (c - 1) * fpnum + g
            uk(m) = uk(m) + terms(g)
          enddo
        endif
      enddo
    enddo
  enddo
end subroutine eenHardcodedWithUk


subroutine allocateTerms()
  integer :: stat

  allocatedTerms = .true.
  allocate(eePowers(ne, ne, -2:max(umax, fmax, 1)), &
           enPowers(nclast, ne, -2:max(xmax, fmax, 1)), &
           rijDeriv(3, ne, ne), &
           raiDeriv(3, nclast, ne), &
           rijSquare(ne, ne), &
           raiSquare(nclast, ne), &
           rijLapl(ne, ne), &
           raiLapl(nclast, ne), &
           stat=stat)
  if(stat /= 0) call abortp("Jastrow allocation failed")
end subroutine allocateTerms


subroutine deallocateTerms()
  if(.not. allocatedTerms) return
  allocatedTerms = .false.
  deallocate(eePowers, enPowers, rijDeriv, raiDeriv, rijSquare, raiSquare, rijLapl, raiLapl)
end subroutine deallocateTerms


subroutine setNumParams()
  if(umax > 0 .and. kmax >= umax) then
    unum = umax - 1
  else if(umax == 0) then
    unum = 0
  else
    call abortp("(inputwf): umax out of range")
  endif

  if(xmax > 0 .and. kmax >= xmax) then
    xpnum = xmax - 1
    xnum = ncdiff * xpnum
  else if(xmax == 0) then
    xpnum = 0
    xnum = 0
  else
    call abortp("(inputwf): xmax out of range")
  endif

  if(fmax > 2 .and. kmax >= fmax) then
    fpnum = numEenTerms()
  else if(fmax >= 0) then
    fpnum = 0
  else
    call abortp("(inputwf): fmax out of range")
  endif
  fnum = ncdiff * fpnum

  numParams = unum + xnum + fnum

  anisoJ1_count = 0
  anisoJ2_count = 0

  if (useAOJasTerms) then
    anisoJ1_count = JasAniso%getNumberofParams_J1()
    anisoJ2_count = JasAniso%getNumberofParams_J2()
    numParams = numParams +  anisoJ1_count + anisoJ2_count
  endif
end subroutine setNumParams

integer function numEenTerms()
  numEenTerms = floor((2*(fmax**3)+3*(fmax**2)-2*fmax+3)/24.0d0) + (fmax**2-5*fmax+4)/2
end function numEenTerms

subroutine jas_diffeecusp_ic(pdiff_ee_cusp)
  logical, intent(in) :: pdiff_ee_cusp
  diffeecusp = pdiff_ee_cusp
end subroutine jas_diffeecusp_ic

end module jastrowIC_m
