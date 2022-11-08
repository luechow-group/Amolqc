! Copyright (C) 2018 Christoph Schulte
!
! SPDX-License-Identifier: GPL-3.0-or-later

module jastrowAniso_m

  use kinds_m, only: r8
  use aosData_m
  use rdataUpdate_m
  implicit none

  ! different optmodes
  !TODO Move to seperate module
  ! only linear parameters, default
  integer, parameter :: OPT_LIN = 1
  ! only linear parameters, numerical
  integer, parameter :: OPT_NUM_LIN = 2
  ! only non-linear parameter (always numerical)
  integer, parameter :: OPT_NONLIN = 3
  ! all parameters (linear analytical, nonlinear numerical)
  integer, parameter :: OPT_ALL = 4
  ! all parameters (numerical)
  integer, parameter :: OPT_NUM_ALL = 5

  type, public :: jasAOLines
  ! allows to save input format for writing wf
  character(len=30)       :: aoTypeEntry = ""
  character(len=120)      :: aoLines(5) = ""
  integer                 :: n = 0
  end type jasAOLines

  type, public :: tsymmAO
    !use symmetry for parameter value for aniso jastrow
    integer :: idx1 = 0
    integer :: idx2 = 0
    real(r8) :: factor = 0.0d0
  end type tsymmAO

  type, public :: tmulitmap
    integer, allocatable :: idx(:)
    real(r8), allocatable :: factor(:)
  end type tmulitmap

  type, public :: JasAnisoENType
    private
    type(tsymmAO), allocatable    :: symmAOen(:)
    type(jasAOLines)              :: curJasAOLines
    logical                       :: new_order = .false.
    integer                       :: gnum=0, gnums=0, gnump=0, gnumd=0, gnumf=0
    integer                       :: goptnum=0
    integer, allocatable          :: goptidx(:)
    type(tmulitmap), allocatable  :: goptidx_reverse(:)
    real(r8), allocatable           :: goptfac(:)
    integer, allocatable          :: gidx(:), old_order(:)
    real(r8), allocatable           :: gl(:)
  contains
    procedure :: eval   => anisoen_eval
    procedure :: init   => anisoen_init
    procedure :: read_coefs => anisoen_readcoefs
    procedure :: getParamCount => anisoen_getcount
    procedure :: CalcOnlyU => anisoen_CalcOnlyU
    procedure :: CalcOnlyUandUk => anisoen_CalcOnlyUandUk
    procedure :: UpdateOnlyU => anisoen_UpdateOnlyU
    procedure :: UpdateOnlyUandUk => anisoen_UpdateOnlyUandUk
  end type JasAnisoENType

  type, public :: JasAnisoEENType
    private
    type(tsymmAO), allocatable :: symmAOeen(:)
    type(jasAOLines)      :: currJaseenAOLines
    real(r8), allocatable   :: geen(:)
    integer, allocatable  :: eenidx(:),geenidx(:,:)
    integer               :: geennum_t=0, geennum=0, geennums=0
    integer               :: geennump=0, geennumd=0, geennumf=0
    integer                       :: geenoptnum=0
    integer, allocatable          :: geenoptidx(:)
    type(tmulitmap), allocatable  :: geenoptidx_reverse(:)
    real(r8), allocatable           :: geenoptfac(:)
  contains
    procedure :: init   => anisoeen_init
    procedure :: precalc => anisoeen_precalc
    procedure :: read_coefs => anisoeen_readcoefs
    procedure :: getParamCount => anisoeen_getcount
    procedure :: eval   => anisoeen_eval
    procedure :: CalcOnlyU => anisoeen_CalcOnlyU
    procedure :: CalcOnlyUandUk => anisoeen_CalcOnlyUandUk
    procedure :: UpdateOnlyU => anisoeen_UpdateOnlyU
    procedure :: UpdateOnlyUandUk => anisoeen_UpdateOnlyUandUk
  end type JasAnisoEENType

  type, public :: JasAnisoEENNType
    private
    type(tsymmAO), allocatable :: symmAOeenn(:)
    type(jasAOLines)      :: currJaseennAOLines
    real(r8), allocatable   :: geenn(:)
    integer, allocatable  :: eennidx(:),geennidx(:,:)
    integer               :: geennnum=0, geennnums=0, geennnump=0, geennnumd=0
    integer               :: geennnumf=0, geennnum_t=0
    integer                       :: geennoptnum=0
    integer, allocatable          :: geennoptidx(:)
    type(tmulitmap), allocatable  :: geennoptidx_reverse(:)
    real(r8), allocatable           :: geennoptfac(:)
  contains
    procedure :: init   => anisoeenn_init
    procedure :: precalc => anisoeenn_precalc
    procedure :: read_coefs => anisoeenn_readcoefs
    procedure :: getParamCount => anisoeenn_getcount
    procedure :: eval   => anisoeenn_eval
    procedure :: CalcOnlyU => anisoeenn_CalcOnlyU
    procedure :: CalcOnlyUandUk => anisoeenn_CalcOnlyUandUk
    procedure :: UpdateOnlyU => anisoeenn_UpdateOnlyU
    procedure :: UpdateOnlyUandUk => anisoeenn_UpdateOnlyUandUk
  end type JasAnisoEENNType


  type, public :: JasAnisoType
    private
    type(JasAnisoENType)   :: en
    type(JasAnisoEENType)  :: een
    type(JasAnisoEENNType) :: eenn
    logical :: m_useAOJasTerms = .false.
  contains
    procedure :: init   => anisoall_init
    procedure :: destroy => anisoall_destroy
    procedure :: add   => anisoall_add
    procedure :: read_coefs => anisoall_readcoefs
    procedure :: zero_coefs => anisoall_zerocoefs
    procedure :: read_symm => anisoall_readsymm
    procedure :: output => anisoall_output
    procedure :: output_param_new => anisoall_output_param_new
    procedure :: output_header_new => anisoall_output_header_new
    procedure :: output_header_new_lines => anisoall_output_header_new_lines
    procedure :: shortoutput => anisoall_shortoutput
    procedure :: getParamCount => anisoall_getcount
    procedure :: getParamVector => anisoall_getParamVector
    procedure :: setParamVector => anisoall_setParamVector
    procedure :: CalcOnlyU => anisoall_CalcOnlyU
    procedure :: CalcOnlyUandUk => anisoall_CalcOnlyUandUk
    procedure :: UpdateOnlyU => anisoall_UpdateOnlyU
    procedure :: UpdateOnlyUandUk => anisoall_UpdateOnlyUandUk
    procedure :: eval   => anisoall_eval
    procedure :: getNumberofParams => anisoall_getNumberofParams
    procedure :: getNumberofParams_J1 => anisoall_getNumberofParams_J1
    procedure :: getNumberofParams_J2 => anisoall_getNumberofParams_J2
  end type JasAnisoType

contains
 ! ------------- DESTROY -------------
  subroutine anisoall_destroy(this)
    class(JasAnisoType), intent(inout) :: this
    !ToDo call destroy of each individual aniso type
    this%en%gnum = 0
    this%een%geennum = 0
    this%eenn%geennnum = 0
    if (allocated(this%en%gidx)) deallocate(this%en%gidx)
    if (allocated(this%en%goptidx)) deallocate(this%en%goptidx)
    if (allocated(this%en%goptfac)) deallocate(this%en%goptfac)
    if (allocated(this%en%goptidx_reverse)) deallocate(this%en%goptidx_reverse)
    if (allocated(this%en%gl)) deallocate(this%en%gl)
    if (allocated(this%en%symmAOen)) deallocate(this%en%symmAOen)

    if (allocated(this%een%geenidx)) deallocate(this%een%geenidx)
    if (allocated(this%een%geenoptidx)) deallocate(this%een%geenoptidx)
    if (allocated(this%een%geenoptfac)) deallocate(this%een%geenoptfac)
    if (allocated(this%een%geenoptidx_reverse)) deallocate(this%een%geenoptidx_reverse)
    if (allocated(this%een%eenidx)) deallocate(this%een%eenidx)
    if (allocated(this%een%geen)) deallocate(this%een%geen)
    if (allocated(this%een%symmAOeen)) deallocate(this%een%symmAOeen)

    if (allocated(this%eenn%geennidx)) deallocate(this%eenn%geennidx)
    if (allocated(this%eenn%geennoptidx)) deallocate(this%eenn%geennoptidx)
    if (allocated(this%eenn%geennoptfac)) deallocate(this%eenn%geennoptfac)
    if (allocated(this%eenn%geennoptidx_reverse)) deallocate(this%eenn%geennoptidx_reverse)
    if (allocated(this%eenn%eennidx)) deallocate(this%eenn%eennidx)
    if (allocated(this%eenn%geenn)) deallocate(this%eenn%geenn)
    if (allocated(this%eenn%symmAOeenn)) deallocate(this%eenn%symmAOeenn)
    this%m_useAOJasTerms = .false.
  end subroutine anisoall_destroy

 ! ------------- INIT -------------
  subroutine anisoall_init(this,lines,idx,useAOJasTerms)
    class(JasAnisoType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    integer, intent(inout) :: idx
    logical, intent(out) :: useAOJasTerms
    integer :: pos, nWords
    character(len=10) :: word(20)

    useAOJasTerms = .false.
    call tokenize(lines(idx),word,nWords)
    pos = 4
    do
      if (pos == nWords+1) exit
      if (word(pos) == "ao" .or. word(pos)=="AO") then
        useAOJasTerms = .true.
        call this%en%init(lines,idx,word,nWords,pos)
      else if (word(pos) == "eenao" .or. word(pos)=="eenAO" .or. word(pos)=="EENAO") then
        useAOJasTerms = .true.
        call this%een%init(lines,idx,word,nWords,pos)
      else if (word(pos) == "eennao" .or. word(pos)=="eennAO" .or. word(pos)=="EENNAO") then
        useAOJasTerms = .true.
        call this%eenn%init(lines,idx,word,nWords,pos)
      else
        call abortp("jastrow_aniso init: wrong format in 1st line")
      end if
    end do

    call this%een%precalc()
    call this%eenn%precalc()
    call this%read_symm("symmAO.dat")

    this%m_useAOJasTerms = useAOJasTerms
  end subroutine anisoall_init

  subroutine anisoall_add(this,lines,nl,useAOJasTerms)
    class(JasAnisoType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    logical, intent(out) :: useAOJasTerms
    integer, intent(in)          :: nl      ! actual # of lines
    integer :: pos, nWords, idx
    character(len=10) :: word(20)
    logical :: tchanged_een,tchanged_eenn

    call assert(index(lines(1),'add_aniso_terms')>0,'$change_jastrow: add_aniso_terms must be in first line')
    call assert(nl>=2,'$change_jastrow: add_aniso_terms requires additional lines')
    call tokenize(lines(2),word,nWords)
    if (nWords<3) call abortp('$change_jastrow: add_aniso_terms: wrong format')

    idx=2
    pos = 1
    do
      if (pos == nWords+1) exit
      if (word(pos) == "ao" .or. word(pos)=="AO") then
        useAOJasTerms = .true.
        call this%en%init(lines,idx,word,nWords,pos)
      else if (word(pos) == "eenao" .or. word(pos)=="eenAO" .or. word(pos)=="EENAO") then
        useAOJasTerms = .true.
        tchanged_een = .true.
        call this%een%init(lines,idx,word,nWords,pos)
      else if (word(pos) == "eennao" .or. word(pos)=="eennAO" .or. word(pos)=="EENNAO") then
        useAOJasTerms = .true.
        tchanged_eenn = .true.
        call this%eenn%init(lines,idx,word,nWords,pos)
      else
        call abortp("add_aniso_terms: wrong format in 3rd line")
      end if
    end do

    if (tchanged_een)  call this%een%precalc()
    if (tchanged_eenn) call this%eenn%precalc()
    call this%read_symm("symmAO.dat")

    this%m_useAOJasTerms = useAOJasTerms

  end subroutine anisoall_add

  subroutine anisoen_init(this,lines,idx,word,nWords,pos)
    class(JasAnisoENType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    integer, intent(inout) :: idx
    character(len=10), intent(in) :: word(20)
    integer, intent(out) :: pos
    integer, intent(in) :: nWords
    integer :: tmp

    this%new_order = .false.
    pos = pos+1
    if (word(pos) == "idx") then          ! index mode: give basis function indices
      pos = pos +1
      idx = idx +1
      read(word(pos),*) this%gnum
      pos = pos +1
      call mread_ao_idx(this%gnum,this%gidx,this%gnums,this%gnump,this%gnumd,this%gnumf, &
                       this%curJasAOLines,lines,idx,this%new_order,this%old_order)
      if (this%new_order) then
        write(iul, "(/A)") 'en - AO idx have been resorted.'
      else
        if (allocated(this%old_order)) deallocate(this%old_order)
      endif
    else if (word(pos) == "nuc") then     ! nuc mode: use all basis functions for given l and nuc
      pos = pos +1
      read(word(pos),*) tmp
      pos = pos +1
      call mread_ao_nuc(this%gnum,this%gidx,this%gnums,this%gnump,this%gnumd,this%gnumf, &
                        this%curJasAOLines,tmp,lines,idx)
    else if (word(pos) == "all") then     ! use all basis functions of given l
      pos = pos+1
      call mread_ao_all(this%gnum,this%gidx,this%gnums,this%gnump,this%gnumd,this%gnumf, &
                        this%curJasAOLines,pos,word,nWords)
    else if (.not. word(pos)=="none") then
      call abortp("jastrow_aniso: wrong format for en ao terms")
    end if
    !allocate parameter
    if (allocated(this%gl)) deallocate(this%gl)
    allocate(this%gl(this%gnum))
    this%gl = 0
  end subroutine anisoen_init

  subroutine anisoeen_init(this,lines,idx,word,nWords,pos)
    class(JasAnisoEENType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    integer, intent(inout) :: idx
    character(len=10), intent(in) :: word(20)
    integer, intent(out) :: pos
    integer, intent(in) :: nWords
    logical  :: new_order
    integer, allocatable :: old_order(:)
    integer :: tmp

    new_order = .false.
    pos = pos + 1
    if (word(pos) == "idx") then
      pos = pos +1
      idx = idx+1
      read(word(pos),*) this%geennum_t
      pos = pos +1
      !geennum_t as this is not the toal number of geen-terms
      call mread_ao_idx(this%geennum_t,this%eenidx,this%geennums,this%geennump,this%geennumd, &
                        this%geennumf,this%currJaseenAOLines,lines,idx,new_order,old_order)
      if (new_order) then
        write(iul, "(/A)") 'een - AO idx were not sorted correctly. Correct order of indices is'
        write(iul,*) this%eenidx
        write(iul,*) "Parameter are used as if given in correct order."
        write(iul,*) "Check calculation and resort parameter if needed!"
        if (allocated(old_order)) deallocate(old_order)
      endif
    else if (word(pos) == "nuc") then     ! nuc mode: use all basis functions for given l and nuc
       pos = pos +1
       read(word(pos),*) tmp
       pos = pos +1
      call mread_ao_nuc(this%geennum_t,this%eenidx,this%geennums,this%geennump,this%geennumd, &
                        this%geennumf,this%currJaseenAOLines,tmp,lines,idx)                ! inner function
    else if (word(pos) == "all") then     ! use all basis functions of given l
      pos = pos+1
      call mread_ao_all(this%geennum_t,this%eenidx,this%geennums,this%geennump,this%geennumd, &
                        this%geennumf,this%currJaseenAOLines,pos,word,nWords)
    else if (.not. word(pos)=="none") then
      call abortp("jastrow_aniso: wrong format for eenAO terms")
    end if
  end subroutine anisoeen_init

  subroutine anisoeenn_init(this,lines,idx,word,nWords,pos)
    class(JasAnisoEENNType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    integer, intent(inout) :: idx
    character(len=10), intent(in) :: word(20)
    integer, intent(out) :: pos
    integer, intent(in) :: nWords
    logical  :: new_order
    integer, allocatable :: old_order(:)
    integer :: tmp

    new_order = .false.
    pos = pos + 1
    if (word(pos) == "idx") then
      pos = pos +1
      idx = idx+1
      read(word(pos),*) this%geennnum_t
      pos = pos +1
      !geennnum_t as this is not the toal number of geenn-terms
      call mread_ao_idx(this%geennnum_t,this%eennidx,this%geennnums,this%geennnump,this%geennnumd, &
                        this%geennnumf,this%currJaseennAOLines,lines,idx,new_order,old_order)
      if (new_order) then
        write(iul, "(/A)") 'eenn - AO idx were not sorted correctly. Correct order of indices is'
        write(iul,*) this%eennidx
        write(iul,*) "Parameter are used as if given in correct order."
        write(iul,*) "Check calculation and resort parameter if needed!"
        if (allocated(old_order)) deallocate(old_order)
      endif
    else if (word(pos) == "nuc") then     ! nuc mode: use all basis functions for given l and nuc
       pos = pos +1
       read(word(pos),*) tmp
       pos = pos +1
      call mread_ao_nuc(this%geennnum_t,this%eennidx,this%geennnums,this%geennnump,this%geennnumd, &
                        this%geennnumf,this%currJaseennAOLines,tmp,lines,idx)
    else if (word(pos) == "all") then     ! use all basis functions of given l
      pos = pos+1
      call mread_ao_all(this%geennnum_t,this%eennidx,this%geennnums,this%geennnump,this%geennnumd, &
                        this%geennnumf,this%currJaseennAOLines,pos,word,nWords)
    else if (.not. word(pos)=="none") then
      call abortp("jastrow_aniso: wrong format for eennAO terms")
    end if
  end subroutine anisoeenn_init

  ! ------------- PRECALC -------------
  subroutine  anisoeen_precalc(this)
    class(JasAnisoEENType), intent(inout) :: this
    integer :: k,l,m,pos,nucl,nuck
    integer :: temp_list(this%geennum_t*this%geennum_t,2)

    this%geennum = 0
      do k=1,this%geennum_t
        pos = 0
        do m=1, ncenter
          pos = pos + AOIdx%nums(m) + 3*AOIdx%nump(m) + 6*AOIdx%numd(m) + 10*AOIdx%numf(m) + 15*AOIdx%numg(m)
          if ( (this%eenidx(k)-pos ) <= 0 ) then
            nuck = m
            exit
          endif
        end do
        do l=k,this%geennum_t
          !Determine Atom Number for k and l
          pos = 0
          do m=1, ncenter
            pos = pos + AOIdx%nums(m) + 3*AOIdx%nump(m) + 6*AOIdx%numd(m) + 10*AOIdx%numf(m) + 15*AOIdx%numg(m)
            if ( (this%eenidx(l)-pos ) <= 0 ) then
              nucl = m
              exit
            endif
          end do
          if (nucl == nuck) then
            this%geennum = this%geennum +1
            temp_list(this%geennum,1) = this%eenidx(k)
            temp_list(this%geennum,2) = this%eenidx(l)
            !if (MASTER) write(*,*) this%eenidx(k),this%eenidx(l)
          end if
        end do
      end do
      if (allocated(this%geenidx)) deallocate(this%geenidx)
      allocate(this%geenidx(this%geennum,2))
      this%geenidx(1:this%geennum,:) = temp_list(1:this%geennum,:)

      if (allocated(this%geen)) deallocate(this%geen)
      allocate(this%geen(this%geennum))
      this%geen = 0.0d0
  end subroutine  anisoeen_precalc

  subroutine  anisoeenn_precalc(this)
    class(JasAnisoEENNType), intent(inout) :: this
    integer :: k,l,m,pos,nucl,nuck
    integer :: temp_list(this%geennnum_t*this%geennnum_t,2)

    this%geennnum = 0
    do k=1,this%geennnum_t
      pos = 0
      do m=1, ncenter
        pos = pos + AOIdx%nums(m) + 3*AOIdx%nump(m) + 6*AOIdx%numd(m) + 10*AOIdx%numf(m) + 15*AOIdx%numg(m)
        if ( (this%eennidx(k)-pos ) <= 0 ) then
          nuck = m
          exit
        endif
      end do
      do l=k,this%geennnum_t
        !Determine Atom Number for k and l
        pos = 0
        do m=1, ncenter
          pos = pos + AOIdx%nums(m) + 3*AOIdx%nump(m) + 6*AOIdx%numd(m) + 10*AOIdx%numf(m) + 15*AOIdx%numg(m)
          if ( (this%eennidx(l)-pos ) <= 0 ) then
            nucl = m
            exit
          endif
        end do
        if (.not.(nucl == nuck)) then
          this%geennnum = this%geennnum +1
          temp_list(this%geennnum,1) = this%eennidx(k)
          temp_list(this%geennnum,2) = this%eennidx(l)
          !if (MASTER) write(*,*) this%eenidx(k),this%eenidx(l)
        end if
      end do
    end do
    if (allocated(this%geennidx)) deallocate(this%geennidx)
    allocate(this%geennidx(this%geennnum,2))
    this%geennidx(1:this%geennnum,:) = temp_list(1:this%geennnum,:)
    if (allocated(this%geenn)) deallocate(this%geenn)
    allocate(this%geenn(this%geennnum))
    this%geenn = 0.0d0
  end subroutine  anisoeenn_precalc

 ! ------------- READ COEFFICIENTS -------------
  subroutine anisoall_readcoefs(this,lines,idx)
    class(JasAnisoType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    integer, intent(inout) :: idx

    if (this%en%gnum > 0) call this%en%read_coefs(lines,idx)
    if (this%een%geennum > 0) call this%een%read_coefs(lines,idx)
    if (this%eenn%geennnum > 0) call this%eenn%read_coefs(lines,idx)

    ! write(*,*) this%en%gl
  end subroutine anisoall_readcoefs

  subroutine anisoall_zerocoefs(this)
    class(JasAnisoType), intent(inout) :: this

    if (this%en%gnum > 0)  this%en%gl = 0.0d0
    if (this%een%geennum > 0)  this%een%geen= 0.0d0
    if (this%eenn%geennnum > 0)  this%eenn%geenn= 0.0d0
  end subroutine anisoall_zerocoefs

  subroutine anisoall_readsymm(this,filename)
    class(JasAnisoType), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer :: tmp,tmp2,i,j,tmpidx1,tmpidx2,idx
    integer, parameter :: ius = 12
    logical :: foundpair
    integer, parameter :: MAXLINES=100
    logical :: symmAOPresent
    integer :: nsymmLines
    character(len=MAXLEN) :: symmLines(MAXLINES)

    if (allocated(this%en%symmAOen)) deallocate(this%en%symmAOen)
    if (allocated(this%een%symmAOeen)) deallocate(this%een%symmAOeen)
    if (allocated(this%eenn%symmAOeenn)) deallocate(this%eenn%symmAOeenn)

    Inquire(file=filename, exist=symmAOPresent)
    if (symmAOPresent) then
      call readFileParallel(mytid,filename,symmLines,nsymmLines)
      idx=1
      !en symm
      read(symmLines(idx),'(I3)') tmp
      idx = idx+1
      allocate(this%en%symmAOen(tmp))
      do i =1, tmp
        read(symmLines(idx),*) this%en%symmAOen(i)%idx1,this%en%symmAOen(i)%idx2,this%en%symmAOen(i)%factor
        idx = idx+1
      enddo
      !een symm
      read(symmLines(idx),'(I3)') tmp
      idx = idx+1
      allocate(this%een%symmAOeen(tmp))
      do i =1, tmp
        read(symmLines(idx),*) this%een%symmAOeen(i)%idx1,this%een%symmAOeen(i)%idx2,this%een%symmAOeen(i)%factor
        !write(*,*) '--> ',this%een%symmAOeen(i)%idx1,this%een%symmAOeen(i)%idx2,this%een%symmAOeen(i)%factor
        idx = idx+1
      enddo
      !eenn symm
      read(symmLines(idx),'(I3)') tmp
      idx = idx+1
      allocate(this%eenn%symmAOeenn(tmp))
      do i =1, tmp
        read(symmLines(idx),*) this%eenn%symmAOeenn(i)%idx1,this%eenn%symmAOeenn(i)%idx2,this%eenn%symmAOeenn(i)%factor
        idx = idx+1
      enddo
    else
      !do not use symmetry, allocate with size=0 to aviod non-standard calls for list creation
      allocate(this%en%symmAOen(0))
      allocate(this%een%symmAOeen(0))
      allocate(this%eenn%symmAOeenn(0))
    endif

    !Create list of unique coefficients when respecting symmetry
    !aniso en
    this%en%goptnum = 0
    allocate(this%en%goptidx(this%en%gnum),this%en%goptfac(this%en%gnum))
    do i=1, this%en%gnum
      !if (MASTER) write(*,*) 'Searching for pair for ',i,' i.e. ao idx ',this%en%gidx(i)
      foundpair = .false.
      do j=1,(i-1)
        do tmp=1,size(this%en%symmAOen)
          if ((this%en%symmAOen(tmp)%idx1==this%en%gidx(i) .and. this%en%symmAOen(tmp)%idx2==this%en%gidx(j)) .or. &
             ((this%en%symmAOen(tmp)%idx2==this%en%gidx(i) .and. this%en%symmAOen(tmp)%idx1==this%en%gidx(j)))) then
            foundpair = .true.
            exit
          endif
        enddo
        if (foundpair) exit
      enddo
      if (foundpair) then
        this%en%goptidx(i) = this%en%goptidx(j)
        this%en%goptfac(i) = this%en%symmAOen(tmp)%factor
      else
        this%en%goptidx(i) = this%en%goptnum + 1
        this%en%goptfac(i) = 1.0d0
        this%en%goptnum = this%en%goptnum + 1
      endif
    enddo
    !create reverse list
    allocate(this%en%goptidx_reverse(this%en%goptnum))
    do i=1,this%en%goptnum
      allocate(this%en%goptidx_reverse(i)%idx(COUNT(this%en%goptidx == i)))
      allocate(this%en%goptidx_reverse(i)%factor(COUNT(this%en%goptidx == i)))
      tmp = 1
      do j=1, this%en%gnum
        if (this%en%goptidx(j) == i) then
          this%en%goptidx_reverse(i)%idx(tmp) = j
          if (tmp == 1 )then
            this%en%goptidx_reverse(i)%factor(tmp) = 1.0d0
          else
            this%en%goptidx_reverse(i)%factor(tmp) = 1.0d0 / this%en%goptfac(j)
          endif
          tmp = tmp+1
        endif
      enddo
    enddo

    !aniso een AO
    this%een%geenoptnum = 0
    allocate(this%een%geenoptidx(this%een%geennum),this%een%geenoptfac(this%een%geennum))
    do i=1, this%een%geennum
      !if (MASTER) write(*,*) 'Searching for pair for ',i,' i.e. ao idx ',this%een%geenidx(i,1),this%een%geenidx(i,2)
      foundpair = .false.
      do tmp=1, size(this%een%symmAOeen)
        if (this%een%geenidx(i,1) == this%een%symmAOeen(tmp)%idx1 .or. &
            this%een%geenidx(i,1) == this%een%symmAOeen(tmp)%idx2) then
          do tmp2=1, size(this%een%symmAOeen)
            if ( this%een%geenidx(i,2) == this%een%symmAOeen(tmp2)%idx1 .or. &
                 this%een%geenidx(i,2) == this%een%symmAOeen(tmp2)%idx2) then
              ! construct pair
              if (this%een%geenidx(i,1) == this%een%symmAOeen(tmp)%idx1) then
                tmpidx1 = this%een%symmAOeen(tmp)%idx2
              else
                tmpidx1 = this%een%symmAOeen(tmp)%idx1
              endif
              if (this%een%geenidx(i,2) == this%een%symmAOeen(tmp2)%idx1) then
                tmpidx2 = this%een%symmAOeen(tmp2)%idx2
              else
                tmpidx2 = this%een%symmAOeen(tmp2)%idx1
              endif
              !if (MASTER) write(*,*) 'Pair for ',i,' is assumed to consist of. ao idx ',tmpidx1,tmpidx2
              foundpair = .false.
              do j=1,(i-1)
                !search for dependend pair in unique pairs
                if ( (this%een%geenidx(j,1) == tmpidx1 .and. this%een%geenidx(j,2) == tmpidx2) &
                   .or. (this%een%geenidx(j,1) == tmpidx2 .and. this%een%geenidx(j,2) == tmpidx1 ) ) then
                  !first hit is always the one used for optimization
                  !if (MASTER) write(*,*) 'Found pair, ID: ',j
                  foundpair = .true.
                  exit
                endif
              enddo
              if (foundpair) exit
            endif
          enddo
          if (foundpair) exit
        endif
      enddo
      if (foundpair) then
        this%een%geenoptidx(i) = this%een%geenoptidx(j)
        this%een%geenoptfac(i) = this%een%symmAOeen(tmp)%factor * this%een%symmAOeen(tmp2)%factor
      else
        this%een%geenoptidx(i) =  this%een%geenoptnum + 1
        this%een%geenoptfac(i) = 1.0d0
        this%een%geenoptnum = this%een%geenoptnum + 1
      endif
    enddo
    !create reverse list
    allocate(this%een%geenoptidx_reverse(this%een%geenoptnum))
    do i=1,this%een%geenoptnum
      allocate(this%een%geenoptidx_reverse(i)%idx(COUNT(this%een%geenoptidx == i)))
      allocate(this%een%geenoptidx_reverse(i)%factor(COUNT(this%een%geenoptidx == i)))
      tmp = 1
      do j=1, this%een%geennum
        if (this%een%geenoptidx(j) == i) then
          this%een%geenoptidx_reverse(i)%idx(tmp) = j
          if (tmp == 1 )then
            this%een%geenoptidx_reverse(i)%factor(tmp) = 1.0d0
          else
            this%een%geenoptidx_reverse(i)%factor(tmp) = 1.0d0 / this%een%geenoptfac(j)
          endif
          tmp = tmp+1
        endif
      enddo
    enddo

    !aniso eenn AO
    this%eenn%geennoptnum = 0
    allocate(this%eenn%geennoptidx(this%eenn%geennnum),this%eenn%geennoptfac(this%eenn%geennnum))
    do i=1, this%eenn%geennnum
      !if (MASTER) write(*,*) 'Searching for pair for ',i,' i.e. ao idx ',this%eenn%geennidx(i,1),this%eenn%geennidx(i,2)
      foundpair = .false.
      do tmp=1, size(this%eenn%symmAOeenn)
        if (this%eenn%geennidx(i,1) == this%eenn%symmAOeenn(tmp)%idx1 .or. &
            this%eenn%geennidx(i,1) == this%eenn%symmAOeenn(tmp)%idx2) then
          do tmp2=1, size(this%eenn%symmAOeenn)
            if ( this%eenn%geennidx(i,2) == this%eenn%symmAOeenn(tmp2)%idx1 .or. &
                 this%eenn%geennidx(i,2) == this%eenn%symmAOeenn(tmp2)%idx2) then
              ! construct pair
              if (this%eenn%geennidx(i,1) == this%eenn%symmAOeenn(tmp)%idx1) then
                tmpidx1 = this%eenn%symmAOeenn(tmp)%idx2
              else
                tmpidx1 = this%eenn%symmAOeenn(tmp)%idx1
              endif
              if (this%eenn%geennidx(i,2) == this%eenn%symmAOeenn(tmp2)%idx1) then
                tmpidx2 = this%eenn%symmAOeenn(tmp2)%idx2
              else
                tmpidx2 = this%eenn%symmAOeenn(tmp2)%idx1
              endif
              !if (MASTER) write(*,*) 'Pair for ',i,' is assumed to consist of. ao idx ',tmpidx1,tmpidx2
              foundpair = .false.
              do j=1,(i-1)
                !search for dependend pair in unique pairs
                if ( (this%eenn%geennidx(j,1) == tmpidx1 .and. this%eenn%geennidx(j,2) == tmpidx2) &
                   .or. (this%eenn%geennidx(j,1) == tmpidx2 .and. this%eenn%geennidx(j,2) == tmpidx1 ) ) then
                  !first hit is always the one used for optimization
                  !if (MASTER) write(*,*) 'Found pair, ID: ',j
                  foundpair = .true.
                  exit
                endif
              enddo
              if (foundpair) exit
            endif
          enddo
          if (foundpair) exit
        endif
      enddo
      if (foundpair) then
        this%eenn%geennoptidx(i) = this%eenn%geennoptidx(j)
        this%eenn%geennoptfac(i) = this%eenn%symmAOeenn(tmp)%factor * this%eenn%symmAOeenn(tmp2)%factor
      else
        this%eenn%geennoptidx(i) = this%eenn%geennoptnum + 1
        this%eenn%geennoptfac(i) = 1.0d0
        this%eenn%geennoptnum = this%eenn%geennoptnum + 1
      endif
    enddo
    !if (this%eenn%geennnum > 0) then
    !  if (MASTER) then
    !    write(*,*) 'Found ',this%eenn%geennoptnum,' different parameters'
    !    write(*,*) this%eenn%geennoptidx
    !  endif
    !endif
    !create reverse list
    allocate(this%eenn%geennoptidx_reverse(this%eenn%geennoptnum))
    do i=1,this%eenn%geennoptnum
      allocate(this%eenn%geennoptidx_reverse(i)%idx(COUNT(this%eenn%geennoptidx == i)))
      allocate(this%eenn%geennoptidx_reverse(i)%factor(COUNT(this%eenn%geennoptidx == i)))
      tmp = 1
      do j=1, this%eenn%geennnum
        if (this%eenn%geennoptidx(j) == i) then
          this%eenn%geennoptidx_reverse(i)%idx(tmp) = j
          if (tmp == 1 )then
            this%eenn%geennoptidx_reverse(i)%factor(tmp) = 1.0d0
          else
            this%eenn%geennoptidx_reverse(i)%factor(tmp) = 1.0d0 / this%eenn%geennoptfac(j)
          endif
          tmp = tmp+1
        endif
      enddo
    enddo
  end subroutine anisoall_readsymm

  subroutine anisoen_readcoefs(this,lines,idx)
    class(JasAnisoENType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    integer, intent(inout) :: idx
    integer :: m
    real(r8), allocatable :: gl_oldorder(:)

    this%gl = 0.0d0
    do m=1,this%gnum
      read(lines(idx),*) this%gl(m)
      idx = idx + 1
    end do
    if (this%new_order) then
      allocate(gl_oldorder(this%gnum))
      gl_oldorder = this%gl
      do m=1, this%gnum
        this%gl(this%old_order(m)) = gl_oldorder(m)
      enddo
      deallocate(gl_oldorder,this%old_order)
      this%new_order = .false.
    endif
  end subroutine anisoen_readcoefs

  subroutine anisoeen_readcoefs(this,lines,idx)
    class(JasAnisoEENType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    integer, intent(inout) :: idx
    integer :: m

    this%geen = 0.0d0
      do m=1,this%geennum
        read(lines(idx),*) this%geen(m)
        idx = idx + 1
      end do
  end subroutine anisoeen_readcoefs

  subroutine anisoeenn_readcoefs(this,lines,idx)
    class(JasAnisoEENNType), intent(inout) :: this
    character(len=*), intent(in) :: lines(:)
    integer, intent(inout) :: idx
    integer :: m

    this%geenn = 0.0d0
      do m=1,this%geennnum
        read(lines(idx),*) this%geenn(m)
        idx = idx + 1
      end do
  end subroutine anisoeenn_readcoefs

 ! ------------- Calc -------------
  subroutine anisoall_CalcOnlyU(this,Rdu,gTerm)
    class(JasAnisoType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    real(r8), intent(out) :: gTerm

    gTerm = 0.0d0
    call this%en%CalcOnlyU(Rdu,gTerm)
    call this%een%CalcOnlyU(Rdu,gTerm)
    call this%eenn%CalcOnlyU(Rdu,gTerm)
  end subroutine anisoall_CalcOnlyU

  subroutine anisoen_CalcOnlyU(this,Rdu,gTerm)
    class(JasAnisoENType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    real(r8), intent(inout) :: gTerm
    integer nn, i,k

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do i=1,ne
      do k=1,this%gnum
        gTerm = gTerm + this%gl(k) * uao(this%gidx(k),i,nn)
        Rdu%Gi(i) = Rdu%Gi(i) + this%gl(k) * uao(this%gidx(k),i,nn)
      end do
    end do
  end subroutine anisoen_CalcOnlyU

  subroutine anisoeen_CalcOnlyU(this,Rdu,gTerm)
    class(JasAnisoEENType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    real(r8), intent(inout) :: gTerm
    integer nn,i,j,k

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do i=1, ne
      do j=i+1,ne
        do k=1, this%geennum
          !geen(k) symmetric w.r.t. basis function
          gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,1),i,nn) * uao(this%geenidx(k,2),j,nn)
          gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,2),i,nn) * uao(this%geenidx(k,1),j,nn)
          Rdu%Fij(j,i) = Rdu%Fij(j,i)  + this%geen(k) * uao(this%geenidx(k,1),i,nn) * uao(this%geenidx(k,2),j,nn)
          Rdu%Fij(j,i) = Rdu%Fij(j,i)  + this%geen(k) * uao(this%geenidx(k,2),i,nn) * uao(this%geenidx(k,1),j,nn)
        enddo
      enddo
    enddo
  end subroutine anisoeen_CalcOnlyU

  subroutine anisoeenn_CalcOnlyU(this,Rdu,gTerm)
    class(JasAnisoEENNType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    real(r8), intent(inout) :: gTerm
    integer nn, i,j,k

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do i=1, ne
      do j=i+1,ne
        do k=1, this%geennnum
          !geenn(k) symmetric w.r.t. basis function
          gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,1),i,nn) * uao(this%geennidx(k,2),j,nn)
          gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,2),i,nn) * uao(this%geennidx(k,1),j,nn)
          Rdu%Fij(j,i) = Rdu%Fij(j,i)  + this%geenn(k) * uao(this%geennidx(k,1),i,nn) * uao(this%geennidx(k,2),j,nn)
          Rdu%Fij(j,i) = Rdu%Fij(j,i)  + this%geenn(k) * uao(this%geennidx(k,2),i,nn) * uao(this%geennidx(k,1),j,nn)
        enddo
      enddo
    enddo
  end subroutine anisoeenn_CalcOnlyU

! ---
  subroutine anisoall_CalcOnlyUandUk(this,Rdu,gTerm,offsetJ1,offsetJ2)
    class(JasAnisoType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    real(r8), intent(out) :: gTerm
    integer, intent(inout) :: offsetJ1, offsetJ2

    gTerm = 0.0d0
    call this%en%CalcOnlyUandUK(Rdu,gTerm,offsetJ1)
    call this%een%CalcOnlyUandUK(Rdu,gTerm,offsetJ2)
    call this%eenn%CalcOnlyUandUK(Rdu,gTerm,offsetJ2)
  end subroutine anisoall_CalcOnlyUandUk

  subroutine anisoen_CalcOnlyUandUk(this,Rdu,gTerm,offset)
    class(JasAnisoENType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    real(r8), intent(inout) :: gTerm
    integer, intent(inout) :: offset
    integer nn, i,k

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do i=1,ne
      do k=1,this%gnum
        gTerm = gTerm + this%gl(k) * uao(this%gidx(k),i,nn)
        Rdu%Gi(i) = Rdu%Gi(i) + this%gl(k) * uao(this%gidx(k),i,nn)
        Rdu%Gki(offset+this%goptidx(k),i) = Rdu%Gki(offset+this%goptidx(k),i) + uao(this%gidx(k),i,nn)*this%goptfac(k)
      end do
    end do
    offset = offset + this%goptnum
  end subroutine anisoen_CalcOnlyUandUk

  subroutine anisoeen_CalcOnlyUandUk(this,Rdu,gTerm,offset)
    class(JasAnisoEENType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    real(r8), intent(inout) :: gTerm
    integer, intent(inout) :: offset
    integer nn,i,j,k,m

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do i=1, ne
      do j=i+1,ne
        do k=1, this%geennum
          !geen(k) symmetric w.r.t. basis function
          gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,1),i,nn) * uao(this%geenidx(k,2),j,nn)
          gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,2),i,nn) * uao(this%geenidx(k,1),j,nn)
          Rdu%Fij(j,i) = Rdu%Fij(j,i)  + this%geen(k) * uao(this%geenidx(k,1),i,nn) * uao(this%geenidx(k,2),j,nn)
          Rdu%Fij(j,i) = Rdu%Fij(j,i)  + this%geen(k) * uao(this%geenidx(k,2),i,nn) * uao(this%geenidx(k,1),j,nn)
          m = offset + this%geenoptidx(k)
          Rdu%Fijk(j,i,m) = Rdu%Fijk(j,i,m) + uao(this%geenidx(k,1),i,nn) * &
                            uao(this%geenidx(k,2),j,nn)* this%geenoptfac(k)
          Rdu%Fijk(j,i,m) = Rdu%Fijk(j,i,m) + uao(this%geenidx(k,2),i,nn) * &
                            uao(this%geenidx(k,1),j,nn)* this%geenoptfac(k)
        enddo
      enddo
    enddo
    offset = offset + this%geenoptnum
  end subroutine anisoeen_CalcOnlyUandUk

  subroutine anisoeenn_CalcOnlyUandUk(this,Rdu,gTerm,offset)
    class(JasAnisoEENNType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    real(r8), intent(inout) :: gTerm
    integer, intent(inout) :: offset
    integer nn, i,j,k,m

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do i=1, ne
      do j=i+1,ne
        do k=1, this%geennnum
          !geenn(k) symmetric w.r.t. basis function
          gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,1),i,nn) * uao(this%geennidx(k,2),j,nn)
          gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,2),i,nn) * uao(this%geennidx(k,1),j,nn)
          Rdu%Fij(j,i) = Rdu%Fij(j,i)  + this%geenn(k) * uao(this%geennidx(k,1),i,nn) * uao(this%geennidx(k,2),j,nn)
          Rdu%Fij(j,i) = Rdu%Fij(j,i)  + this%geenn(k) * uao(this%geennidx(k,2),i,nn) * uao(this%geennidx(k,1),j,nn)
          m = offset + this%geennoptidx(k)
          Rdu%Fijk(j,i,m) = Rdu%Fijk(j,i,m) + uao(this%geennidx(k,1),i,nn) * &
                            uao(this%geennidx(k,2),j,nn)* this%geennoptfac(k)
          Rdu%Fijk(j,i,m) = Rdu%Fijk(j,i,m) + uao(this%geennidx(k,2),i,nn) * &
                            uao(this%geennidx(k,1),j,nn)* this%geennoptfac(k)
        enddo
      enddo
    enddo
    offset = offset + this%geennoptnum
  end subroutine anisoeenn_CalcOnlyUandUk

! ---
  subroutine anisoall_UpdateOnlyU(this,Rdu,ie)
    class(JasAnisoType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    integer, intent(in) :: ie !update only ie

    call this%en%UpdateOnlyU(Rdu,ie)
    call this%een%UpdateOnlyU(Rdu,ie)
    call this%eenn%UpdateOnlyU(Rdu,ie)
  end subroutine anisoall_UpdateOnlyU

  subroutine anisoen_UpdateOnlyU(this,Rdu,ie)
    class(JasAnisoENType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    integer, intent(in) :: ie !update only ie
    integer nn, k

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do k=1,this%gnum
       !gTerm = gTerm + this%gl(k) * uao(this%gidx(k),i,nn)
       Rdu%Gi(ie) = Rdu%Gi(ie) + this%gl(k) * uao(this%gidx(k),ie,nn)
    end do
  end subroutine anisoen_UpdateOnlyU

  subroutine anisoeen_UpdateOnlyU(this,Rdu,ie)
    class(JasAnisoEENType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    integer, intent(in) :: ie !update only ie
    integer nn,j,k

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do j = 1, ie-1
      do k=1, this%geennum
        !gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,1),ie,nn) * uao(this%geenidx(k,2),j,nn)
        !gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,2),ie,nn) * uao(this%geenidx(k,1),j,nn)
        !geen(k) symmetric w.r.t. basis function
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geen(k) * uao(this%geenidx(k,1),ie,nn) * &
                        uao(this%geenidx(k,2),j,nn)
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geen(k) * uao(this%geenidx(k,2),ie,nn) * &
                        uao(this%geenidx(k,1),j,nn)
      enddo
    enddo
    do j = ie+1, ne
      do k=1, this%geennum
        !gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,1),ie,nn) * uao(this%geenidx(k,2),j,nn)
        !gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,2),ie,nn) * uao(this%geenidx(k,1),j,nn)
        !geen(k) symmetric w.r.t. basis function
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geen(k) * uao(this%geenidx(k,1),ie,nn) * &
                          uao(this%geenidx(k,2),j,nn)
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geen(k) * uao(this%geenidx(k,2),ie,nn) * &
                          uao(this%geenidx(k,1),j,nn)
      enddo
    enddo
  end subroutine anisoeen_UpdateOnlyU

  subroutine anisoeenn_UpdateOnlyU(this,Rdu,ie)
    class(JasAnisoEENNType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculation
    integer, intent(in) :: ie !update only ie
    integer nn,j,k

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do j = 1, ie-1
      do k=1, this%geennnum
        !gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,1),ie,nn) * uao(this%geennidx(k,2),j,nn)
        !gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,2),ie,nn) * uao(this%geennidx(k,1),j,nn)
        !geenn(k) symmetric w.r.t. basis function
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geenn(k) * uao(this%geennidx(k,1),ie,nn) * &
                        uao(this%geennidx(k,2),j,nn)
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geenn(k) * uao(this%geennidx(k,2),ie,nn) * &
                        uao(this%geennidx(k,1),j,nn)
      enddo
    enddo
    do  j = ie+1, ne
      do k=1, this%geennnum
        !gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,1),ie,nn) * uao(this%geennidx(k,2),j,nn)
        !gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,2),ie,nn) * uao(this%geennidx(k,1),j,nn)
        !geenn(k) symmetric w.r.t. basis function
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geenn(k) * uao(this%geennidx(k,1),ie,nn) * &
                        uao(this%geennidx(k,2),j,nn)
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geenn(k) * uao(this%geennidx(k,2),ie,nn) * &
                        uao(this%geennidx(k,1),j,nn)
      enddo
    enddo
  end subroutine anisoeenn_UpdateOnlyU

!----
  subroutine anisoall_UpdateOnlyUandUk(this,Rdu,offsetJ1,offsetJ2,ie)
    class(JasAnisoType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    integer, intent(in) :: ie !update only ie
    integer, intent(inout) :: offsetJ1, offsetJ2

    call this%en%UpdateOnlyUandUK(Rdu,offsetJ1,ie)
    call this%een%UpdateOnlyUandUK(Rdu,offsetJ2,ie)
    call this%eenn%UpdateOnlyUandUK(Rdu,offsetJ2,ie)
  end subroutine anisoall_UpdateOnlyUandUk

  subroutine anisoen_UpdateOnlyUandUk(this,Rdu,offset,ie)
    class(JasAnisoENType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    integer, intent(inout) :: offset
    integer, intent(in) :: ie !update only ie
    integer :: nn,k

    nn = 1
    do k=1,this%gnum
       !gTerm = gTerm + this%gl(k) * uao(this%gidx(k),i,nn)
       Rdu%Gi(ie) = Rdu%Gi(ie) +this%gl(k) * uao(this%gidx(k),ie,nn)
       Rdu%Gki(offset+this%goptidx(k),ie) = Rdu%Gki(offset+this%goptidx(k),ie) + &
                                            uao(this%gidx(k),ie,nn)*this%goptfac(k)
    end do
    offset = offset + this%goptnum
  end subroutine anisoen_UpdateOnlyUandUk

  subroutine anisoeen_UpdateOnlyUandUk(this,Rdu,offset,ie)
    class(JasAnisoEENType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    integer, intent(inout) :: offset
    integer, intent(in) :: ie !update only ie
    integer nn,j,k,m

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do j = 1, ie-1
      do k=1, this%geennum
        !gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,1),ie,nn) * uao(this%geenidx(k,2),j,nn)
        !gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,2),ie,nn) * uao(this%geenidx(k,1),j,nn)
        !geen(k) symmetric w.r.t. basis function
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geen(k) * uao(this%geenidx(k,1),ie,nn) * &
                        uao(this%geenidx(k,2),j,nn)
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geen(k) * uao(this%geenidx(k,2),ie,nn) * &
                        uao(this%geenidx(k,1),j,nn)

        m = offset + this%geenoptidx(k)
        Rdu%Fijk(j,ie,m) = Rdu%Fijk(j,ie,m) + uao(this%geenidx(k,1),ie,nn) * &
                            uao(this%geenidx(k,2),j,nn)* this%geenoptfac(k)
        Rdu%Fijk(j,ie,m) = Rdu%Fijk(j,ie,m) + uao(this%geenidx(k,2),ie,nn) * &
                            uao(this%geenidx(k,1),j,nn)* this%geenoptfac(k)
      enddo
    enddo
    do j = ie+1, ne
      do k=1, this%geennum
        !gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,1),ie,nn) * uao(this%geenidx(k,2),j,nn)
        !gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,2),ie,nn) * uao(this%geenidx(k,1),j,nn)
        !geen(k) symmetric w.r.t. basis function
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geen(k) * uao(this%geenidx(k,1),ie,nn) * &
                        uao(this%geenidx(k,2),j,nn)
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geen(k) * uao(this%geenidx(k,2),ie,nn) * &
                        uao(this%geenidx(k,1),j,nn)

        m = offset + this%geenoptidx(k)
        Rdu%Fijk(j,ie,m) = Rdu%Fijk(j,ie,m) + uao(this%geenidx(k,1),ie,nn) *&
                           uao(this%geenidx(k,2),j,nn)* this%geenoptfac(k)
        Rdu%Fijk(j,ie,m) = Rdu%Fijk(j,ie,m) + uao(this%geenidx(k,2),ie,nn) * &
                           uao(this%geenidx(k,1),j,nn)* this%geenoptfac(k)
      enddo
    enddo
    offset = offset + this%geenoptnum
  end subroutine anisoeen_UpdateOnlyUandUk

  subroutine anisoeenn_UpdateOnlyUandUk(this,Rdu,offset,ie)
    class(JasAnisoEENNType), intent(in) :: this
    type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
    integer, intent(inout) :: offset
    integer, intent(in) :: ie !update only ie
    integer nn,j,k,m

    nn = 1 ! electron configuration idx, only nn=1 allowed
    do j = 1, ie-1
      do k=1, this%geennnum
        !gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,1),ie,nn) * uao(this%geennidx(k,2),j,nn)
        !gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,2),ie,nn) * uao(this%geennidx(k,1),j,nn)
        !geenn(k) symmetric w.r.t. basis function
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geenn(k) * uao(this%geennidx(k,1),ie,nn) * &
                        uao(this%geennidx(k,2),j,nn)
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geenn(k) * uao(this%geennidx(k,2),ie,nn) * &
                        uao(this%geennidx(k,1),j,nn)

        m = offset + this%geennoptidx(k)
        Rdu%Fijk(j,ie,m) = Rdu%Fijk(j,ie,m) + uao(this%geennidx(k,1),ie,nn) * &
                           uao(this%geennidx(k,2),j,nn)* this%geennoptfac(k)
        Rdu%Fijk(j,ie,m) = Rdu%Fijk(j,ie,m) + uao(this%geennidx(k,2),ie,nn) * &
                           uao(this%geennidx(k,1),j,nn)* this%geennoptfac(k)
      enddo
    enddo
    do  j = ie+1, ne
      do k=1, this%geennnum
        !gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,1),ie,nn) * uao(this%geennidx(k,2),j,nn)
        !gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,2),ie,nn) * uao(this%geennidx(k,1),j,nn)
         !geenn(k) symmetric w.r.t. basis function
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geenn(k) * uao(this%geennidx(k,1),ie,nn) * &
                        uao(this%geennidx(k,2),j,nn)
        Rdu%Fij(j,ie) = Rdu%Fij(j,ie)  + this%geenn(k) * uao(this%geennidx(k,2),ie,nn) * &
                        uao(this%geennidx(k,1),j,nn)

        m = offset + this%geennoptidx(k)
        Rdu%Fijk(j,ie,m) = Rdu%Fijk(j,ie,m) + uao(this%geennidx(k,1),ie,nn) * &
                           uao(this%geennidx(k,2),j,nn)* this%geennoptfac(k)
        Rdu%Fijk(j,ie,m) = Rdu%Fijk(j,ie,m) + uao(this%geennidx(k,2),ie,nn) * &
                           uao(this%geennidx(k,1),j,nn)* this%geennoptfac(k)
      enddo
    enddo
    offset = offset + this%geennoptnum
  end subroutine anisoeenn_UpdateOnlyUandUk

! ---
  subroutine anisoall_eval(this,gTerm,gDeriv,gLapli,deriveParams,nn,gk,gkgrad,gklapli,gklapl)
    class(JasAnisoType), intent(inout) :: this
    real(r8), intent(inout) :: gTerm, gDeriv(3*ne), gLapli(ne)
    logical, intent(in) :: deriveParams
    integer, intent(in) :: nn
    real(r8), optional, intent(out) :: gk(:),gkgrad(:,:),gklapli(:,:),gklapl(:)
    integer :: offset

    offset = 0
    gTerm = 0.0d0
    gDeriv = 0.0d0
    gLapli = 0.0d0
    if (deriveParams) then
      gk = 0.0d0
      gkgrad = 0.0d0
      gklapli = 0.0d0
      gklapl = 0.0d0
      call this%en%eval(gTerm,gDeriv,gLapli,deriveParams,offset,nn,gk,gkgrad,gklapli,gklapl)
      call this%een%eval(gTerm,gDeriv,gLapli,deriveParams,offset,nn,gk,gkgrad,gklapli,gklapl)
      call this%eenn%eval(gTerm,gDeriv,gLapli,deriveParams,offset,nn,gk,gkgrad,gklapli,gklapl)
    else
      call this%en%eval(gTerm,gDeriv,gLapli,deriveParams,offset,nn)
      call this%een%eval(gTerm,gDeriv,gLapli,deriveParams,offset,nn)
      call this%eenn%eval(gTerm,gDeriv,gLapli,deriveParams,offset,nn)
    endif
  end subroutine anisoall_eval

  subroutine anisoen_eval(this,gTerm,gDeriv,gLapli,deriveParams,offset,nn,gk,gkgrad,gklapli,gklapl)
    class(JasAnisoENType), intent(inout) :: this
    real(r8), intent(inout) :: gTerm, gDeriv(3*ne), gLapli(ne)
    logical, intent(in) :: deriveParams
    integer, intent(inout) :: offset
    integer, intent(in) :: nn
    real(r8), optional, intent(inout) :: gk(:),gkgrad(:,:),gklapli(:,:),gklapl(:)
    integer :: i,k,m

    do i=1,ne
      do k=1,this%gnum
        gTerm = gTerm + this%gl(k) * uao(this%gidx(k),i,nn)
        gDeriv(3*i-2) = gDeriv(3*i-2) + this%gl(k) * uxao(this%gidx(k),i,nn)
        gDeriv(3*i-1) = gDeriv(3*i-1) + this%gl(k) * uyao(this%gidx(k),i,nn)
        gDeriv(3*i-0) = gDeriv(3*i-0) + this%gl(k) * uzao(this%gidx(k),i,nn)
        gLapli(i) = gLapli(i) + this%gl(k) * u2ao(this%gidx(k),i,nn)
        if (deriveParams) then
          m = offset + this%goptidx(k)
          gk(m) = gk(m) + uao(this%gidx(k),i,nn) * this%goptfac(k)
          gkgrad(3*i-2,m) = gkgrad(3*i-2,m) + uxao(this%gidx(k),i,nn) * this%goptfac(k)
          gkgrad(3*i-1,m) = gkgrad(3*i-1,m) + uyao(this%gidx(k),i,nn) * this%goptfac(k)
          gkgrad(3*i-0,m) = gkgrad(3*i-0,m) + uzao(this%gidx(k),i,nn) * this%goptfac(k)
          gklapli(i,m) = gklapli(i,m) + u2ao(this%gidx(k),i,nn) * this%goptfac(k)
        end if
      end do
    end do

    if (deriveParams) then
      do k=offset+1,offset+this%goptnum
        gklapl(k) = sum(gklapli(:,k))
      end do
      offset = offset+this%goptnum
    end if

  end subroutine anisoen_eval

  subroutine anisoeen_eval(this,gTerm,gDeriv,gLapli,deriveParams,offset,nn,gk,gkgrad,gklapli,gklapl)
    class(JasAnisoEENType), intent(inout) :: this
    real(r8), intent(inout) :: gTerm, gDeriv(3*ne), gLapli(ne)
    logical, intent(in) :: deriveParams
    integer, intent(inout) :: offset
    integer, intent(in) :: nn
    real(r8), optional, intent(inout) :: gk(:),gkgrad(:,:),gklapli(:,:),gklapl(:)
    integer :: i,j,k,m

    do i=1, ne
      do j=i+1,ne
        do k=1, this%geennum
          gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,1),i,nn) * uao(this%geenidx(k,2),j,nn)

          gDeriv(3*i-2) = gDeriv(3*i-2) + this%geen(k) * uao(this%geenidx(k,2),j,nn) * &
                          uxao(this%geenidx(k,1),i,nn)
          gDeriv(3*i-1) = gDeriv(3*i-1) + this%geen(k) * uao(this%geenidx(k,2),j,nn) * &
                          uyao(this%geenidx(k,1),i,nn)
          gDeriv(3*i-0) = gDeriv(3*i-0) + this%geen(k) * uao(this%geenidx(k,2),j,nn) * &
                          uzao(this%geenidx(k,1),i,nn)

          gDeriv(3*j-2) = gDeriv(3*j-2) + this%geen(k) * uao(this%geenidx(k,1),i,nn) * &
                          uxao(this%geenidx(k,2),j,nn)
          gDeriv(3*j-1) = gDeriv(3*j-1) + this%geen(k) * uao(this%geenidx(k,1),i,nn) * &
                          uyao(this%geenidx(k,2),j,nn)
          gDeriv(3*j-0) = gDeriv(3*j-0) + this%geen(k) * uao(this%geenidx(k,1),i,nn) * &
                          uzao(this%geenidx(k,2),j,nn)

          gLapli(i) = gLapli(i) + this%geen(k) * u2ao(this%geenidx(k,1),i,nn) * &
                      uao(this%geenidx(k,2),j,nn)
          gLapli(j) = gLapli(j) + this%geen(k) * u2ao(this%geenidx(k,2),j,nn) * &
                      uao(this%geenidx(k,1),i,nn)

          !this%geen(k) symmetric w.r.t. basis function
          gTerm = gTerm + this%geen(k) * uao(this%geenidx(k,2),i,nn) * uao(this%geenidx(k,1),j,nn)

          gDeriv(3*i-2) = gDeriv(3*i-2) + this%geen(k) * uao(this%geenidx(k,1),j,nn) * &
                          uxao(this%geenidx(k,2),i,nn)
          gDeriv(3*i-1) = gDeriv(3*i-1) + this%geen(k) * uao(this%geenidx(k,1),j,nn) * &
                          uyao(this%geenidx(k,2),i,nn)
          gDeriv(3*i-0) = gDeriv(3*i-0) + this%geen(k) * uao(this%geenidx(k,1),j,nn) * &
                          uzao(this%geenidx(k,2),i,nn)

          gDeriv(3*j-2) = gDeriv(3*j-2) + this%geen(k) * uao(this%geenidx(k,2),i,nn) * &
                          uxao(this%geenidx(k,1),j,nn)
          gDeriv(3*j-1) = gDeriv(3*j-1) + this%geen(k) * uao(this%geenidx(k,2),i,nn) * &
                          uyao(this%geenidx(k,1),j,nn)
          gDeriv(3*j-0) = gDeriv(3*j-0) + this%geen(k) * uao(this%geenidx(k,2),i,nn) * &
                          uzao(this%geenidx(k,1),j,nn)

          gLapli(i) = gLapli(i) + this%geen(k) * u2ao(this%geenidx(k,2),i,nn) * uao(this%geenidx(k,1),j,nn)
          gLapli(j) = gLapli(j) + this%geen(k) * u2ao(this%geenidx(k,1),j,nn) * uao(this%geenidx(k,2),i,nn)

        enddo
        if (deriveParams) then
          do k=1, this%geennum
            m = offset + this%geenoptidx(k)

            gk(m) = gk(m) +  uao(this%geenidx(k,1),i,nn) * uao(this%geenidx(k,2),j,nn) * this%geenoptfac(k)
            gkgrad(3*i-2,m) = gkgrad(3*i-2,m) + uao(this%geenidx(k,2),j,nn) * &
                              uxao(this%geenidx(k,1),i,nn)* this%geenoptfac(k)
            gkgrad(3*i-1,m) = gkgrad(3*i-1,m) + uao(this%geenidx(k,2),j,nn) * &
                              uyao(this%geenidx(k,1),i,nn)* this%geenoptfac(k)
            gkgrad(3*i-0,m) = gkgrad(3*i-0,m) + uao(this%geenidx(k,2),j,nn) * &
                              uzao(this%geenidx(k,1),i,nn)* this%geenoptfac(k)

            gkgrad(3*j-2,m) = gkgrad(3*j-2,m) + uao(this%geenidx(k,1),i,nn) * &
                              uxao(this%geenidx(k,2),j,nn)* this%geenoptfac(k)
            gkgrad(3*j-1,m) = gkgrad(3*j-1,m) + uao(this%geenidx(k,1),i,nn) * &
                              uyao(this%geenidx(k,2),j,nn)* this%geenoptfac(k)
            gkgrad(3*j-0,m) = gkgrad(3*j-0,m) + uao(this%geenidx(k,1),i,nn) * &
                              uzao(this%geenidx(k,2),j,nn)* this%geenoptfac(k)

            gklapli(i,m) = gklapli(i,m) + u2ao(this%geenidx(k,1),i,nn) * &
                            uao(this%geenidx(k,2),j,nn)* this%geenoptfac(k)
            gklapli(j,m) = gklapli(j,m) + u2ao(this%geenidx(k,2),j,nn) * &
                            uao(this%geenidx(k,1),i,nn)* this%geenoptfac(k)
            !geen(k) symmetric w.r.t. basis function
            gk(m) = gk(m) +  uao(this%geenidx(k,2),i,nn) * uao(this%geenidx(k,1),j,nn)* this%geenoptfac(k)
            gkgrad(3*i-2,m) = gkgrad(3*i-2,m) + uao(this%geenidx(k,1),j,nn) * &
                              uxao(this%geenidx(k,2),i,nn)* this%geenoptfac(k)
            gkgrad(3*i-1,m) = gkgrad(3*i-1,m) + uao(this%geenidx(k,1),j,nn) * &
                              uyao(this%geenidx(k,2),i,nn)* this%geenoptfac(k)
            gkgrad(3*i-0,m) = gkgrad(3*i-0,m) + uao(this%geenidx(k,1),j,nn) * &
                              uzao(this%geenidx(k,2),i,nn)* this%geenoptfac(k)

            gkgrad(3*j-2,m) = gkgrad(3*j-2,m) + uao(this%geenidx(k,2),i,nn) * &
                              uxao(this%geenidx(k,1),j,nn)* this%geenoptfac(k)
            gkgrad(3*j-1,m) = gkgrad(3*j-1,m) + uao(this%geenidx(k,2),i,nn) * &
                              uyao(this%geenidx(k,1),j,nn)* this%geenoptfac(k)
            gkgrad(3*j-0,m) = gkgrad(3*j-0,m) + uao(this%geenidx(k,2),i,nn) * &
                              uzao(this%geenidx(k,1),j,nn)* this%geenoptfac(k)

            gklapli(i,m) = gklapli(i,m) + u2ao(this%geenidx(k,2),i,nn) * &
                            uao(this%geenidx(k,1),j,nn)* this%geenoptfac(k)
            gklapli(j,m) = gklapli(j,m) + u2ao(this%geenidx(k,1),j,nn) * &
                            uao(this%geenidx(k,2),i,nn)* this%geenoptfac(k)
          enddo
        endif
      enddo
    enddo
    if (deriveParams) then
      do m=offset+1,offset + this%geenoptnum
        gklapl(m) = sum(gklapli(:,m))
      end do
      offset = offset+this%geenoptnum
    end if

  end subroutine anisoeen_eval

  subroutine anisoeenn_eval(this,gTerm,gDeriv,gLapli,deriveParams,offset,nn,gk,gkgrad,gklapli,gklapl)
    class(JasAnisoEENNType), intent(inout) :: this
    real(r8), intent(inout) :: gTerm, gDeriv(3*ne), gLapli(ne)
    logical, intent(in) :: deriveParams
    integer, intent(inout) :: offset
    integer, intent(in) :: nn
    real(r8), optional, intent(inout) :: gk(:),gkgrad(:,:),gklapli(:,:),gklapl(:)
    integer :: i,j,k,m

    do i=1, ne
      do j=i+1,ne
        do k=1, this%geennnum
          gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,1),i,nn) * uao(this%geennidx(k,2),j,nn)

          gDeriv(3*i-2) = gDeriv(3*i-2) + this%geenn(k) * uao(this%geennidx(k,2),j,nn) * &
                          uxao(this%geennidx(k,1),i,nn)
          gDeriv(3*i-1) = gDeriv(3*i-1) + this%geenn(k) * uao(this%geennidx(k,2),j,nn) * &
                          uyao(this%geennidx(k,1),i,nn)
          gDeriv(3*i-0) = gDeriv(3*i-0) + this%geenn(k) * uao(this%geennidx(k,2),j,nn) * &
                          uzao(this%geennidx(k,1),i,nn)

          gDeriv(3*j-2) = gDeriv(3*j-2) + this%geenn(k) * uao(this%geennidx(k,1),i,nn) * &
                          uxao(this%geennidx(k,2),j,nn)
          gDeriv(3*j-1) = gDeriv(3*j-1) + this%geenn(k) * uao(this%geennidx(k,1),i,nn) * &
                          uyao(this%geennidx(k,2),j,nn)
          gDeriv(3*j-0) = gDeriv(3*j-0) + this%geenn(k) * uao(this%geennidx(k,1),i,nn) * &
                          uzao(this%geennidx(k,2),j,nn)

          gLapli(i) = gLapli(i) + this%geenn(k) * u2ao(this%geennidx(k,1),i,nn) * uao(this%geennidx(k,2),j,nn)
          gLapli(j) = gLapli(j) + this%geenn(k) * u2ao(this%geennidx(k,2),j,nn) * uao(this%geennidx(k,1),i,nn)

          !this%geenn(k) symmetric w.r.t. basis function
          gTerm = gTerm + this%geenn(k) * uao(this%geennidx(k,2),i,nn) * uao(this%geennidx(k,1),j,nn)

          gDeriv(3*i-2) = gDeriv(3*i-2) + this%geenn(k) * uao(this%geennidx(k,1),j,nn) * &
                          uxao(this%geennidx(k,2),i,nn)
          gDeriv(3*i-1) = gDeriv(3*i-1) + this%geenn(k) * uao(this%geennidx(k,1),j,nn) * &
                          uyao(this%geennidx(k,2),i,nn)
          gDeriv(3*i-0) = gDeriv(3*i-0) + this%geenn(k) * uao(this%geennidx(k,1),j,nn) * &
                          uzao(this%geennidx(k,2),i,nn)

          gDeriv(3*j-2) = gDeriv(3*j-2) + this%geenn(k) * uao(this%geennidx(k,2),i,nn) * &
                          uxao(this%geennidx(k,1),j,nn)
          gDeriv(3*j-1) = gDeriv(3*j-1) + this%geenn(k) * uao(this%geennidx(k,2),i,nn) * &
                          uyao(this%geennidx(k,1),j,nn)
          gDeriv(3*j-0) = gDeriv(3*j-0) + this%geenn(k) * uao(this%geennidx(k,2),i,nn) * &
                          uzao(this%geennidx(k,1),j,nn)

          gLapli(i) = gLapli(i) + this%geenn(k) * u2ao(this%geennidx(k,2),i,nn) * uao(this%geennidx(k,1),j,nn)
          gLapli(j) = gLapli(j) + this%geenn(k) * u2ao(this%geennidx(k,1),j,nn) * uao(this%geennidx(k,2),i,nn)

        enddo
        if (deriveParams) then
          do k=1, this%geennnum
            m = offset + this%geennoptidx(k)

            gk(m) = gk(m) +  uao(this%geennidx(k,1),i,nn) * uao(this%geennidx(k,2),j,nn) * this%geennoptfac(k)
            gkgrad(3*i-2,m) = gkgrad(3*i-2,m) + uao(this%geennidx(k,2),j,nn) * &
                              uxao(this%geennidx(k,1),i,nn)* this%geennoptfac(k)
            gkgrad(3*i-1,m) = gkgrad(3*i-1,m) + uao(this%geennidx(k,2),j,nn) * &
                              uyao(this%geennidx(k,1),i,nn)* this%geennoptfac(k)
            gkgrad(3*i-0,m) = gkgrad(3*i-0,m) + uao(this%geennidx(k,2),j,nn) * &
                              uzao(this%geennidx(k,1),i,nn)* this%geennoptfac(k)

            gkgrad(3*j-2,m) = gkgrad(3*j-2,m) + uao(this%geennidx(k,1),i,nn) * &
                              uxao(this%geennidx(k,2),j,nn)* this%geennoptfac(k)
            gkgrad(3*j-1,m) = gkgrad(3*j-1,m) + uao(this%geennidx(k,1),i,nn) * &
                              uyao(this%geennidx(k,2),j,nn)* this%geennoptfac(k)
            gkgrad(3*j-0,m) = gkgrad(3*j-0,m) + uao(this%geennidx(k,1),i,nn) * &
                              uzao(this%geennidx(k,2),j,nn)* this%geennoptfac(k)

            gklapli(i,m) = gklapli(i,m) + u2ao(this%geennidx(k,1),i,nn) * &
                            uao(this%geennidx(k,2),j,nn)* this%geennoptfac(k)
            gklapli(j,m) = gklapli(j,m) + u2ao(this%geennidx(k,2),j,nn) * &
                            uao(this%geennidx(k,1),i,nn)* this%geennoptfac(k)
            !geenn(k) symmetric w.r.t. basis function
            gk(m) = gk(m) +  uao(this%geennidx(k,2),i,nn) * uao(this%geennidx(k,1),j,nn)* this%geennoptfac(k)
            gkgrad(3*i-2,m) = gkgrad(3*i-2,m) + uao(this%geennidx(k,1),j,nn) * &
                              uxao(this%geennidx(k,2),i,nn)* this%geennoptfac(k)
            gkgrad(3*i-1,m) = gkgrad(3*i-1,m) + uao(this%geennidx(k,1),j,nn) * &
                              uyao(this%geennidx(k,2),i,nn)* this%geennoptfac(k)
            gkgrad(3*i-0,m) = gkgrad(3*i-0,m) + uao(this%geennidx(k,1),j,nn) * &
                              uzao(this%geennidx(k,2),i,nn)* this%geennoptfac(k)

            gkgrad(3*j-2,m) = gkgrad(3*j-2,m) + uao(this%geennidx(k,2),i,nn) * &
                              uxao(this%geennidx(k,1),j,nn)* this%geennoptfac(k)
            gkgrad(3*j-1,m) = gkgrad(3*j-1,m) + uao(this%geennidx(k,2),i,nn) * &
                              uyao(this%geennidx(k,1),j,nn)* this%geennoptfac(k)
            gkgrad(3*j-0,m) = gkgrad(3*j-0,m) + uao(this%geennidx(k,2),i,nn) * &
                              uzao(this%geennidx(k,1),j,nn)* this%geennoptfac(k)

            gklapli(i,m) = gklapli(i,m) + u2ao(this%geennidx(k,2),i,nn) * &
                            uao(this%geennidx(k,1),j,nn)* this%geennoptfac(k)
            gklapli(j,m) = gklapli(j,m) + u2ao(this%geennidx(k,1),j,nn) * &
                            uao(this%geennidx(k,2),i,nn)* this%geennoptfac(k)
          enddo
        endif
      enddo
    enddo
    if (deriveParams) then
      do m=offset+1,offset + this%geennoptnum
        gklapl(m) = sum(gklapli(:,m))
      end do
       offset = offset+this%geennoptnum
    end if

  end subroutine anisoeenn_eval

   ! -------------   internal subroutines -------------
  subroutine mread_ao_idx(num,tidx,tnums,tnump,tnumd,tnumf,tJasAOLines,lines,idx,resorted,tidx_temp)
    integer, intent(out) :: tnums,tnump,tnumd,tnumf
    integer, intent(out), allocatable :: tidx(:), tidx_temp(:)
    logical,  intent(out) :: resorted
    integer :: tpos,offset,tnum, tnWords, i,j,k
    character(len=10) :: tword(20)
    type(jasAOLines), intent(inout):: tJasAOLines
    integer, intent(inout) :: num, idx
    character(len=*), intent(in) :: lines(:)! lines array
    character(len=120) :: temp_input
    character(len=3) :: temp_number

    if (allocated(tidx_temp)) deallocate(tidx_temp)
    allocate(tidx_temp(num))
    tnums = 0
    tnump = 0
    tnumd = 0
    tnumf = 0
    resorted = .false.
    write(tword(1),'(a,i3)') "idx " , num
    tJasAOLines%aoTypeEntry = tword(1)  ! save format for output
    if (allocated(tidx)) deallocate(tidx)
    allocate(tidx(num))
    tJasAOLines%n = 0
    tnum=0
    do
      call tokenize(lines(idx),tword,tnWords)
      if (tnum+tnWords > num) call abortp("$jastrow: illegal format for anisotropic terms")
      do k=1,tnWords
        read(tword(k),*) tidx_temp(tnum+k)
      end do
      tnum = tnum + tnWords
      tJasAOLines%n = tJasAOLines%n + 1
      !Change size in jasAOLines-type if needed
      if (tJasAOLines%n > SIZE(tJasAOLines%aoLines)) then
        call abortp("Too many input lines in anisotropic terms defintion.")
      endif
      tJasAOLines%aoLines(tJasAOLines%n) = lines(idx)
      if (tnum==num) exit
      idx = idx+1
    end do
    !increase counter for s, p, d or f terms
    !determine p,d or f by index
    !sort anisotropic terms w.r.t. s then p then d then f, whereas each sorted w.r.t.
    ! center number
    ! i.e. first all p terms, then all d terms, then all f terms
    tpos = 1
    do k=1,ncenter
      offset = 1
      do j=1,k-1
        offset = offset + AOIdx%nums(j) + 3*AOIdx%nump(j) + 6*AOIdx%numd(j) + &
                  10*AOIdx%numf(j) + 15*AOIdx%numg(j)
      enddo
      do i=offset,offset+AOIdx%nums(k)-1
        do j=1,size(tidx_temp)
          !if (tidx_temp(j)==i) then
          !  call abortp('mread_ao_idx: no s function supported as anisotropic function')
          !endif
          if (tidx_temp(j)==i) then
            tnums = tnums +1
            tidx(tpos) = tidx_temp(j)
            tpos = tpos +1
          endif
        enddo
      enddo
    enddo
    !p terms
    do k=1,ncenter
      offset = 1
      do j=1,k-1
        offset = offset + AOIdx%nums(j) + 3*AOIdx%nump(j) + 6*AOIdx%numd(j) + &
                    10*AOIdx%numf(j) + 15*AOIdx%numg(j)
      enddo
      offset = offset+AOIdx%nums(k)
      do i=offset,offset+3*AOIdx%nump(k)-1
        do j=1,size(tidx_temp)
          if (tidx_temp(j)==i) then
            tnump = tnump +1
            tidx(tpos) = tidx_temp(j)
            tpos = tpos +1
          endif
        enddo
      enddo
    enddo
    !d terms
    do k=1,ncenter
      offset = 1
      do j=1,k-1
        offset = offset + AOIdx%nums(j) + 3*AOIdx%nump(j) + 6*AOIdx%numd(j) + &
                  10*AOIdx%numf(j) + 15*AOIdx%numg(j)
      enddo
      offset = offset+AOIdx%nums(k)+3*AOIdx%nump(k)
      do i=offset,offset+6*AOIdx%numd(k)-1
        do j=1,size(tidx_temp)
          if (tidx_temp(j)==i) then
            tnumd = tnumd +1
            tidx(tpos) = tidx_temp(j)
            tpos = tpos +1
          endif
        enddo
      enddo
    enddo
    !f terms
    do k=1,ncenter
      offset = 1
      do j=1,k-1
        offset = offset + AOIdx%nums(j) + 3*AOIdx%nump(j) + 6*AOIdx%numd(j) + &
                  10*AOIdx%numf(j) + 15*AOIdx%numg(j)
      enddo
      offset = offset+AOIdx%nums(k)+3*AOIdx%nump(k)+6*AOIdx%numd(k)
      do i=offset,offset+10*AOIdx%numf(k)-1
        do j=1,size(tidx_temp)
          if (tidx_temp(j)==i) then
            tnumf = tnumf +1
            tidx(tpos) = tidx_temp(j)
            tpos = tpos +1
          endif
        enddo
      enddo
    enddo
    if (MASTER) then
      if (.not.( all(tidx==tidx_temp) )) then
        resorted = .true.
        do k=1 , (SIZE(tidx) / 10) +1
          tJasAOLines%n = k
          temp_input = ''
          if (tJasAOLines%n > SIZE(tJasAOLines%aoLines)) then
            call abortp("Too many input lines after rewriting anisotropic terms defintion.")
          endif
          do i=1, min(k*10,SIZE(tidx)-((k-1)*10))
            write(temp_number,'(i3)') tidx((k-1)*10+i)
            temp_input = trim(temp_input) // trim(temp_number)
          enddo
          tJasAOLines%aoLines(tJasAOLines%n) = temp_input
        enddo
      endif
    endif
    call assert((num==tpos-1),&
      'read_ao terms: Missmatch while reading ao idx')

    !Check for consistency
    call assert((num==tnums+tnump+tnumd+tnumf),&
      'read_ao terms: Inconsistent number of anistropic parameters')
  end subroutine mread_ao_idx

  subroutine mread_ao_nuc(num,tidx,tnums,tnump,tnumd,tnumf,tJasAOLines,nlin,lines,idx)
    integer, intent(out) :: num,tnums,tnump,tnumd,tnumf
    integer, intent(out), allocatable :: tidx(:)
    integer :: tnWords, i,j,k,l,m
    character(len=10) :: tword(20)
    type(jasAOLines), intent(inout):: tJasAOLines
    character(len=*), intent(in) :: lines(:)! lines array
    integer, intent(inout) :: idx
    integer, intent(in) :: nlin

    write(tword(1),'(a,i3)') "nuc ", nlin
    tJasAOLines%aoTypeEntry = tword(1)   ! save format for output
    do i=1,nlin
      tJasAOLines%aoLines(i) = lines(idx+i)
    end do
    tJasAOLines%n = nLin
    num = 0
    tnums = 0
    tnump = 0
    tnumd = 0
    tnumf = 0
    do i=1,nlin
      call tokenize(lines(idx+i),tword,tnWords)
      if (tword(1)=="s") then
        do j=2,tnWords
          read(tword(j),*) k
          num = num + AOIdx%nums(k)
          !increase counter for s type ao term
          tnums = tnums + AOIdx%nums(k)
        end do
     else if (tword(1)=="p") then
        do j=2,tnWords
          read(tword(j),*) k
          num = num + 3*AOIdx%nump(k)
          !increase counter for p type ao term
          tnump = tnump + 3*AOIdx%nump(k)
        end do
      else if (tword(1)=="d") then
        do j=2,tnWords
          read(tword(j),*) k
          num = num + 6*AOIdx%numd(k)
          !increase counter for d type ao term
          tnumd = tnumd + 6*AOIdx%numd(k)
        end do
      else if (tword(1)=="f") then
        do j=2,tnWords
          read(tword(j),*) k
          num = num + 10*AOIdx%numf(k)
          !increase counter for f type ao term
          tnumf = tnumf + 10*AOIdx%numf(k)
        end do
      else
        call abortp("Unknown AO type in mread_ao_nuc")
      end if
    end do
    if (allocated(tidx)) deallocate(tidx)
    allocate(tidx(num))
    m = 0
    do i=1,nlin
      call tokenize(lines(idx+i),tword,tnWords)
      if (tword(1)=="s") then
        do j=2,tnWords
          read(tword(j),*) k
          do l=1,AOIdx%nums(k)
            tidx(m+1) = AOIdx%s(k,l)
            m = m+1
          end do
        end do
      else if (tword(1)=="p") then
        do j=2,tnWords
          read(tword(j),*) k
          do l=1,AOIdx%nump(k)
            tidx(m+1) = AOIdx%p(k,l)
            tidx(m+2) = AOIdx%p(k,l) + 1
            tidx(m+3) = AOIdx%p(k,l) + 2
            m = m+3
          end do
        end do
      else if (tword(1)=="d") then
        do j=2,tnWords
          read(tword(j),*) k
          do l=1,AOIdx%numd(k)
            tidx(m+1) = AOIdx%d(k,l)
            tidx(m+2) = AOIdx%d(k,l) + 1
            tidx(m+3) = AOIdx%d(k,l) + 2
            tidx(m+4) = AOIdx%d(k,l) + 3
            tidx(m+5) = AOIdx%d(k,l) + 4
            tidx(m+6) = AOIdx%d(k,l) + 5
            m = m+6
          end do
        end do
      else if (tword(1)=="f") then
        do j=2,tnWords
          read(tword(j),*) k
          do l=1,AOIdx%numf(k)
            tidx(m+1) = AOIdx%f(k,l)
            tidx(m+2) = AOIdx%f(k,l) + 1
            tidx(m+3) = AOIdx%f(k,l) + 2
            tidx(m+4) = AOIdx%f(k,l) + 3
            tidx(m+5) = AOIdx%f(k,l) + 4
            tidx(m+6) = AOIdx%f(k,l) + 5
            tidx(m+7) = AOIdx%f(k,l) + 6
            tidx(m+8) = AOIdx%f(k,l) + 7
            tidx(m+9) = AOIdx%f(k,l) + 8
            tidx(m+10)= AOIdx%f(k,l) + 9
            m = m+10
          end do
        end do
      end if
    end do
    idx = idx + nlin
    !sort anisotropic terms w.r.t. p then d then f, whereas each sorted w.r.t.
    ! center number
    ! i.e. first all p terms, then all d terms, then all f terms
    !Check for consistency
    call assert((num==tnums+tnump+tnumd+tnumf),'read_ao terms: Inconsistent number of anistropic parameters')
  end subroutine mread_ao_nuc

  subroutine mread_ao_all(num,tidx,tnums,tnump,tnumd,tnumf,tJasAOLines,tpos,tword,tnWords)
    integer, intent(out) :: num,tnums,tnump,tnumd,tnumf
    integer, intent(out), allocatable :: tidx(:)
    integer :: tend, i,k,l,m
    character(len=*), intent(in) :: tword(20)
    type(jasAOLines), intent(inout):: tJasAOLines
    integer, intent(inout) :: tpos
    integer, intent(in) :: tnWords

    ! the entries after the 5th=="all" are: "p"|"d"|"f"
    tJasAOLines%n = 0
    tJasAOLines%aoTypeEntry = "all"
    tend = 0
    do i=tpos,tnWords
      if ( .not. ( tword(i)=="s" .or. tword(i)=="p" .or. tword(i)=="d" .or. tword(i)=="f")) exit
      tJasAOLines%aoTypeEntry = trim(tJasAOLines%aoTypeEntry) // " " // trim(tword(i))
      tend = tend +1
    end do
    num = 0
    tnums = 0
    tnump = 0
    tnumd = 0
    tnumf = 0
    do i=tpos,tpos+tend-1
      if (tword(i)=="s") then
        do k=1,ncenter
          num = num + AOIdx%nums(k)
          !increase counter for s type ao term
          tnums = tnums + AOIdx%nums(k)
        end do
      else if (tword(i)=="p") then
        do k=1,ncenter
          num = num + 3*AOIdx%nump(k)
          !increase counter for p type ao term
          tnump = tnump + 3*AOIdx%nump(k)
        end do
      else if (tword(i)=="d") then
        do k=1,ncenter
          num = num + 6*AOIdx%numd(k)
          !increase counter for d type ao term
          tnumd = tnumd + 6*AOIdx%numd(k)
        end do
      else if (tword(i)=="f") then
        do k=1,ncenter
          num = num + 10*AOIdx%numf(k)
          !increase counter for f type ao term
          tnumf = tnumf + 10*AOIdx%numf(k)
        end do
      else
        call abortp("Unknown AO type in mread_ao_all")
      end if
    end do
    if (allocated(tidx)) deallocate(tidx)
    allocate(tidx(num))
    m = 0
    do i=tpos,tpos+tend-1
      if (tword(i)=="s") then
        do k=1,ncenter
          do l=1,AOIdx%nums(k)
            tidx(m+1) = AOIdx%s(k,l)
            m = m+1
          end do
        end do
      else if (tword(i)=="p") then
        do k=1,ncenter
          do l=1,AOIdx%nump(k)
            tidx(m+1) = AOIdx%p(k,l)
            tidx(m+2) = AOIdx%p(k,l) + 1
            tidx(m+3) = AOIdx%p(k,l) + 2
            m = m+3
          end do
        end do
      else if (tword(i)=="d") then
        do k=1,ncenter
          do l=1,AOIdx%numd(k)
            tidx(m+1) = AOIdx%d(k,l)
            tidx(m+2) = AOIdx%d(k,l) + 1
            tidx(m+3) = AOIdx%d(k,l) + 2
            tidx(m+4) = AOIdx%d(k,l) + 3
            tidx(m+5) = AOIdx%d(k,l) + 4
            tidx(m+6) = AOIdx%d(k,l) + 5
            m = m+6
          end do
        end do
      else if (tword(i)=="f") then
        do k=1,ncenter
          do l=1,AOIdx%numf(k)
            tidx(m+1) = AOIdx%f(k,l)
            tidx(m+2) = AOIdx%f(k,l) + 1
            tidx(m+3) = AOIdx%f(k,l) + 2
            tidx(m+4) = AOIdx%f(k,l) + 3
            tidx(m+5) = AOIdx%f(k,l) + 4
            tidx(m+6) = AOIdx%f(k,l) + 5
            tidx(m+7) = AOIdx%f(k,l) + 6
            tidx(m+8) = AOIdx%f(k,l) + 7
            tidx(m+9) = AOIdx%f(k,l) + 8
            tidx(m+10)= AOIdx%f(k,l) + 9
            m = m+10
          end do
        end do
      end if
    end do
    !Check for consistency
    call assert((num==tnums+tnump+tnumd+tnumf),'read_ao terms: Inconsistent number of anisotropic parameters')
    tpos = tpos + tend
  end subroutine mread_ao_all

  ! ------------- get Count -------------
  subroutine anisoen_getcount(this,count)
    class(JasAnisoENType), intent(in) :: this
    integer,intent(out) :: count
    count = this%goptnum
  end subroutine anisoen_getcount

  subroutine anisoeen_getcount(this,count)
    class(JasAnisoEENType), intent(in) :: this
    integer,intent(out) :: count
    count = this%geenoptnum
  end subroutine anisoeen_getcount

  subroutine anisoeenn_getcount(this,count)
    class(JasAnisoEENNType), intent(in) :: this
    integer,intent(out) :: count
    count = this%geennoptnum
  end subroutine anisoeenn_getcount

  subroutine anisoall_getcount(this,optMode,np_one,np_two)
    class(JasAnisoType), intent(in) :: this
    integer, intent(in) :: optMode
    integer, intent(out) :: np_one,np_two
    integer :: temp

    np_one = 0
    np_two = 0
    select case(optMode)
    case (OPT_LIN, OPT_NUM_LIN,OPT_ALL, OPT_NUM_ALL)
      call this%en%getParamCount(temp)
      np_one = temp
      call this%een%getParamCount(temp)
      np_two = np_two + temp
      call this%eenn%getParamCount(temp)
      np_two = np_two + temp
    case (OPT_NONLIN)
      np_one = 0
      np_two = 0
    case default
      call abortp("anisoall_getcount: optMode not implemented")
   end select

  end subroutine anisoall_getcount

  function anisoall_getNumberofParams(this) result(nParams)
    class(JasAnisoType), intent(in) :: this
    integer :: nParams
    integer :: temp, tnParams
    tnParams = 0
    temp = 0
    call this%en%getParamCount(temp)
    tnParams = temp
    call this%een%getParamCount(temp)
    tnParams = tnParams + temp
    call this%eenn%getParamCount(temp)
    tnParams = tnParams + temp
    nParams = tnParams
  end function anisoall_getNumberofParams

  function anisoall_getNumberofParams_J1(this) result(nParams)
    class(JasAnisoType), intent(in) :: this
    integer :: nParams
    integer :: temp, tnParams
    tnParams = 0
    temp = 0
    call this%en%getParamCount(temp)
    tnParams = temp
    nParams = tnParams
  end function anisoall_getNumberofParams_J1

  function anisoall_getNumberofParams_J2(this) result(nParams)
    class(JasAnisoType), intent(in) :: this
    integer :: nParams
    integer :: temp, tnParams
    tnParams = 0
    temp = 0
    call this%een%getParamCount(temp)
    tnParams = tnParams + temp
    call this%eenn%getParamCount(temp)
    tnParams = tnParams + temp
    nParams = tnParams
  end function anisoall_getNumberofParams_J2


  ! ------------- various output routines -------------
  subroutine anisoall_output(this)
    class(JasAnisoType), intent(in) :: this
    integer i
    write(iul, "(/A,I4)") "No. of en AO Jastrow terms = ", this%en%gnum
    do i=1, this%en%gnum
      write(iul,*) this%en%gl(i)
    end do
    write(iul, "(/A,I4)") "No. of een AO Jastrow terms = ", this%een%geennum
    do i=1,this%een%geennum
      write(iul,*) this%een%geen(i)
    end do
    write(iul, "(/A,I4)") "No. of eenn AO Jastrow terms = ", this%eenn%geennnum
    do i=1,this%eenn%geennnum
      write(iul,*) this%eenn%geenn(i)
    end do
  end subroutine anisoall_output


  subroutine anisoall_shortoutput(this,iu)
    class(JasAnisoType), intent(in) :: this
    integer, intent(in) :: iu

    write(iu,'(i4,a)') this%en%gnum,' anisotropic en (AO) terms'
    write(iu,'(i4,a)') this%een%geennum,' anisotropic een (AO) terms'
    write(iu,'(i4,a)') this%eenn%geennnum,' anisotropic eenn (AO) terms'
    if (this%en%gnum /=  this%en%goptnum ) then
      write(iu,'(a,i4,a)') 'using pairwise symmetry of ',size(this%en%symmAOen),' functions for anisotropic en terms'
      write(iu,'(a,i4,a)') '  leading to ',this%en%goptnum,' independent anisotropic en parameters'
    endif
    if (this%een%geennum /=  this%een%geenoptnum ) then
      write(iu,'(a,i4,a)') 'using pairwise symmetry of ',size(this%een%symmAOeen),' functions for anisotropic een terms'
      write(iu,'(a,i4,a)') '  leading to ',this%een%geenoptnum,' independent anisotropic een parameters'
    endif
    if (this%eenn%geennnum /= this%eenn%geennoptnum ) then
      write(iu,'(a,i4,a)') 'using pairwise symmetry of ',size(this%eenn%symmAOeenn),' functions for anisotropic eenn terms'
      write(iu,'(a,i4,a)') '  leading to ',this%eenn%geennoptnum,' independent anisotropic eenn parameters'
    endif

  end subroutine anisoall_shortoutput

  subroutine anisoall_output_param_new(this,iu)
    class(JasAnisoType), intent(in) :: this
    integer, intent(in) :: iu
    integer :: i

    if (this%en%gnum > 0) then
      do i=1,this%en%gnum
        write(iu, "(ES15.7)") this%en%gl(i)
      end do
    end if
    if (this%een%geennum > 0 ) then
      do i=1,this%een%geennum
        write(iu, "(ES15.7)") this%een%geen(i)
      end do
    end if
    if (this%eenn%geennnum > 0) then
      do i=1,this%eenn%geennnum
        write(iu, "(ES15.7)") this%eenn%geenn(i)
      end do
    end if

  end subroutine anisoall_output_param_new

  subroutine anisoall_output_header_new(this,tline)
    class(JasAnisoType), intent(in) :: this
    character(len=120),intent(inout) :: tline

    if (this%en%gnum>0) then
      write(tline, "(3a)") trim(tline), " ao ",trim(this%en%curJasAOLines%aoTypeEntry)
    endif
    if (this%een%geennum>0) then
        write(tline, "(3a)") trim(tline)," eenao ",trim(this%een%currJaseenAOLines%aoTypeEntry)
    end if
    if (this%eenn%geennnum>0) then
      write(tline, "(3a)") trim(tline)," eennao ",trim(this%eenn%currJaseennAOLines%aoTypeEntry)
    end if

  end subroutine anisoall_output_header_new

  subroutine anisoall_output_header_new_lines(this,iu)
    class(JasAnisoType), intent(in) :: this
    integer, intent(in) :: iu
    integer :: i

    if (this%en%gnum>0) then
      do i=1,this%en%curJasAOLines%n
        write(iu, "(a)") this%en%curJasAOLines%aoLines(i)
      end do
    end if
    if (this%een%geennum>0) then
      do i=1,this%een%currJaseenAOLines%n
        write(iu, "(a)") this%een%currJaseenAOLines%aoLines(i)
      end do
    end if
    if (this%eenn%geennnum>0) then
      do i=1,this%eenn%currJaseennAOLines%n
        write(iu, "(a)") this%eenn%currJaseennAOLines%aoLines(i)
      end do
    end if
  end subroutine anisoall_output_header_new_lines

  ! ------------- get Params -------------
  subroutine anisoall_getParamVector(this,optMode,p,length)
    class(JasAnisoType), intent(in) :: this
    integer, intent(in) :: optMode
    real(r8), allocatable, intent(out) :: p(:)
    integer, intent(out) :: length
    integer :: start, aostart, t1, t2

    call this%getParamCount(optMode,t1,t2)
    length = t1+t2
    if (allocated(p)) deallocate(p)
    allocate(p(length))
    start = 1
    aostart = 1
    if (this%m_useAOJasTerms) then
      !enAO -  all functions
      do aostart=0, (this%en%goptnum-1)
        p(start+aostart) = this%en%gl(this%en%goptidx_reverse(aostart+1)%idx(1)) * &
                            this%en%goptidx_reverse(aostart+1)%factor(1)
        !if (MASTER) write(*,*) this%en%gl(this%en%goptidx_reverse(aostart+1)%idx(1)) * this%en%goptidx_reverse(aostart+1)%factor(1), this%en%gl(aostart+1)
      enddo
      !if (MASTER) write(*,*) 'en: anisoall_getParamVector done'
      start = start+this%en%goptnum
      !eenAO
      do aostart=0, (this%een%geenoptnum-1)
        p(start+aostart) = this%een%geen(this%een%geenoptidx_reverse(aostart+1)%idx(1)) * &
                            this%een%geenoptidx_reverse(aostart+1)%factor(1)
        !if (MASTER) write(*,*) this%een%geen(this%een%geenoptidx_reverse(aostart+1)%idx(1)) * this%een%geenoptidx_reverse(aostart+1)%factor(1), this%een%geen(aostart+1)
      enddo
      !if (MASTER) write(*,*) 'een: anisoall_getParamVector done'
      start = start+this%een%geenoptnum
      !eennAO
      do aostart=0, (this%eenn%geennoptnum-1)
        p(start+aostart) = this%eenn%geenn(this%eenn%geennoptidx_reverse(aostart+1)%idx(1)) * &
                            this%eenn%geennoptidx_reverse(aostart+1)%factor(1)
        !if (MASTER) write(*,*) this%eenn%geenn(this%eenn%geennoptidx_reverse(aostart+1)%idx(1)) * this%eenn%geennoptidx_reverse(aostart+1)%factor(1), this%eenn%geenn(aostart+1)
      enddo
      !if (MASTER) write(*,*) 'eenn: anisoall_getParamVector done'
      start = start+this%eenn%geennoptnum
    endif
  end subroutine anisoall_getParamVector

  ! ------------- set Params -------------
  subroutine anisoall_setParamVector(this,optMode,p)
    class(JasAnisoType), intent(inout) :: this
    integer, intent(in) :: optMode
    real(r8), intent(in) :: p(:)
    integer :: start,aostart
    !nedded for symmetry adaption
    integer             ::  i,j

    call this%getParamCount(optmode,i,j) !i,j used to add npJ1 and npJ2
    call assert(size(p)==(i+j),"size missmatch in anisoall_setParamVector")
    start = 1
    aostart = 1
    if (this%m_useAOJasTerms) then
      !enAO - all functions
      do aostart=0, (this%en%goptnum-1)
        !if (MASTER) write(*,*) aostart+1,'th parameter used for ',this%en%goptidx_reverse(aostart+1)%idx,' value: ',p(start+aostart)
        do i=1, SIZE(this%en%goptidx_reverse(aostart+1)%idx)
          this%en%gl(this%en%goptidx_reverse(aostart+1)%idx(i)) =  p(start+aostart) * &
                      this%en%goptidx_reverse(aostart+1)%factor(i)
        enddo
      enddo
      start = start + this%en%goptnum
      !een AOTerms
      do aostart=0, (this%een%geenoptnum-1)
        !if (MASTER) write(*,*) aostart+1,'th parameter used for ',this%een%geenoptidx_reverse(aostart+1)%idx,' value: ',p(start+aostart)
        do i=1, SIZE(this%een%geenoptidx_reverse(aostart+1)%idx)
          this%een%geen(this%een%geenoptidx_reverse(aostart+1)%idx(i)) =  p(start+aostart) * &
                          this%een%geenoptidx_reverse(aostart+1)%factor(i)
        enddo
      enddo
      start = start+this%een%geenoptnum
      !eenn AOTerms
      do aostart=0, (this%eenn%geennoptnum-1)
        !if (MASTER) write(*,*) aostart+1,'th parameter used for ',this%eenn%geennoptidx_reverse(aostart+1)%idx,' value: ',p(start+aostart)
        do i=1, SIZE(this%eenn%geennoptidx_reverse(aostart+1)%idx)
          this%eenn%geenn(this%eenn%geennoptidx_reverse(aostart+1)%idx(i)) =  p(start+aostart) * &
                      this%eenn%geennoptidx_reverse(aostart+1)%factor(i)
        enddo
      enddo
      start = start+this%eenn%geennoptnum
    endif

  end subroutine anisoall_setParamVector

end module jastrowAniso_m
