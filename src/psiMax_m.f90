! Copyright (C) 2012-2015, 2017-2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

! module to perform the maximization of psi**2
! implemented as minimization of -ln(psi**2)
!
! new implementation with abstract minimizer classes
!
! parallel version:
!   both old and new codes run optimizer locally and collect local maxima on MASTER (in add_to_list/addTOList)
!   ONLY the MASTER keeps/stores/sorts references
!   new code: ONLY MASTER has/needs referenceContainer object allocated (mRC_p)

module psiMax_m

   use kinds_m, only: r8
   use error_m
   use parsing_m
   use global_m, only: getNElec, iul, iull, MASTER, logmode, getMyTaskId
   use eloc_m, only: eloc
   use eConfigs_m, only: eConfigArray, eConfigArray_new, eConfigArray_set
   use elocData_m, only: eloc_getStopAtSingularity, eloc_setStopAtSingularity, elEkini
   use wfData_m, only: atoms, getNNuc, getNAlpha
   use atom_m, only: atoms_getPositionMatrix
   use fPsi2_m
   use minimizer_w_sing_module
   use minimizer_ws_factory_module, only: create_ws_minimizer
   use randomWalker_m, only: RandomWalker, pos, resetTo

   use maxRawData_m, only: maxraw_isInitialized, maxraw_destroy, maxraw_add, maxraw_writeResults

   use maxAnalysis_m, only: maxana_isInitialized, maxana_destroy, maxana_add, maxana_writeResults, &
       maxana_getFirstF, maxana_getDiffMax
   use maxBasins_m, only: maxbas_isInitialized, maxbas_destroy, maxbas_add, maxbas_writeResults
   use mpiInterface_m, only: myMPIReduceSumDouble

   implicit none
   private
   public :: psimax, RAWDATA, NOANALYSIS
   !!!public :: psimax_init, psimax_destroy, psimax_opt, psimax_writeResults, psimax_writeParams, psimax_doMaxAnalysis, &
   !!!          psimax_doMaxSearch, psimax_getDiffMax, psimax_getFirstF, rebuildvec, writeCoordsToVec

   integer, parameter :: NOANALYSIS = 0, MAXIMA = 1, BASINS = 2, RAWDATA = 3

   type psimax
      private
      class(minimizer_w_sing), pointer :: minimizer_p => null()
      type(fctn_psi2)                  :: fg
      integer                          :: analyze_mode_ = NOANALYSIS
      logical                          :: allow_singularities_ = .true.
      integer                          :: verbose_ = 0
      logical                          :: save_opt_paths_ = .false.
      integer                          :: path_count_ = 0
   contains
      procedure :: init => psimax_init
      procedure :: destroy => psimax_destroy
      procedure :: do_max_search => psimax_doMaxSearch
      procedure :: do_max_analysis => psimax_doMaxAnalysis
      procedure :: get_diff_max => psimax_getDiffMax
      procedure :: get_first_f => psimax_getFirstF
      procedure :: write_params => psimax_writeParams
      procedure :: opt => psimax_opt
      procedure :: is_converged => psimax_is_converged
      procedure :: function_evaluations => psimax_function_evaluations
      procedure :: iterations => psimax_iterations
      procedure :: write_results => psimax_writeResults
      procedure :: get_analyze_mode => psimax_getAnalyzeMode
   end type psimax

   ! counters variables:
   real(r8)       :: sCounter(0:9) = 0.d0  ! counters: exit codes 0:7, 8 all, 9 successful minimizer calls

contains

   subroutine psimax_init(this, lines, nl)
      class(psimax), intent(inout)    :: this
      character(len=120), intent(in) :: lines(:)
      integer, intent(in)            :: nl
      logical l1, l2, l3
      integer iflag, i, iu

      this%fg = fctn_psi2()

      this%minimizer_p => create_ws_minimizer(lines)

      call getinta(lines, nl, "verbose=", this%verbose_, iflag)
      call this%minimizer_p%set_verbose(this%verbose_)
      call this%minimizer_p%set_verbose_unit(iull)

      call this%minimizer_p%set_singularities(atoms_getPositionMatrix(atoms), &
                                              last_a=getNAlpha(), &
                                              scalings=REAL(atoms%Get_atomic_number(), r8))

      l1 = maxana_isInitialized()
      l2 = maxbas_isInitialized()
      l3 = maxraw_isInitialized()


      if ( (l1.and.(l2.or.l3)) .or. (l2.and.(l1.or.l3)) .or. (l3.and.(l1.or.l2)) ) then
         call error(' $init_max_search cannot do multiple analyses at once')
      endif

      if (l1) then
         this%analyze_mode_ = MAXIMA
      else if (l2) then
         this%analyze_mode_ = BASINS
      else if (l3) then
         this%analyze_mode_ = RAWDATA
      else
         this%analyze_mode_ = NOANALYSIS
      endif
      if (finda(lines, nl, "save_opt_paths")) then
         this%save_opt_paths_ = .true.
         this%path_count_ = 0
         iu = 300 + getMyTaskId()
         write(iu,'(i5,a)') getNNuc(), " nuclei:"
         do i = 1, getNNuc()
            write(iu,'(i4,1x,a2,1x,3f14.7)') i, atoms(i)%elem, atoms(i)%cx, atoms(i)%cy, atoms(i)%cz
         enddo
         write(iu,'(i5,a,2i5,a,i5,a,i5)') 0, " MINPATH ", 0, 0,&
         " nAlpha ", getNAlpha()," nBeta ", getNElec()-getNAlpha()
      end if

      if (finda(lines, nl, "no_allow_sing")) then
         this%allow_singularities_ = .false.
      end if

      if (logmode >= 2) then
         call this%write_params(iul)
      end if
   end subroutine psimax_init


   subroutine psimax_destroy(this)
      class(psimax), intent(in)    :: this
      !sCounter = 0.d0
      select case(this%analyze_mode_)
      case (RAWDATA)
         call maxraw_destroy()
      case (MAXIMA)
         call maxana_destroy()
      case (BASINS)
         call maxbas_destroy()
      end select
   end subroutine psimax_destroy

   function psimax_getAnalyzeMode(this) result(mode)
      class(psimax), intent(in) :: this
      integer :: mode
      mode = this%analyze_mode_
   end function psimax_getAnalyzeMode

   function psimax_doMaxSearch(this) result(res)
      class(psimax), intent(in)    :: this
      logical :: res
      res = (this%analyze_mode_ == RAWDATA .or. this%analyze_mode_ == MAXIMA .or. this%analyze_mode_ == BASINS)
   end function psimax_doMaxSearch


   function psimax_doMaxAnalysis(this) result(res)
      class(psimax), intent(in)    :: this
      logical :: res
      res = (this%analyze_mode_ == MAXIMA)
   end function psimax_doMaxAnalysis


   function psimax_getFirstF(this) result(res)
      class(psimax), intent(in)    :: this
      real(r8) :: res
      res = 0.d0
      if (this%analyze_mode_ == MAXIMA) then
         res = maxana_getFirstF()
      endif
   end function psimax_getFirstF


   function psimax_getDiffMax(this) result(res)
      class(psimax), intent(in)    :: this
      integer :: res
      res = 0
      if (this%analyze_mode_ == MAXIMA) then
         res = maxana_getDiffMax()
      endif
   end function psimax_getDiffMax


   subroutine psimax_writeParams(this, iu)
      class(psimax), intent(in) :: this
      integer, intent(in) :: iu

      call this%minimizer_p%write_params(iu)
   end subroutine psimax_writeParams


   subroutine psimax_opt(this, rw, update_walker)
      class(psimax), intent(inout)      :: this
      type(RandomWalker), intent(inout) :: rw
      logical, intent(in)               :: update_walker
      real(r8)  :: x(getNElec()), y(getNElec()), z(getNElec())
      real(r8)  :: xx(3*getNElec())
      real(r8)  :: f
      real(r8) :: maximum(3*getNElec())
      real(r8) :: sample(3*getNElec())
      real(r8) :: kinetic_energies(getNElec())
      integer :: i, n, w, iflag
      integer :: idx(getNElec())
      logical :: lold
      type(eConfigArray) :: sample_ec

      n = getNElec()

      if (this%allow_singularities_) then
         lold = eloc_getStopAtSingularity()
         call eloc_setStopAtSingularity(.false.)
      end if

      call pos(rw, x, y, z)

      if (this%verbose_ >= 2) then
         write(iull,*) 'before minimize:'
         do i = 1, n
            write(iull,'(i5,3f15.6)') i, x(i), y(i), z(i)
         end do
      end if

      do i = 1, n
         xx(3*i-2) = x(i)
         xx(3*i-1) = y(i)
         xx(3*i)   = z(i)
      end do
      call eConfigArray_new(sample_ec,n,1)
      call eConfigArray_set(sample_ec,1,x,y,z)
      sample = xx


      call this%minimizer_p%reset()
      if (this%save_opt_paths_) then
         this%path_count_ = this%path_count_ + 1
         call this%minimizer_p%write_opt_path(this%path_count_)
      end if
      call this%minimizer_p%minimize(this%fg, xx)

      f = this%minimizer_p%value()

      do i = 1, n
         x(i) = xx(3*i-2)
         y(i) = xx(3*i-1)
         z(i) = xx(3*i)
      end do
      maximum = xx

      if (this%verbose_ >= 2) then
         write(iull,*) 'after minimize:'
         do i = 1, n
            write(iull,'(i5,3f15.6)') i, x(i), y(i), z(i)
         end do
      end if


      sCounter(8) = sCounter(8) + 1.d0
      !!sCounter(iflag) = sCounter(iflag) + 1.d0
      if (this%minimizer_p%is_converged()) then
         sCounter(9) = sCounter(9) + 1.d0
         sCounter(7) = sCounter(7) + this%minimizer_p%iterations()
         sCounter(6) = sCounter(6) + this%minimizer_p%function_evaluations()
         iflag = 0
      else
         iflag = 1
      end if

      idx = (/ (i, i = 1, n) /)

      ! note: rw must contain input positions for maxbas_add
      select case (this%analyze_mode_)
      case (RAWDATA)
         if (this%minimizer_p%is_converged()) then
            call eloc(0,sample_ec,'none')
            kinetic_energies=elEkini(:,1)
            call maxraw_add(sample,kinetic_energies,maximum,f)
         end if
      case (MAXIMA)
         call maxana_add(x, y, z, f, iflag)
      case (BASINS)
         call maxbas_add(x, y, z, f, rw, idx, iflag)
      end select

      if (this%allow_singularities_) then
         call eloc_setStopAtSingularity(lold)
      end if

      if (update_walker) then
         call resetTo(rw, x, y, z)   ! this recalculates eloc, save but unnecessary
         !!! call resetTo_without_calc(rw, x, y, z)
      end if
   end subroutine psimax_opt

   function psimax_function_evaluations(this) result(res)
      class(psimax), intent(in)    :: this
      integer                      :: res
      res = this%minimizer_p%function_evaluations()
   end function psimax_function_evaluations

   function psimax_iterations(this) result(res)
      class(psimax), intent(in)    :: this
      integer                      :: res
      res = this%minimizer_p%iterations()
   end function psimax_iterations

   function psimax_is_converged(this) result(res)
      class(psimax), intent(in)    :: this
      logical                      :: res
      res = this%minimizer_p%is_converged()
   end function psimax_is_converged

   subroutine psimax_writeResults(this)
      class(psimax), intent(in)    :: this
      integer, parameter :: iu=12
      real(r8)          :: min_save, rcvCount(0:9)
      integer         :: i,k
      real(r8)          :: F

      call myMPIReduceSumDouble(sCounter, rcvCount, 10)

      if (MASTER) then

         write(iul,*)
         write(iul,*) "Summary for maxima search:"
         write(iul,*)


         write(iul,'(a41,i8)') " # minimizer calls                      :", nint(rcvCount(8))
         write(iul,'(a41,i8)') " # maxima analyzed (converged)          :", nint(rcvCount(9))
         write(iul,'(a41,i8)') " average # iterations                   :", nint(rcvCount(7) / rcvCount(9))
         write(iul,'(a41,i8)') " average # function/gradient evaluations:", nint(rcvCount(6) / rcvCount(9))
        ! write(iul,'(a)') "  exit code percentages of minimizer:"
         write(iul,*)
      endif

      select case (this%analyze_mode_)
      case (RAWDATA)
         call maxraw_writeResults()
      case (MAXIMA)
         call maxana_writeResults()
      case (BASINS)
         call maxbas_writeResults()
      end select

   end subroutine psimax_writeResults

end module psiMax_m
