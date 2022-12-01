! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module singularityCorrection_m

   ! this data structure is suitable only for collections of 3d particles (dim=3n)
   ! with non-smooth minima (singularities) at certain particle positions 
   ! (array of 3d positions)
   ! By non-smooth minima we mean discontinuous derivatives in 3 of 3n gradient
   ! components.
   ! The type allows calculation of distances and vectors to the nearest 3d (singular)
   ! position as well as setting particles at a singular position.
   ! The position data structure is assumed to be r = (x_1, y_1, z_1, x_2, ... z_n)

   use kinds_m, only: r8
   use error_m, only: assert, asserts
   use singularityParticles_m, only: singularity_particles
   use verbosity_m, only: Verbosity_t
   implicit none

   private
   public :: singularity_correction

   integer, parameter, public :: NONE = 1, CUTSTEP = 2, UMRIGAR = 3
   real(r8), parameter, private :: pi = ACOS(-1._r8)

   type, extends(Verbosity_t) :: singularity_correction
      private
      real(r8), allocatable   :: singularities_(:,:)
      real(r8), allocatable   :: scalings_(:)
      real(r8)                :: sing_thresh_ = 0._r8 ! set particle at singularity if distance smaller
      real(r8)                :: corr_thresh_ = 0._r8 ! do Umrigar/cut step correction if distance smaller
      integer                     :: last_a_ = 1         ! last particle of kind a
      integer                     :: corr_mode_ = NONE   ! correction mode
      logical                     :: scaling_ = .false.
   contains
      procedure                                :: set_singularities
      procedure                                :: set_scalings
      procedure                                :: singularities
      procedure                                :: n_singularities
      procedure                                :: set_last_a
      procedure                                :: last_a
      procedure                                :: set_correction_parameters
      procedure                                :: correction_parameters
      procedure                                :: get_nearest_singularity
      procedure                                :: update_singularity_list
      procedure                                :: set_to_singularity
      procedure                                :: correct_for_singularities
      procedure                                :: write_params_singularity_correction
      procedure                                :: write_params => write_params_singularity_correction
   end type singularity_correction

   interface singularity_correction
      module procedure constructor1
   end interface singularity_correction

contains
   function constructor1(sing_thresh, corr_thresh, corr_mode, scaling)
      real(r8), intent(in)               :: sing_thresh, corr_thresh
      integer, intent(in)                    :: corr_mode
      logical, intent(in), optional          :: scaling
      type(singularity_correction), pointer  :: constructor1
      allocate(constructor1)
      constructor1%sing_thresh_ = sing_thresh
      constructor1%corr_thresh_ = corr_thresh
      constructor1%corr_mode_   = corr_mode

      if (PRESENT(scaling) .and. scaling) constructor1%scaling_ = scaling
   end function constructor1

   subroutine set_singularities(this, s)
      class(singularity_correction), intent(inout) :: this
      real(r8), intent(in)    :: s(:,:)
      if (allocated(this%singularities_)) deallocate (this%singularities_)
      allocate (this%singularities_, source=s)
   end subroutine set_singularities

   subroutine set_scalings(this, s)
      class(singularity_correction), intent(inout) :: this
      real(r8), intent(in)    :: s(:)
      if (allocated(this%scalings_)) deallocate (this%scalings_)
      allocate (this%scalings_, source=s)
   end subroutine set_scalings

   function singularities(this) result (s)
      class(singularity_correction), intent(in)  :: this
      real(r8), allocatable                  :: s(:,:)
      allocate (s, source=this%singularities_)
   end function singularities

   pure function n_singularities(this) result (n)
      class(singularity_correction), intent(in)  :: this
      integer                                    :: n
      n = size(this%singularities_, dim=2)
   end function n_singularities

   subroutine set_last_a(this, last_a)
      class(singularity_correction), intent(inout)  :: this
      integer, intent(in)                           :: last_a
      this%last_a_ = last_a
   end subroutine set_last_a

   pure function last_a(this) result (res)
      class(singularity_correction), intent(in)  :: this
      integer                                    :: res
      res = this%last_a_
   end function last_a

   subroutine set_correction_parameters(this, sing_thresh, corr_thresh, corr_mode)
      class(singularity_correction), intent(inout) :: this
      real(r8), intent(in)                     :: sing_thresh, corr_thresh
      integer, intent(in)                    :: corr_mode
      this%sing_thresh_ = sing_thresh
      this%corr_thresh_ = corr_thresh
      this%corr_mode_   = corr_mode
   end subroutine set_correction_parameters

   subroutine correction_parameters(this, sing_thresh, corr_thresh, corr_mode)
      class(singularity_correction), intent(in)  :: this
      real(r8), intent(inout)                :: sing_thresh, corr_thresh
      integer, intent(inout)               :: corr_mode
      sing_thresh = this%sing_thresh_
      corr_thresh = this%corr_thresh_
      corr_mode   = this%corr_mode_
   end subroutine correction_parameters

   subroutine get_nearest_singularity(this, x, d, x_diff, x_sing, n_sing)
      class(singularity_correction) :: this
      real(r8), intent(in) :: x(3)       ! 3d position
      real(r8), intent(out) :: d          ! scaled distance
      real(r8), intent(out) :: x_diff(3)  ! distance vector x_sing - x
      real(r8), intent(out) :: x_sing(3)  ! 3d position of singularity
      integer, intent(out) :: n_sing     ! idx of singularity
      real(r8) xx_diff(3), dd
      integer i
      ! check vectorizability, forall? d as vector then minval?
      d = HUGE(d)
      do i = 1, size(this%singularities_, dim=2)
         xx_diff = this%singularities_(1:3, i) - x
         dd = NORM2(xx_diff)
         if (this%scaling_) then
            if (asserts) call assert(ALLOCATED(this%scalings_), &
                    'get_nearest_singularity: scaling, but scalings_ not allocated')
            dd = dd * this%scalings_(i)
         end if
         if (dd < d) then
            d = dd
            x_diff = xx_diff
            x_sing = this%singularities_(1:3, i)
            n_sing = i
         end if
      end do
   end subroutine get_nearest_singularity

   subroutine update_singularity_list(this, x, thresh, sp)
      ! slist contains index of singularity closer than thresh for each particle
      ! note: tests only particles with slist(n)==0 (i.e. not yet near singularity)
      class(singularity_correction), intent(in)  :: this
      real(r8), intent(in)                   :: x(:)
      real(r8), intent(in)                   :: thresh
      type(singularity_particles), intent(inout) :: sp
      real(r8) d(size(this%singularities_, dim=2))
      real(r8) xx_diff(3), thresh2
      integer i, n 

      if (asserts) call assert(mod(size(x),3) == 0 .and. size(sp%slist) == size(x)/3, &
       "update_singularity_list: illegal size")
      thresh2 = thresh * thresh
      do n = 1, size(x)/3
         if (sp%slist(n) > 0) cycle
         do i = 1, size(this%singularities_, dim=2)
            xx_diff = this%singularities_(1:3, i) - x(3*n-2:3*n)
            d(i) = dot_product(xx_diff, xx_diff)
         end do
         if (minval(d) < thresh2) sp%slist(n) = minloc(d, dim=1)
      end do
   end subroutine update_singularity_list

   subroutine set_to_singularity(this, x, sp, is_corrected, x_sing, idx_sing, i)
      class(singularity_correction), intent(in)  :: this
      real(r8), intent(inout) :: x(:)  ! current particle position vector
      type(singularity_particles), intent(inout) :: sp
      logical, intent(inout) :: is_corrected
      real(r8), intent(in) :: x_sing(3)
      integer, intent(in) :: idx_sing, i

      integer :: verbose, iul, spin

      verbose = this%verbose()
      iul = this%verbose_unit()

      if (i <= this%last_a_) then
         spin = 2  ! alpha
      else
         spin = 3  ! beta
      end if

      if (MOD(sp%nuc_list(idx_sing),spin) == 0) then
         is_corrected = .true.
         x(3*i-2:3*i) = x_sing   ! set to position of singularity
         sp%slist(i) = idx_sing
         sp%nuc_list(idx_sing) = sp%nuc_list(idx_sing) + spin + 1
         if (verbose > 3) write(iul,"(i6,a,i6)") i, " (a) set at singularity ", idx_sing
      else
         if (verbose > 2) write(iul,"(i6,a,i6)") i, " (a) NOT SET AT SINGULARITY", idx_sing
      end if

   end subroutine set_to_singularity

   subroutine umrigar_path_correction(r, v, tau, r_sing)
      ! particle position r, direction v with scaling tau (unmodified step: tau*v)
      ! vector z from nearest singular position r_sing
      ! notation like in Umrigar et al. JCP 1993, 99, 2865
      real(r8), intent(in) :: r(:)  ! r(3)!
      real(r8), intent(inout) :: v(:)  ! v(3)!
      real(r8), intent(inout) :: tau
      real(r8), intent(in) :: r_sing(:)

      real(r8) :: z(3), z0(3), norm_z, rho(3), rho0(3), v_rho, v_z, z11, rho11, r_new(3), norm_v, scaling

      if (asserts) call assert(size(r)==3 .and. size(v)==3 .and. size(r_sing)==3, &
         "umrigar_path_correction: illegal sizes")

      z = r - r_sing
      norm_z = NORM2(z)
      z0 = z / norm_z
      v_z = DOT_PRODUCT(v, z0)
      rho = v - v_z*z0
      v_rho = NORM2(rho)
      rho0 = rho / v_rho
      z11 = MAX(norm_z + v_z*tau, 0.d0)
      rho11 = 2 * v_rho * tau * z11 / (norm_z + z11)
      r_new = r_sing + rho11*rho0 + z11*z0

      ! determining new v and tau
      norm_v = NORM2(v)
      v = r_new - r
      ! making corrected v as long as uncorrected v and adapting tau
      scaling = norm_v / NORM2(v)
      v = v * scaling
      tau = tau / scaling
   end subroutine umrigar_path_correction

   subroutine cut_step_correction(r, v, tau, r_sing)
      ! particle position r, direction v with scaling tau (unmodified step: tau*v)
      ! vector z from nearest singular position r_sing
      ! cut step at nearest position to r_sing, i.e. prevent overshooting
      real(r8), intent(in) :: r(:)  ! r(3)!
      real(r8), intent(in) :: v(:)  ! v(3)!
      real(r8), intent(inout) :: tau
      real(r8), intent(in) :: r_sing(:)

      if (asserts) call assert(size(r)==3 .and. size(v)==3 .and. size(r_sing)==3, &
              "cut_step_path_correction: illegal sizes")

      tau = MIN( tau, - DOT_PRODUCT(r - r_sing, v) / DOT_PRODUCT(v, v) )
   end subroutine cut_step_correction

   subroutine correct_for_singularities(this, x, delta_x, sp, is_corrected, correction_only)
      ! updates x based on proposed step x_delta using singularity data
      ! updates the list of particles at singularity (slist, n_sing)
      class(singularity_correction), intent(in) :: this
      real(r8), intent(inout) :: x(:) ! current particle position vector
      real(r8), intent(in) :: delta_x(:) ! proposed step
      type(singularity_particles), intent(inout) :: sp
      logical, intent(inout) :: is_corrected
      logical, intent(in), optional :: correction_only ! if .true. do not move particles above corr_thresh

      real(r8) :: dd, x_diff(3), x_sing(3), tau(SIZE(x)/3), v(SIZE(delta_x))
      real(r8) :: tau_0 = 1._r8
      integer :: i, idx_sing, verbose, iul
      logical :: do_step

      if (asserts) call assert(MOD(SIZE(x),3) == 0 .and. SIZE(delta_x) == SIZE(x) .and. SIZE(sp%slist) == SIZE(x)/3, &
       "singularity_correction: illegal size")

      v = delta_x
      tau = tau_0
      !print*, "DBG:sing_corr"
      verbose = this%verbose()
      iul = this%verbose_unit()
      is_corrected = .false.
      do_step = .true.
      if (PRESENT(correction_only)) then
         do_step = .not. correction_only
      end if

      if (verbose > 4) write(iul,"(a,18i3)") " singularity correction with slist:", sp%slist
      do i = 1, size(x)/3
         if (sp%slist(i) /= 0) cycle  ! electron i is already placed at a singularity
         call this%get_nearest_singularity(x(3*i-2:3*i), dd, x_diff, x_sing, idx_sing)
         if (verbose > 4) write(iul,"(i4,a,i3,7g13.5)") i, " dist:", idx_sing, dd, x_diff, x_sing
         if (dd < this%sing_thresh_) then
            call this%set_to_singularity(x, sp, is_corrected, x_sing, idx_sing, i)
            cycle  ! electron i was just placed at a singularity
         end if
         if (dd < this%corr_thresh_) then
            if (MOD(this%corr_mode_,10) == UMRIGAR) then
               is_corrected = .true.
               call umrigar_path_correction(x(3*i-2:3*i), v(3*i-2:3*i), tau(i), x_sing)
               if (verbose > 3) write(iul,"(i6,a,2g13.5,i6)") i, " umrigar", dd, tau(i), idx_sing
            else if (MOD(this%corr_mode_,10) == CUTSTEP) then
               call cut_step_correction(x(3*i-2:3*i), v(3*i-2:3*i), tau(i), x_sing)
               if (verbose > 3) write(iul,"(i6,a,2g13.5,i6,l6)") i, " cutstep ", dd, tau(i), idx_sing
            end if
         end if
      end do

      is_corrected = is_corrected .or. ANY(tau /= tau_0)

      if (do_step) then
         if (.not. this%corr_mode_ < 10) then
            tau = MINVAL(tau)
         end if
         do i = 1, size(x)/3
            if (sp%slist(i) /= 0) cycle
            x(3*i-2:3*i) = x(3*i-2:3*i) + v(3*i-2:3*i) * tau(i)
            ! setting to singularities
            call this%get_nearest_singularity(x(3*i-2:3*i), dd, x_diff, x_sing, idx_sing)
            if (dd < this%sing_thresh_) then
               call this%set_to_singularity(x, sp, is_corrected, x_sing, idx_sing, i)
            end if
         end do
      end if
   end subroutine correct_for_singularities

   subroutine write_params_singularity_correction(this, iu)
      class(singularity_correction), intent(in) :: this
      integer, intent(in), optional :: iu
      integer iull

      if (present(iu)) then
         iull = iu
      else
         iull = this%verbose_unit()
      end if
      write(iull, "(a)") ""
      write(iull, "(a)") " singularity correction parameters:"
      write(iull, "(a,g13.5,a,l1)")  " sing_thresh         =", this%sing_thresh_
      write(iull, "(a,g13.5)")       " corr_thresh         =", this%corr_thresh_
      write(iull, "(a,i3)")          " corr_mode           =", this%corr_mode_
      write(iull, "(a,l3)")          " scaling             =" , this%scaling_
   end subroutine write_params_singularity_correction

end module singularityCorrection_m
