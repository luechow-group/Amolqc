! Copyright (C) 2013 Alexander Sturm
! Copyright (C) 2018 Arne Luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module random_m
  use mrg_m, only: init_mrgran, mrg_ran, mrg_gran, mrg_is_initialized
  use mt_m, only: init_mtran, mt_ran, mt_gran, mt_is_initialized

! module to make explicit all interfaces of subroutines and functions
! in the utils library (that are not included in modules)
! use random_m has to be included in routines that use routines from
! utils

  interface myran
#ifdef MT
    module procedure mt_ran
#else
    module procedure mrg_ran
#endif
  end interface myran

  interface mygran
#ifdef MT
    module procedure mt_gran
#else
    module procedure mrg_gran
#endif
  end interface mygran

  interface init_ran
#ifdef MT
    module procedure init_mtran
#else
    module procedure init_mrgran
#endif
  end interface init_ran

interface myran_is_initialized
#ifdef MT
    module procedure mt_is_initialized
#else
    module procedure mrg_is_initialized
#endif
end interface myran_is_initialized

end module
