# compilers
ifeq ($(COMPILER),intel)
  FC = ifort
else ifeq ($(COMPILER),gnu)
  FC = gfortran
else ifeq ($(COMPILER),pgi)
  FC = pgfortran
else ifeq ($(COMPILER),nag)
  FC = nagfor
else
  $(error compiler not known)
endif

ifeq ($(MPI),yes)
  ifeq ($(MPICOMP),mpiifort)
    MPIFC = $(MPICOMP) -DMPI
  else ifeq ($(MPICOMP),mpifort)
    MPIFC = $(MPICOMP) -DMPI
  else ifeq ($(MPICOMP),mpif90)
    $(error mpif90 should no longer be used)
  else
    $(error mpi compiler not known)
  endif
  FC = $(MPIFC)
endif

# compile cache
ifdef CCACHE
  FC := $(CCACHE) $(FC)
endif

# compiler flags
include make/flags/$(COMPILER).mk

ifeq ($(WARNINGS),yes)
    DIALECT += $(WARNFLAGS)
else
    DIALECT += $(NOWARNFLAGS)
endif

ifeq ($(OPTLEVEL),no)
    OPTFLAGS = $(NOOPTFLAGS)
else ifeq ($(OPTLEVEL),debug)
    OPTFLAGS = $(DEBUGFLAGS)
else ifeq ($(OPTLEVEL),std)
    OPTFLAGS = $(STDFLAGS)
else ifeq ($(OPTLEVEL),high)
    OPTFLAGS = $(HIGHFLAGS)
else ifeq ($(OPTLEVEL),fast)
    OPTFLAGS = $(FASTFLAGS)
endif

FNOFLAGS += $(DIALECT) $(NOOPTFLAGS)
F90FLAGS += $(DIALECT) $(OPTFLAGS)

ifeq ($(MPI),yes)
    LDFLAGS += $(OMPFLAG)
endif

ifneq ($(RNG),)
FNOFLAGS += -D$(RNG)
F90FLAGS += -D$(RNG)
endif

ifeq ($(LINEINFO),yes)
  DIALECT += $(LINEINFOFLAGS)
else
  ifeq ($(OPTLEVEL),debug)
	 OPTFLAGS += $(LINEINFOFLAGS)
  endif
endif

ifeq ($(MLIBS),)
  ifeq ($(LAPACK),mkl)
    ifeq ($(MKLROOT),)
       $(error LAPACK mkl is given, but MKLROOT is not set)
    endif
    MKLLIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
    MKLPATH = $(MKLROOT)/lib/intel64
    MLIBS = -L$(MKLPATH) $(MKLLIBS)
  else ifeq ($(LAPACK),compiled)
    ifeq ($(MATHLIBS),)
       $(error LAPACK compiled is given, but MATHLIBS is not set)
    endif
    MLIBS = -L$(MATHLIBS) -llapack -lblas
  endif
endif
