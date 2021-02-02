# https://software.intel.com/en-us/node/677967

OMPFLAG = -fopenmp
AR = xiar
RANLIB = xiar -s
NOOPTFLAGS = -O0
DEBUGFLAGS = -check all
STDFLAGS = -O2
HIGHFLAGS = -O3
FASTFLAGS = -fast -fp-model source
LINEINFOFLAGS = -g -traceback
WARNFLAGS = -warn all
NOWARNFLAGS = -warn none
DIALECT = -fpp -module $(MODDIR)

ifneq ($(FLAGS_INTEL_FAST_NO_FPOPT),)
  FASTFLAGS = $(FLAGS_INTEL_FAST_NO_FPOPT)
endif
