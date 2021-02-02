# https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html

OMPFLAG = -fopenmp
AR = ar
RANLIB = ranlib
NOOPTFLAGS = -O0
DEBUGFLAGS = -fcheck=all
STDFLAGS = -O2
HIGHFLAGS = -O3
LINEINFOFLAGS = -g
WARNFLAGS = -Wall
NOWARNFLAGS = -Wno-all
DIALECT = -cpp -J$(MODDIR)
