# https://www.nag.co.uk/nagware/np/r62_doc/manual/compiler_2_4.html

OMPFLAG = -openmp
AR = ar
RANLIB = ranlib
NOOPTFLAGS = -O0
DEBUGFLAGS = -C -C=alias -C=dangling
STDFLAGS = -O2
HIGHFLAGS = -O3
LINEINFOFLAGS = -g -gline
WARNFLAGS = -w=x95 -w=x77
NOWARNFLAGS = -w=all
DIALECT = -fpp -mdir $(MODDIR) -I $(MODDIR) -DNAG -colour
