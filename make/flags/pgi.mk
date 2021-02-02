# https://www.pgroup.com/resources/docs/17.10/x86/pgi-ref-guide/index.htm#cmdln-options-ref

OMPFLAG = -mp
AR = ar
RANLIB = ranlib
NOOPTFLAGS = -O0
DEBUGFLAGS = -C -Minfo
STDFLAGS = -O2
HIGHFLAGS = -O3
LINEINFOFLAGS = -g
DIALECT = -cpp -module $(MODDIR)
