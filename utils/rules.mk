#
#  Makefile for amolqc/utils
#

utilsinc := $(utilsdir)
UTILSOBJSSUB :=
# include subdirectories
basicsdir := $(utilsdir)/basics
include $(basicsdir)/rules.mk

genericdir := $(utilsdir)/generic
include $(genericdir)/rules.mk
UTILSOBJSSUB += $(GENERICOBJS)

randomdir := $(utilsdir)/random
include $(randomdir)/rules.mk
UTILSOBJSSUB += $(RANDOMOBJS)

optimizersdir := $(utilsdir)/optimizers
include $(optimizersdir)/rules.mk
UTILSOBJSSUB += $(OPTIMIZERSOBJS)

dddadir := $(utilsdir)/ddda
include $(dddadir)/rules.mk
UTILSOBJSSUB += $(DDDAOBJS)

statisticsdir := $(utilsdir)/statistics
include $(statisticsdir)/rules.mk
UTILSOBJSSUB += $(STATISTICSOBJS)

utilssub = $(basicsdir):$(randomdir):$(optimizersdir):$(dddadir):$(statisticsdir):$(genericdir)
utilsinc =$(utilsdir):$(utilssub)


# common object files
dir := $(utilsdir)
UTILSOBJ_ := cspline_m.o parsing_m.o atom_m.o \
             utils_m.o bins3D_m.o gcube_m.o sorting_m.o \
             numderivs_m.o string_utility_module.o compilerStrings_m.o \
             linAlg_m.o

# binaries
UTILSLIB_ := libutils.a

FNOTARGETS += $(dir)/lbfgsb3_m.o

UTILSOBJS := $(addprefix $(dir)/,$(UTILSOBJ_))
UTILSOBJS += $(UTILSOBJSSUB)
UTILSLIB := $(addprefix $(dir)/,$(UTILSLIB_))

.PHONY: libutils
libutils: $(UTILSLIB)

$(dir)/libutils.a: $(UTILSOBJS)
	$(AR) -r $@ $(UTILSOBJS)
# we have to call ranlib (i.e. ar -s) separately to work around a bug
# on OS X 10.7, where ar -s is a no-op
	$(RANLIB) $@

.PHONY: clean-utils
clean-utils: clean-utils-self clean-basics clean-random clean-optimizers clean-ddda clean-statistics

.PHONY: clean-utils-self
clean-utils-self:
	@echo "Cleaning utils"
	@-rm -f $(UTILSLIB)
	@-rm -f $(UTILSOBJS)

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d -I$(utilssub) $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
