#
#  Makefile for amolqc/utils/ddda
#

dir := $(dddadir)

DDDAOBJ_ := boltzmann_m.o decorrelation_m.o parallelDDDA_m.o statistic_m.o

DDDAOBJS := $(addprefix $(dir)/,$(DDDAOBJ_))

.PHONY: clean-ddda
clean-ddda:
	@echo "  Cleaning ddda"
	@-rm -f $(DDDAOBJS)

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d -I$(basicsdir):$(optimizersdir):$(genericdir) $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
