#
#  Makefile for amolqc/utils/generic
#

dir := $(genericdir)

GENERICOBJ_ := blockAllocator_m.o genericFilter_m.o plateau_m.o

GENERICOBJS := $(addprefix $(dir)/,$(GENERICOBJ_))

.PHONY: clean-generic
clean-generic:
	@echo "  Cleaning generic"
	@-rm -f $(GENERICOBJS)

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d -I$(basicsdir) $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
