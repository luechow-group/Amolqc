#
#  Makefile for amolqc/utils/random
#

dir := $(randomdir)

RANDOMOBJ_ := mrg_m.o mt_m.o random_m.o normal_m.o

RANDOMOBJS := $(addprefix $(dir)/,$(RANDOMOBJ_))

.PHONY: clean-random
clean-random:
	@echo "  Cleaning random"
	@-rm -f $(RANDOMOBJS)

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d -I$(basicsdir) $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
