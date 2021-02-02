#
#  Makefile for amolqc/utils/optimizers
#

dir := $(optimizersdir)

OPTIMIZERSOBJ_ := lbfgsb3_m.o nl2sol.o nl2sol_i.o

OPTIMIZERSOBJS := $(addprefix $(dir)/,$(OPTIMIZERSOBJ_))

.PHONY: clean-optimizers
clean-optimizers:
	@echo "  Cleaning optimizers"
	@-rm -f $(OPTIMIZERSOBJS)

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d -I$(basicsdir) $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
