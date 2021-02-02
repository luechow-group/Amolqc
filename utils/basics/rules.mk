#
#  Makefile for amolqc/utils/basics
#

dir := $(basicsdir)

BASICSOBJ_ := error_m.o globalUtils_m.o kinds_m.o verbosity_m.o

BASICSOBJS := $(addprefix $(dir)/,$(BASICSOBJ_))

.PHONY: clean-basics
clean-basics:
	@echo "  Cleaning basics"
	@-rm -f $(BASICSOBJS)

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
