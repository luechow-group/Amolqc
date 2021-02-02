#
#  Makefile for amolqc/utils/statistics
#

dir := $(statisticsdir)

STATISTICSOBJ_ := intList_m.o statistics_m.o newStatistics_m.o vectorStatistics_m.o

STATISTICSOBJS := $(addprefix $(dir)/,$(STATISTICSOBJ_))

.PHONY: clean-statistics
clean-statistics:
	@echo "  Cleaning statistics"
	@-rm -f $(STATISTICSOBJS)

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d -I$(basicsdir) $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
