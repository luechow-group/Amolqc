#
#  Makefile for amolqc/minimizers
#

dir := $(mindir)

MINOBJ_ :=  fctn_m.o line_search_simple_m.o line_search_weak_wolfe_m.o minimizer_m.o \
            minimizer_steep_desc_m.o minimizer_bfgs_m.o minimizer_fire_m.o minimizer_factory_m.o \
            singularityCorrection_m.o line_search_ws_simple_m.o minimizer_w_sing_m.o \
            minimizer_ws_steep_desc_m.o minimizer_ws_fire_m.o minimizer_ws_bfgs_m.o minimizer_ws_bfgst_m.o \
            minimizer_ws_newton_m.o minimizer_ws_factory_m.o singularityParticles_m.o line_search_m.o \
            minimizer_ws_none_m.o
MINOBJS := $(addprefix $(dir)/,$(MINOBJ_))

.PHONY: clean-minimizers
clean-minimizers:
	@echo "  Cleaning minimizers"
	@-rm -f $(MINOBJS)

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d -I$(utilsinc) $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
