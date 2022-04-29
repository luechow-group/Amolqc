#
#  Makefile for amolqc/src
#

# include subdirectories
mindir := src/minimizers
include $(mindir)/rules.mk

dir := src
INC := -I$(utilsdir)

# common object files
SRCOBJ0_ := global_m.o machine_m.o init_m.o mainLoop_m.o subloop_m.o mpiInterface_m.o

# elocal and dependencies
SRCOBJ1_ := waveFunction_m.o wfData_m.o jastrow_m.o jastrowSM_m.o jastrowDTN_m.o jastrowIC_m.o jastrowParamData_m.o \
            jastrowAniso_m.o \
            aos_m.o aosData_m.o cuspOpt_m.o ecpIo_m.o ecp_m.o elocData_m.o eloc_m.o multiDet_m.o multiDetParam_m.o moParam_m.o \
            mos_m.o aoCut.o aoMoTask_m.o aoMo_m.o aoMoCut.o eConfigs_m.o linkedList_m.o \
            refUtils_m.o refBase_m.o refSimple_m.o refStr_m.o refStrPerm_m.o refVal_m.o refValStr_m.o refPos_m.o \
            refCtr_m.o refVector_m.o refADT_m.o refVList_m.o refListVList_m.o refDomain_m.o  \
            references_m.o rdataupdate_m.o sphericalIntGrids_m.o coulombDensity_m.o

# qmc
SRCOBJ2_ := qmc_m.o qmcSample_m.o propagator_m.o randomWalker_m.o  rWSample_m.o \
            rWStatistics_m.o elocTest_m.o reconf_m.o properties_m.o initialPositions_m.o

# optimization
SRCOBJ3_ := optimizeParams_m.o wfParameters_m.o optDerivsTest_m.o

SRCOBJ4_ := elocAndPsiTermsBase_m.o elocAndPsiTermsLin_m.o elocAndPsiTermsLM_m.o elocAndPsiTermsWEBFGS_m.o \
            elocAndPsiTermsENR_m.o elocAndPsiTermsEBFGS_m.o elocAndPsiTermsGen_m.o \
            optParamsBFGS_m.o optParamsELin_m.o optParamsENR_m.o \
            optParamsLBFGS_m.o optParamsVarmin_m.o optParamsVNL2SOL_m.o optParamsWBFGS_m.o optParamsPOpt_m.o

# psi^2 analysis
SRCOBJ5_ := psiMax_m.o assign_m.o findNucElecs_m.o energyPart_m.o \
            hungarian_m.o maxBasins_m.o maxAnalysis_m.o maxRawData_m.o  \
            fPsi2_m.o maximizeSample_m.o electronDensity_m.o rhoMax_m.o rhoData_m.o rhoGrid_m.o \
            maximizeSampleRho_m.o partition_m.o moMax_m.o posList_m.o vxc_m.o eigenVectAnalysis_m.o

FOPENMPTARGETS += $(dir)/aoMoTask_m.o

# add version pp-flag to init.f90
$(dir)/init_m.o: F90FLAGS += -DVERSION=\"$(VERSION)\"

# linker flags
LIBS := -L$(utilsdir) -lutils

SRCOBJ0 := $(addprefix $(dir)/,$(SRCOBJ0_))
SRCOBJ1 := $(addprefix $(dir)/,$(SRCOBJ1_))
SRCOBJ2 := $(addprefix $(dir)/,$(SRCOBJ2_))
SRCOBJ3 := $(addprefix $(dir)/,$(SRCOBJ3_))
SRCOBJ4 := $(addprefix $(dir)/,$(SRCOBJ4_))
SRCOBJ5 := $(addprefix $(dir)/,$(SRCOBJ5_))
SRCOBJS := $(SRCOBJ0) $(SRCOBJ1) $(SRCOBJ2) $(SRCOBJ3) $(SRCOBJ4) $(SRCOBJ5)


.PHONY: src
src: $(BINARY)

$(BINARY): $(dir)/amolqc_p.o $(SRCOBJS) $(MINOBJS) $(BASICSOBJS) $(utilsdir)/libutils.a
	$(FC) $(LDFLAGS) -o $@ $(dir)/amolqc_p.o $(SRCOBJS) $(MINOBJS) $(BASICSOBJS) $(LIBS) $(MLIBS)

.PHONY: clean-src
clean-src: clean-src-self clean-minimizers

.PHONY: clean-src-self
clean-src-self:
	@echo "Cleaning src"
	@-rm -f $(SRCOBJS)
	@-rm -f $(dir)/amolqc_p.o

# automatic dependency generation
$(dir)/depends.mk: $(dir)/*.f90
	$(PYTHONBIN) $(MODDEPENDS) -d -I$(utilsinc):$(mindir) $(<D) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dir)/depends.mk
endif
