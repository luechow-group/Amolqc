#
# recursive Makefile for amolqc
#

# create mod directory
MODDIR = mod
$(shell mkdir -p $(MODDIR))

# create binary directory
BINDIR = bin
$(shell mkdir -p $(BINDIR))

# set binary
BINARY := $(BINDIR)/amolqc

VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)

# set default target to all
.PHONY: all
all: $(BINARY)

utilsdir := utils
FNOTARGETS :=
FOPENMPTARGETS :=

# user config
PYTHONBIN = python
MODDEPENDS = make/moduledepends.py

# default rng is MRG, alternative rng is MT (Mersenne-Twister)
RNG =

# include config
include make/config.mk

# rules
include make/vars.mk
include $(utilsdir)/rules.mk
include src/rules.mk

.PHONY: veryclean clean tests testp
veryclean: clean remove-config remove-testlog
clean: clean-src clean-utils clean-depends clean-mod clean-bin clean-testsuite clean-chmod
chmod: clean-chmod
tests: testserial
testp: testparallel

.PHONY: clean-mod clean-bin clean-depends clean-testsuite remove-config remove-testlog
clean-mod:
	@echo "Cleaning modules"
	@-rm -rf $(MODDIR)
clean-bin:
	@echo "Cleaning binary"
	@-rm -rf $(BINDIR)
clean-depends:
	@echo "Cleaning depends"
	@-find . -name "depends.mk" -exec rm {} \;
clean-testsuite:
	@echo "Cleaning testsuite"
	@-git clean -fd testsuite > /dev/null
clean-chmod:
	@echo "Resetting access permissions"
	@-chmod -x testsuite/run/tests.sh
	@-chmod -x make/testsuite.sh
	@-chmod -x make/moduledepends.py
	@-chmod +x configure.sh
	@-chmod +x tools/wfgen.py
	@-chmod +x tools/gamess2wf.rb
	@-chmod +x tools/gaussian2wf.py
remove-config:
	@echo "Removing config"
	@-rm -f make/config.mk
remove-testlog:
	@echo "Removing testlog"
	@-rm -f testsuite/run/test.log

test =

.PHONY: testserial testparallel
testserial:
	@bash ./make/testsuite.sh $(test)

testparallel:
	@bash ./make/testsuite.sh PARALLEL $(test)

# Remove implict rules
.SUFFIXES:

# add openmp flags for files that require it
$(FOPENMPTARGETS): F90FLAGS += $(LDFLAGS)
$(FOPENMPTARGETS): FNOFLAGS += $(LDFLAGS)

# Pattern rules for compilation
%.o: %.c
	gcc -std=c99 -O3 -c $< -o $@

%.o: %.f90
	$(FC) -c $(F90FLAGS) $(INC) $< -o $@

# Targets that should not be optimized
$(FNOTARGETS): %.o: %.f90
	$(FC) -c $(FNOFLAGS) $(INC) $< -o $@
