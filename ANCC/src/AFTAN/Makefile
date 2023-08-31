# Makefile for ANCC

# Select your compiler
MPIFC=mpiifort # MPI Fortran compiler binding with ifort [recommended]
#MPIFC=mpif90   # MPI Fortran compiler binding with gfortran or intel
                # you should check mpif90 --version
                # It binds with the GNU compiler if it outputs GNU Fortran (GCC)
                # It binds with the Intel compiler if it outputs ifort (IFORT)
#MPIFC=pgf90    # PGI Fortran compiler with -Mmpi=mpich2



SUBDIRS = \
	AFTAN \
	ANCC \
	TF_PWS \
	$(EMPTY_MACRO)

SUBDIRS := $(addprefix ./src/,$(SUBDIRS))


。PHONY: config build

config:
	@./config.guess $(MPIFC)
include ./config.h
ifeq ($(wildcard $(INST_DIR)),)
	mkdir $(INST_DIR)
endif

define BUILD
$(MAKE) -C $(1) all
endef
define CLEAN
$(MAKE) -C $(1) clean
endef
define CLEANALL
$(MAKE) -C $(1) cleanall
endef


build:
	$(foreach dir, $(SUBDIRS), $(call BUILD,$(dir));)

clean:
	$(foreach dir, $(SUBDIRS), $(call CLEAN,$(dir));)

cleanall: clean
	$(foreach dir, $(SUBDIRS), $(call CLEANALL,$(dir));)
