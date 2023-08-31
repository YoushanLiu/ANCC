# Makefile for ANCC


############################################################################################################
#       Option to be set
############################################################################################################
# Select your compiler
MPIFC=mpiifort # MPI Fortran compiler binding with ifort [recommended]
#MPIFC=mpifort  # It may bound to the Intel compiler or to the GNU compiler. Don't fooled by its name.
                # You should check it using mpifort --version
#MPIFC=mpif90   # MPI Fortran compiler bound to gfortran or intel
                # You should check it using mpif90 --version
                # It binds to the GNU compiler if the above command outputs GNU Fortran (GCC)
                # It binds to the Intel compiler if the above command outputs ifort (IFORT)
#MPIFC=pgf90    # PGI Fortran compiler with -Mmpi=mpich2
############################################################################################################



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
