# Makefile for ANCC


############################################################################################################
#       Option to be set
############################################################################################################
# Select your compiler
MPIFC=mpiifort # MPI Fortran compiler binding with ifort [recommended]
#MPIFC=mpiifx   # MPI Fortran compiler binding with ifx
#MPIFC=mpifort  # It may bound to the Intel compiler or to the GNU compiler. Don't fooled by its name.
                # You should check it using mpifort --version
#MPIFC=mpif90   # MPI Fortran compiler bound to gfortran or intel
                # You should check it using mpif90 --version
                # It binds to the GNU compiler if the above command outputs GNU Fortran (GCC)
                # It binds to the Intel compiler if the above command outputs ifort (IFORT)
#MPIFC=pgf90    # PGI Fortran compiler with -Mmpi=mpich2
############################################################################################################



SUBDIRS = \
	ANCC \
	TF_PWS \
	AFTAN \
	$(EMPTY_MACRO)

SUBDIRS := $(addprefix ./src/,$(SUBDIRS))


all: config build clean

config:
	@./config.guess $(MPIFC)

build:
	$(call create_bindir)
	$(call COMPILE,$(SUBDIRS),$@)

clean:
	$(call COMPILE,$(SUBDIRS),$@)

cleanall: clean
	$(call COMPILE,$(SUBDIRS),$@)


define create_bindir:
include ./config.h
ifeq ($(wildcard $(INST_DIR)),)
	@mkdir -p $(DIST_DIR)
endif
endef
define COMMAND
$(MAKE) -C $(1) $(2)
endef
define COMPILE
$(foreach dir, $(1), $(call COMMAND,$(dir),$(2));)
endef
