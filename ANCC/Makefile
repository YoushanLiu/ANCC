# Makefile for ANCC

c



SUBDIRS = \
	AFTAN \
	ANCC \
	TF_PWS \
	$(EMPTY_MACRO)

SUBDIRS := $(addprefix ./src/,$(SUBDIRS))


all: config build

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
