################################################
#   Install directory
################################################

INST_DIR = ~/bin

################################################
#   Intel ifort
################################################

CC = icc
FC = ifort
MPIFC = mpiifort
CFLAGS = -O3 -xHost -I./
FFLAGS = -O3 -xHost -fpp -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -diag-disable 6477 -implicitnone -gen-interfaces -mcmodel=medium -shared-intel
##FFLAGS = -O3 -check nobounds -fpp -xAVX -ftz -assume buffered_io -assume byterecl -implicitnone -warn truncated_source -warn -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel

################################################
#   GNU gfortran
################################################

#CC = gcc
#FC = gfortran
#MPIFC = mpif90
#CFLAGS = -O3 -I./
#FFLAGS = -O3 -fimplicit-none -cpp
##FFLAGS = -O3 -std=gnu -fimplicit-none -frange-check -pedantic -pedantic-errors -Waliasing -Wampersand -Wline-truncation -Wsurprising -Wunderflow -fbounds-check

