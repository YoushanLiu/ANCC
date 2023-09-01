include ../../config.h

BIN = AFTAN

LIBS = -lfftw3

OBJS = misc.o swapn.o taper.o dispio.o ftfilt.o fmax.o trigger.o tapers.o tgauss.o mspline.o pred_cur.o phtovel.o aftanpg.o aftanipg.o AFTAN.o


# dependencies
misc.o: swapn.o
dispio.o: misc.o


build: $(BIN)
$(BIN): $(OBJS)
	-rm -rf $(INST_DIR)/$(BIN)
ifeq ($(wildcard $(INST_DIR)),)
	mkdir $(INST_DIR)
endif
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(INST_DIR)/$(BIN)

clean:
	rm -rf $(OBJS) *.mod *__genmod.f90

cleanall:
	rm -rf $(INST_DIR)/$(BIN)


# C
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Fortran
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@
