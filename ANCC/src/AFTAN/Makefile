include ../../config.h

BIN = AFTAN

LIBS = -lfftw3

OBJS = misc.o swapn.o taper.o dispio.o ftfilt.o fmax.o trigger.o tapers.o tgauss.o mspline.o pred_cur.o phtovel.o aftanpg.o aftanipg.o AFTAN.o


all: $(OBJS)
%.o: %.c
	$(CC) $(CFLAGS) $(LIBS) -c $< -o $@
%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

install:
ifeq ($(wildcard $(INST_DIR)),)
	mkdir $(INST_DIR)
endif
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(INST_DIR)/$(BIN)

clean:
	rm -rf $(OBJS) *.mod *__genmod.f90

cleanall:
	rm -rf $(INST_DIR)/$(BIN)
