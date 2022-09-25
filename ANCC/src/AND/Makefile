include ../../config.h

BIN = AND

LIBS = -lfftw3f


OBJS = db_m.o string.o sac_io.o date_time.o bindata_io.o math.o xcc.o AND.o butter.o


all: $(OBJS)
%.o: %.f90
	$(MPIFC) $(FFLAGS) $(LIBS) -c $< -o $@
%.o: %.for
	$(MPIFC) $(FFLAGS) $(LIBS) -c $< -o $@
%.o: %.f
	$(MPIFC) $(FFLAGS) $(LIBS) -c $< -o $@


install:
ifeq ($(wildcard $(INST_DIR)),)
	mkdir $(INST_DIR)
endif
	$(MPIFC) $(FFLAGS) $(OBJS) $(LIBS) -o $(INST_DIR)/$(BIN)

clean:
	rm -rf *.o *.mod *_genmod.f90

cleanall: clean
	rm -rf $(INST_DIR)/$(BIN)
