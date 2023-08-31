include ../../config.h

BIN = ANCC

LIBS = -lfftw3f


OBJS = db_m.o string.o sac_io.o date_time.o bindata_io.o math.o xcc.o ANCC.o butter.o


# dependencies
ANCC.O: db_m.o xcc.o
math.o: db_m.o
sac_io.o: db_m.o
date_time.o: db_m.o
bindata_io.o: db_m.o sac_io.o
xcc.o: db_m.o math.o sac_io.o string.o date_time.o bindata_io.o


all: $(BIN)
$(BIN): $(OBJS)
	-rm -rf $(INST_DIR)/$(BIN)
	$(MPIFC) $(FFLAGS) $(OBJS) $(LIBS) -o $(INST_DIR)/$(BIN)

clean:
	rm -rf *.o *.mod *_genmod.f90

cleanall: clean
	rm -rf $(INST_DIR)/$(BIN)


%.o: %.f90
	$(MPIFC) $(FFLAGS) -c $< -o $@
%.o: %.for
	$(MPIFC) $(FFLAGS) -c $< -o $@
%.o: %.f
	$(MPIFC) $(FFLAGS) -c $< -o $@
