include ../../config.h

LIBS = -lfftw3f


ifeq ($(findstring gcc, $(CC)), gcc)
	LIBS += -lm
endif


BIN = TF_PWS

OBJS = sacrw.o st.o tf_pws.o


# dependencies
tf_pws.o: sacrw.o st.o


build: $(BIN)
$(BIN): $(OBJS)
	-rm -rf $(INST_DIR)/$(BIN)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(INST_DIR)/$(BIN)


clean:
	rm -rf *.o

cleanall: clean
	rm -rf $(INST_DIR)/$(BIN)


%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
