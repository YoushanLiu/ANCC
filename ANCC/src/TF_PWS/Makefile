include ../../config.h

LIBS = -lfftw3f


ifeq ($(findstring gcc, $(CC)), gcc)
	LIBS += -lm
endif


BIN = TF_PWS

OBJS = sacrw.o st.o tf_pws.o



all: $(OBJS)
%.o: %.c
	$(CC) $(CFLAGS) $(LIBS) -c $< -o $@


install:
ifeq ($(wildcard $(INST_DIR)),)
	mkdir $(INST_DIR)
endif
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(INST_DIR)/$(BIN)


clean:
	rm -rf *.o


cleanall: clean
	rm -rf $(INST_DIR)/$(BIN)
