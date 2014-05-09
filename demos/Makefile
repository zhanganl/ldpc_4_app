include ../Makefile.common

BINDIR = ../bin
LIBRARIES = -lc -lm -lpthread

LDPC_LIBRARY = $(BINDIR)/libldpc.a

CODE_FILES = simple_coder.c
DEC_FILES = simple_decoder.c
CODE_OBJ = $(BINDIR)/simple_coder
DEC_OBJ = $(BINDIR)/simple_decoder

all: $(CODE_OBJ) $(DEC_OBJ)

$(CODE_OBJ):$(CODE_FILES)
	@$(CC) $(CFLAGS) $(CODE_FILES) $(LIBRARIES) $(LDPC_LIBRARY) -o $(CODE_OBJ)
$(DEC_OBJ):$(DEC_FILES)
	@$(CC) $(CFLAGS) $(DEC_FILES) $(LIBRARIES) $(LDPC_LIBRARY) -o $(DEC_OBJ)

clean :
	@rm -rf *~

cleanall : clean
	@rm -rf $(CODE_OBJ) $(DEC_OBJ)
