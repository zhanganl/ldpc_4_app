include ../Makefile.common

LIBRARIES = -lc -lm -lpthread 
BINDIR = ../bin
LIB_OBJ = $(BINDIR)/libldpc.a

SRCFILES  = ldpc_create_pchk.c ldpc_fec.c ldpc_fec_iterative_decoding.c ldpc_matrix_sparse.c ldpc_group.c
OFILES = $(SRCFILES:.c=.o)

all: lib

%.o:%.c
	@$(CC) -c $(CFLAGS) -o $@ $<

lib:$(OFILES)
	@$(AR) -rcsv $(LIB_OBJ) $(OFILES)

clean :
	@rm -rf *~
	@rm -rf $(OFILES)

cleanall : clean
	@rm -rf $(LIB_OBJ)
