SUBDIRS = spmv

.PHONY: all clean

all: cg

GCC=gcc
CC=cc
CFLAGS=-Wall -Wno-unused-variable -O3 #-funroll-all-loops -fopenmp -m64
DEFS=-DCPU_MHZ_SH=\"$(MHZ_SH)\" -DSPM_NUMA #-Dpar #-Dcheck #-DCPU_CORE -DCPU_BITS=64 -DPAUSE
OTHER_FLAGS=#CFLAGS+=$(DEFS) $(OTHER_FLAGS)

LIBDIR=
LIBS=-lblas -lm -lpthread -lnuma

dynarray_dir = ../lib/dynarray
spm_dir      = ../lib/spm
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts

INC= -I$(spm_dir)

SPM_OBJS=$(spm_dir)/libspmv.o
OBJS=cg_vector.o cg.o $(SPM_OBJS)
	 
cg_vector.o: cg_vector.c cg_vector.h
	$(CC) $(LIBS) $(CFLAGS) $(INC) -c $< -o $@
	
cg.o: cg.c cg.h
	$(CC) $(LIBS) $(CFLAGS) $(INC) -c $< -o $@
	
cg: $(OBJS) main.c
	$(GCC) $(OBJS) main.c -o cg -L/$(LIBDIR) $(LIBS) $(INC) $(CFLAGS)
    
clean:
	rm -f *.o cg
