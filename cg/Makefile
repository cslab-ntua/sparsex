SUBDIRS = spmv

.PHONY: all clean

all: cg

GCC=g++
CC=cc
CFLAGS=-Wall -O3 -Wdisabled-optimization #-funroll-all-loops -fopenmp -m64
DEFS=-DCPU_MHZ_SH=\"$(MHZ_SH)\" -DSPM_NUMA #-Dpar #-Dcheck #-DCPU_CORE -DCPU_BITS=64 -DPAUSE
OTHER_FLAGS=#CFLAGS+=$(DEFS) $(OTHER_FLAGS)

LIBDIR=
LIBS=-lm -lpthread -lnuma

dynarray_dir = ../lib/dynarray
spm_dir      = ../lib/spm
prfcnt_dir   = ../lib/prfcnt
csx_dir      = ../csx

INC= -I$(spm_dir) -I$(csx_dir) -I$(dynarray_dir)

SPM_OBJS=$(spm_dir)/libspmv.o
CSX_OBJS=$(csx_dir)/libcsx.so

OBJS=cg_vector.o cg.o $(SPM_OBJS) $(CSX_OBJS)
	 
cg_vector.o: cg_vector.cc cg_vector.h $(spm_dir)/libspmv.o
	$(CC) -L/$(LIBDIR) $(LIBS) $(CFLAGS) $(INC) -c $< -o $@
	
cg.o: cg.cc cg.h
	$(CC) -L/$(LIBDIR) $(LIBS) $(CFLAGS) $(INC) -c $< -o $@
	
cg: $(OBJS) main.cc
	$(GCC) $(OBJS) main.cc -o cg -L/$(LIBDIR) $(LIBS) $(INC) $(CFLAGS)
    
clean:
	rm -f *.o cg
