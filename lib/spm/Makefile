SHELL = /bin/bash
SUBDIRS = bcsr vbl

.PHONY: all clean $(SUBDIRS)

all: $(SUBDIRS) libspmv.o spmv
clean: $(SUBDIRS)

## Configuration
dynarray_dir = ../dynarray
prfcnt_dir   = ../prfcnt
script_dir   = ../../scripts

dynarray_dep = $(dynarray_dir)/dynarray.o

# This is ugly and non-portable
# (we need a script to return CPU MhZs)
MHZ_SH      = $(shell pwd)/$(prfcnt_dir)/scripts/cpu_mhz.sh
CC          = gcc
LD          = ld
CFLAGS      = -g -Wall -Winline -O3 -Wdisabled-optimization -fPIC
MAKE_LDFLAGS = 

# Always add the following flags
MAKE_CPPFLAGS = -DCPU_MHZ_SH=\"$(MHZ_SH)\"
MAKE_CPPFLAGS += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE -D__STDC_FORMAT_MACROS
MAKE_CPPFLAGS += -I$(dynarray_dir) -I$(prfcnt_dir)
ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	MAKE_CPPFLAGS += -DSPM_NUMA
endif

override CPPFLAGS += $(MAKE_CPPFLAGS)
override LDFLAGS += $(MAKE_LDFLAGS)

# Support for performance counters
# override CPPFLAGS    += -DCPU_CORE -DSPMV_PRFCNT # performance counters for an intel core
LIBS         = -lm -lpthread -lnuma
COMPILE      = $(CC) $(CPPFLAGS) $(CFLAGS)

#ifeq ($(CPU),NIAGARA)
#	CFLAGS += -m64
#	LD = ld -melf64_sparc
#endif
#
export CC
export CFLAGS
export MAKE_CPPFLAGS
export LD
export MAKE_LDFLAGS

$(SUBDIRS):
	 $(MAKE) -C $@ $(MAKECMDGOALS)

# do this, so that it doesn't match default %.c - >%.o rules
$(dynarray_dep):

deps         = $(dynarray_dep)
libspmv_deps = vector.o method.o mt_lib.o mmf.o  \
               spm_crs.o spm_crs_mt.o            \
               bcsr/libbcsr.o \
               vbl/libvbl.o \
               spmv_loops.o spmv_loops_mt.o spmv_loops_mt_numa.o \
               spmv_matvec_mt.o \
               numa_util.o

libspmv.o: $(libspmv_deps) $(deps) $(SUBDIRS)
	$(LD) -i --allow-multiple-definition $(libspmv_deps) $(deps) -o libspmv.o

vector.o: vector.c vector.h
	for t in double float; do                                 \
           $(COMPILE) -DELEM_TYPE=$$t  -c $< -o vector_$${t}.o;   \
        done
	$(LD) -i vector_{double,float}.o -o vector.o

mmf.o: mmf.c mmf.h
	$(COMPILE) -c $< -o $@

method.o: method.c method.h
	$(COMPILE) -c $< -o $@

mt_lib.o: mt_lib.c mt_lib.h
	$(COMPILE) -c $< -o $@

spm_crs.o:  spm_crs.c spm_crs.h
	for t in double float; do                          \
	  for ci in 32 64; do                              \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t \
	               -o spm_crs$${ci}_$${t}.o -c $< ;    \
	  done                                             \
	done
	$(LD) -i spm_crs{64,32}_{double,float}.o -o spm_crs.o

spm_crs_mt.o:  spm_crs_mt.c spm_crs_mt.h spm_mt.h
	for t in double float; do                             \
	  for ci in 32 64; do                                 \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t    \
	               -o spm_crs$${ci}_$${t}_mt.o -c $< ;    \
	  done                                                \
	done
	$(LD) -i spm_crs{64,32}_{double,float}_mt.o -o spm_crs_mt.o

spmv_loops.o: spmv_loops.c spmv_method.h vector.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o spmv_loops_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spmv_loops_double.o
	$(LD) -i spmv_loops_{float,double}.o -o spmv_loops.o

spmv_loops_mt.o: spmv_loops_mt.c spmv_loops_mt.h spmv_method.h vector.h spm_mt.h ../prfcnt/*.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o spmv_loops_mt_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spmv_loops_mt_double.o
	$(LD) -i spmv_loops_mt_{float,double}.o -o spmv_loops_mt.o

spmv_loops_mt_numa.o: spmv_loops_mt_numa.c spmv_loops_mt_numa.h spmv_method.h vector.h spm_mt.h ../prfcnt/*.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o spmv_loops_mt_numa_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spmv_loops_mt_numa_double.o
	$(LD) -i spmv_loops_mt_numa_{float,double}.o -o spmv_loops_mt_numa.o

spmv_matvec_mt.o: spmv_matvec_mt.c spmv_matvec_mt.h vector.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o spmv_matvec_mt_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spmv_matvec_mt_double.o
	$(LD) -i spmv_matvec_mt_{float,double}.o -o spmv_matvec_mt.o

spmv.o: spmv.c
	$(COMPILE) -c $< -o $@

%.o: %.c %.h
	$(COMPILE) -c $< -o $@

spmv: spmv.o libspmv.o
	$(CC) $(LDFLAGS) $(LIBS) $^ -o $@

## Usefull fallbacks
%.s: %.c
	$(COMPILE) -S -fverbose-asm $<
%.o: %.c
	$(COMPILE) -c $<
%.i: %.c
	$(COMPILE) -E $< | indent -kr > $@

clean:
	rm -f *.s *.o *.i spmv
