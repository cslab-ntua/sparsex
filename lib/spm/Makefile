SHELL = /bin/bash

.PHONY: all clean

all: libspmv.o spmv

## Configuration
debug_dir    = ../lib/debug
cpu_dir      = $(shell pwd)/../cpu
dynarray_dir = ../dynarray
phash_dir    = ../phash
prfcnt_dir   = ../prfcnt

dynarray_dep = $(dynarray_dir)/dynarray.o
phash_dep    = $(phash_dir)/phash.o

.PHONY: $(phash_dep) $(dynarray_dep)

#MHZ_SH      ?= $(cpu_dir)/cpu_mhz.sh
#CL_BYTES    ?= $(shell $(cpu_dir)/cl_bytes.sh)
#CACHE_BYTES ?= $(shell $(cpu_dir)/cache_bytes.sh)
#CPU         ?= $(shell $(cpu_dir)/cpu_info.sh)
GCC         ?= llvm-gcc
CFLAGS      ?= -Wall -Winline -O3 -Wdisabled-optimization -fPIC #-DCPU_CORE #-DSPMV_PRFCNT
CFLAGS      += -g
#DEFS         += -DCPU_MHZ_SH=\"$(MHZ_SH)\"
DEFS        += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
LIBS         = -lm -lpthread
INC          = -I$(dynarray_dir) -I$(phash_dir) -I$(prfcnt_dir)
COMPILE      = $(GCC) $(CFLAGS) $(INC) $(DEFS)
COMPILE_UR   = $(COMPILE) -funroll-loops
PYLIBS       = $(shell python2.5-config --ldflags)
PYCFLAGS     = $(shell python2.5-config --cflags)

ifeq ($(CPU),NIAGARA)
	CFLAGS += -m64
	LD = ld -melf64_sparc
endif

#ifeq ($(shell $(cpu_dir)/numa_lib.sh FOO),FOO)
#	DEFS += -DSPM_NUMA
#	LIBS += -lnuma
#endif

deps         = $(dynarray_dep) $(phash_dep)
libspmv_deps = vector.o method.o mt_lib.o mmf.o spm_parse.o \
               spm_crs.o spm_crs_mt.o spm_bcsr.o            \
               spm_csrdu_mul.o spm_csrdu.o                  \
               spm_crsvi_utils.o spm_crsvi.o spm_crsvi_mt.o \
               spm_csrdu_vi.o spm_csrdu_vi_mul.o            \
               spmv_loops.o spmv_loops_mt.o

libspmv.o: $(libspmv_deps) $(deps)
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

spm_parse.o: spm_parse.c spm_parse.h
	$(COMPILE) -c $< -o $@

spm_crs.o:  spm_crs.c spm_crs.h
	for t in double float; do                          \
	  for ci in 32 64; do                              \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t \
	               -o spm_crs$${ci}_$${t}.o -c $< ;    \
	  done                                             \
	done
	$(LD) -i spm_crs{64,32}_{double,float}.o -o spm_crs.o

spm_bcsr.o: spm_bcsr.c
	$(COMPILE) -DSPM_CRS_BITS=32 -DELEM_TYPE=double -o spm_bcsr32_double.o -c $<
	$(LD) -i spm_bcsr32_double.o -o spm_bcsr.o

spm_crs_mt.o:  spm_crs_mt.c spm_crs_mt.h
	for t in double float; do                             \
	  for ci in 32 64; do                                 \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t    \
	               -o spm_crs$${ci}_$${t}_mt.o -c $< ;    \
	  done                                                \
	done
	$(LD) -i spm_crs{64,32}_{double,float}_mt.o -o spm_crs_mt.o

spm_crsvi.o:  spm_crs_vi.c spm_crs_vi.h
	for t in double float; do                       \
	   for ci in 32 64; do                          \
	      for vi in 32 16 8; do                     \
	         $(COMPILE)                             \
		      -DELEM_TYPE=$$t                   \
		      -DSPM_CRSVI_CI_BITS=$$ci          \
		      -DSPM_CRSVI_VI_BITS=$$vi          \
		      -o spm_crs$${ci}_vi$${vi}_$${t}.o \
		      -c $<;                            \
	      done                                      \
	   done                                         \
	done
	$(LD) -i spm_crs{64,32}_vi{32,16,8}_{double,float}.o -o spm_crsvi.o

spm_crsvi_utils.o:  spm_crsvi_utils.c
	for t in double float; do                       \
		for vi in 32 16 8; do                   \
		  $(COMPILE)                            \
		    -DELEM_TYPE=$$t                     \
		    -DSPM_CRSVI_VI_BITS=$$vi            \
		    -o spm_crsvi$${vi}_$${t}_utils.o    \
		    -c $<;                              \
		done                                    \
	done
	$(LD) -i spm_crsvi{32,16,8}_{double,float}_utils.o -o spm_crsvi_utils.o

spm_csrdu_vi.o: spm_csrdu_vi.c
	for t in double float; do                       \
		for vi in 32 16 8; do                   \
		  $(COMPILE)                            \
		    -DELEM_TYPE=$$t                     \
		    -DSPM_CRSVI_VI_BITS=$$vi            \
		    -o spm_csrdu_vi$${vi}_$${t}.o       \
		    -c $<;                              \
		done                                    \
	done
	$(LD) -i spm_csrdu_vi{32,16,8}_{double,float}.o -o spm_csrdu_vi.o

spm_csrdu_vi_mul.o: spm_csrdu_vi_mul.c
	for t in double float; do                               \
		for vi in 32 16 8; do                           \
		  $(COMPILE)                                    \
		    -DELEM_TYPE=$$t                             \
		    -DSPM_CRSVI_VI_BITS=$$vi                    \
		    -o spm_csrdu_vi$${vi}_$${t}_mul.o           \
		    -c $<;                                      \
		  $(COMPILE)                                    \
		    -DCSRDU_ALIGNED                             \
		    -DELEM_TYPE=$$t                             \
		    -DSPM_CRSVI_VI_BITS=$$vi                    \
		    -o spm_csrdu_aligned_vi$${vi}_$${t}_mul.o   \
		    -c $<;                                      \
		done                                            \
	done
	$(LD) -i spm_csrdu{,_aligned}_vi{32,16,8}_{double,float}_mul.o \
	      -o spm_csrdu_vi_mul.o

spm_crsvi_mt.o:  spm_crs_vi_mt.c spm_crs_vi_mt.h
	for t in double float; do                          \
	   for ci in 32 64; do                             \
	      for vi in 32 16 8; do                        \
	         $(COMPILE)                                \
		      -DELEM_TYPE=$$t                      \
		      -DSPM_CRSVI_CI_BITS=$$ci             \
		      -DSPM_CRSVI_VI_BITS=$$vi             \
		      -o spm_crs$${ci}_vi$${vi}_mt_$${t}.o \
		      -c $<;                               \
	      done                                         \
	   done                                            \
	done
	$(LD) -i spm_crs{64,32}_vi{32,16,8}_mt_{double,float}.o -o spm_crsvi_mt.o

spmv_loops.o: spmv_loops.c spmv_method.h vector.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o spmv_loops_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spmv_loops_double.o
	$(LD) -i spmv_loops_{float,double}.o -o spmv_loops.o

spmv_loops_mt.o: spmv_loops_mt.c spmv_method.h vector.h ../prfcnt/*.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o spmv_loops_mt_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spmv_loops_mt_double.o
	$(LD) -i spmv_loops_mt_{float,double}.o -o spmv_loops_mt.o

spm_csrdu.o: spm_csrdu.c spm_csrdu.h spm_csrdu_mul.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spm_csrdu_double.o
	$(COMPILE) -DELEM_TYPE=float -c $< -o spm_csrdu_float.o
	$(LD) -i spm_csrdu_{double,float}.o spm_csrdu_mul.o -o spm_csrdu.o

spm_csrdu_mul.o: spm_csrdu_mul.c spm_csrdu.h
	$(COMPILE) -DELEM_TYPE=double -c $< -o spm_csrdu_mul_double.o
	$(COMPILE) -DELEM_TYPE=float -c $< -o spm_csrdu_mul_float.o
	$(COMPILE) -DCSRDU_ALIGNED -DELEM_TYPE=double -c $< -o spm_csrdu_aligned_mul_double.o
	$(COMPILE) -DCSRDU_ALIGNED -DELEM_TYPE=float -c $< -o spm_csrdu_aligned_mul_float.o
	$(LD) -i spm_csrdu{,_aligned}_mul_{double,float}.o -o spm_csrdu_mul.o

spmv.o: spmv.c
	$(COMPILE) -c $< -o $@

spmv: spmv.o libspmv.o
	$(COMPILE) $(LIBS) $^ -o $@
## Usefull fallbacks
%.s: %.c
	$(COMPILE) -S -fverbose-asm $<
%.o: %.c
	$(COMPILE) -c $<
%.i: %.c
	$(COMPILE) -E $< | indent -kr > $@

clean:
	rm -f *.s *.o *.i spmv
