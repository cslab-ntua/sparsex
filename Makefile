.PHONY: all clean

all: spmv_crs spmv_crsvi spmv_crs64
#all: spmv spmv-noxmiss dmv vxv spm_crsr_test
#all: spmv dmv vxv spmv_check spmv_lib.o

CACHE_BYTES ?= $(shell cpu/cache_bytes.sh)
CPU         ?= $(shell cpu/cpu_info.sh)
MHZ         ?= $(shell cpu/cpu_mhz.sh)
CL_BYTES    ?= $(shell cpu/cl_bytes.sh)
GCC         ?= gcc-4.2
CFLAGS      ?= -Wall -Winline -O3 -Wdisabled-optimization -fPIC
CFLAGS      += -g
#CFLAGS      += -funroll-all-loops #-march=nocona
DEFS        += -DCACHE_BYTES="$(CACHE_BYTES)" -DCL_BYTES=$(CL_BYTES)
DEFS        += -DCPU_$(CPU) -DCPU_MHZ=$(MHZ)
DEFS        += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
LIBS         = -lm
INC          = -Iprfcnt
COMPILE      = $(GCC) $(CFLAGS) $(INC) $(DEFS)
COMPILE_UR   = $(COMPILE) -funroll-loops
PYLIBS       = $(shell python2.5-config --ldflags)
PYCFLAGS     = $(shell python2.5-config --cflags) 

spmv_deps    = method.o mmf.o spm_parse.o spm_crs.o spm_delta.o spm_delta_vec.o #spmv_ur.o spm_crsr.o matrix.o
libspmv_deps = vector.o mmf.o method.o spm_parse.o spm_crs.o spm_crsvi.o spmv_loops.o spm_delta.o spm_delta_cv.o phash.o

vector.o: vector.c vector.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o vector_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o vector_double.o
	$(LD) -i vector_float.o vector_double.o -o vector.o

mmf.o: mmf.c mmf.h
	$(COMPILE) -c $< -o $@

method.o: method.c method.h 
	$(COMPILE) -c $< -o $@

spm_parse.o: spm_parse.c spm_parse.h 
	$(COMPILE) -c $< -o $@

spm_crs.o:  spm_crs.c spm_crs.h
	$(COMPILE) -DSPM_CRS_BITS=64 -DELEM_TYPE=float  -o spm_crs64_float.o -c $<
	$(COMPILE) -DSPM_CRS_BITS=32 -DELEM_TYPE=float  -o spm_crs32_float.o -c $<
	$(COMPILE) -DSPM_CRS_BITS=64 -DELEM_TYPE=double -o spm_crs64_double.o -c $<
	$(COMPILE) -DSPM_CRS_BITS=32 -DELEM_TYPE=double -o spm_crs32_double.o -c $<
	$(LD) -i spm_crs{64,32}_{double,float}.o -o spm_crs.o

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

spmv_loops.o: spmv_loops.c spmv_method.h vector.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o spmv_loops_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spmv_loops_double.o
	$(LD) -i spmv_loops_{float,double}.o -o spmv_loops.o

spm_delta.o: spm_delta_mul.c spm_delta.h vector.h
	$(COMPILE_UR) -DELEM_TYPE=float  -c $< -o spm_delta_mul_float.o
	$(COMPILE_UR) -DELEM_TYPE=double -c $< -o spm_delta_mul_double.o
	$(LD) -i spm_delta_mul_{float,double}.o -o spm_delta.o

spm_delta_cv.o: spm_delta_cv_mul.c spm_delta_cv.h spm_delta.h vector.h
	$(COMPILE_UR) -DELEM_TYPE=float  -c $< -o spm_delta_cv_mul_float.o
	$(COMPILE_UR) -DELEM_TYPE=double -c $< -o spm_delta_cv_mul_double.o
	$(LD) -i spm_delta_cv_mul_{float,double}.o -o spm_delta_cv.o

libspmv.o: $(libspmv_deps)
	$(LD) -i $(libspmv_deps) -o libspmv.o

dynarray.o: dynarray.c dynarray.h
	$(COMPILE) -c $< -o $@

phash.o: phash.c phash.h
	$(COMPILE) -c $< -o $@

ext_prog.o: ext_prog.c ext_prog.h
	$(COMPILE) -c $< -o $@

spmv_crsvi: libspmv.o dynarray.o spmv_crsvi.o
	$(COMPILE) libspmv.o dynarray.o spmv_crsvi.o -o $@

spmv_crs: libspmv.o dynarray.o spmv_crs.o
	$(COMPILE) libspmv.o dynarray.o spmv_crs.o -o $@

spmv_crs64: libspmv.o dynarray.o spmv_crs64.o
	$(COMPILE) libspmv.o dynarray.o spmv_crs64.o -o $@

cv_simple: cv_simple.o libspmv.o dynarray.o
	$(COMPILE_UR)  cv_simple.o libspmv.o dynarray.o -o $@

cv1: cv1.o libspmv.o dynarray.o
	$(COMPILE_UR)  cv1.o libspmv.o dynarray.o -o $@

cv1_d: cv1_d.o libspmv.o dynarray.o
	$(COMPILE_UR)  cv1_d.o libspmv.o dynarray.o -o $@

vals_idx: vals_idx.c
	$(COMPILE) $(PYCFLAGS) $(PYLIBS) $< -o $@

%.s: %.c
	$(COMPILE) -S -fverbose-asm $<
%.o: %.c
	$(COMPILE) -c $<
%.i: %.c
	$(COMPILE) -E $< | indent -kr > $@

clean:
	rm -rf *.s *.o *.i spmv_crs spmv_crs64 spmv_crsvi
