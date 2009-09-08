SHELL = /bin/bash

.PHONY: all clean

all: scripts/calc_bits spmv
#all: spmv_crs spmv_crsvi spmv_crs64 spmv_crsvi_check spmv_crs_mt spmv_crs_mt_check spmv_crsvi_mt spmv_crsvh spmv_crsvh_check spmv_crsvh_mt
#all: spmv spmv-noxmiss dmv vxv spm_crsr_test
#all: spmv dmv vxv spmv_check spmv_lib.o

bitutils_dir  = $(shell rsrc resource 'bitutils')
bitutils_dep  = $(bitutils_dir)/bitutils.o
dynarray_dir = $(shell rsrc resource 'dynarray')
dynarray_dep = $(dynarray_dir)/dynarray.o
phash_dir    = $(shell rsrc resource 'phash')
phash_dep    = $(phash_dir)/phash.o
huffman_dir  = $(shell rsrc resource 'huffman')
huffman_dep  = $(huffman_dir)/libhuffman.o
deps         = $(dynarray_dep) $(phash_dep) $(huffman_dep) $(bitutils_dep)

.PHONY: $(huffman_dep) $(phash_dep) $(bitutils_dep) $(dynarray_dep)

CACHE_BYTES ?= $(shell $(shell rsrc resource cache_bytes.sh))
CPU         ?= $(shell $(shell rsrc resource cpu_info.sh))
MHZ         ?= $(shell $(shell rsrc resource cpu_mhz.sh))
CL_BYTES    ?= $(shell $(shell rsrc resource cl_bytes.sh))
GCC         ?= gcc-4.3
CFLAGS      ?= -Wall -Winline -O3 -Wdisabled-optimization -fPIC
CFLAGS      += -g
#CFLAGS      += -funroll-all-loops #-march=nocona
DEFS        += -DCACHE_BYTES="$(CACHE_BYTES)" -DCL_BYTES=$(CL_BYTES)
DEFS        += -DCPU_$(CPU) -DCPU_MHZ=$(MHZ)
DEFS        += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
LIBS         = -lm -lpthread
INC          = -I$(shell rsrc resource 'prfcnt') -I$(dynarray_dir) -I$(phash_dir)
INC         += -I$(bitutils_dir) -I$(huffman_dir)
COMPILE      = $(GCC) $(CFLAGS) $(INC) $(DEFS)
COMPILE_UR   = $(COMPILE) -funroll-loops
PYLIBS       = $(shell python2.5-config --ldflags)
PYCFLAGS     = $(shell python2.5-config --cflags)

ifeq ($(CPU),NIAGARA)
	CFLAGS += -m64
	LD = ld -melf64_sparc
endif

ifeq ($(shell $(shell rsrc resource numa_lib.sh) FOO),FOO)
	DEFS += -DSPM_NUMA
	LIBS += -lnuma
endif

spmv_deps    = method.o mmf.o spm_parse.o spm_crs.o #spm_delta.o spm_delta_vec.o spmv_ur.o spm_crsr.o matrix.o
libspmv_deps = vector.o mmf.o method.o spm_parse.o spm_crs.o spm_crsvi.o spmv_loops.o spm_crs_mt.o spmv_loops_mt.o mt_lib.o spm_crsvi_mt.o spm_csrdu.o spm_crsvi_utils.o spm_csrdu_vi.o spm_csrdu_vi_mul.o spm_bcsr.o

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

spm_crsvh.o:  spm_crs_vh.c spm_crs_vh.h
	for t in double float; do                 \
	   for ci in 32 64; do                    \
	       $(COMPILE)                         \
	        -DELEM_TYPE=$$t                   \
	        -DSPM_CRSVH_CI_BITS=$$ci          \
	        -o spm_crs$${ci}_vh_$${t}.o       \
	        -c $<;                            \
	   done                                   \
	done
	$(LD) -i spm_crs{64,32}_vh_{double,float}.o -o spm_crsvh.o

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

spm_crsvh_mt.o:  spm_crs_vh_mt.c spm_crs_vh_mt.h
	for t in double float; do                 \
	   for ci in 32 64; do                    \
	       $(COMPILE)                         \
	        -DELEM_TYPE=$$t                   \
	        -DSPM_CRSVH_CI_BITS=$$ci          \
	        -o spm_crs$${ci}_vh_mt_$${t}.o    \
	        -c $<;                            \
		if [[ $$? != 0 ]]; then           \
			exit;                     \
		fi                                \
	   done                                   \
	done
	$(LD) -i spm_crs{64,32}_vh_mt_{double,float}.o -o spm_crsvh_mt.o

spmv_loops.o: spmv_loops.c spmv_method.h vector.h
	$(COMPILE) -DELEM_TYPE=float  -c $< -o spmv_loops_float.o
	$(COMPILE) -DELEM_TYPE=double -c $< -o spmv_loops_double.o
	$(LD) -i spmv_loops_{float,double}.o -o spmv_loops.o

spmv_loops_mt.o: spmv_loops_mt.c spmv_method.h vector.h
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

spm_csrdu_test.o: spm_csrdu_test.c spm_csrdu.h
	$(COMPILE) -c $< -o $@

spm_csrdu_test: spm_csrdu.o spm_csrdu_test.o $(dynarray_dep) mmf.o method.o mt_lib.o
	$(COMPILE) $^ -o $@

#spm_delta.o: spm_delta_mul.c spm_delta.h vector.h
#	$(COMPILE_UR) -DELEM_TYPE=float  -c $< -o spm_delta_mul_float.o
#	$(COMPILE_UR) -DELEM_TYPE=double -c $< -o spm_delta_mul_double.o
#	$(LD) -i spm_delta_mul_{float,double}.o -o spm_delta.o
#
#spm_delta_mt.o: spm_delta_mt_part.c spm_delta_mt.h spm_delta_mt_mul.c
#	$(COMPILE_UR) -DELEM_TYPE=float  -c spm_delta_mt_mul.c -o spm_delta_mt_mul_float.o
#	$(COMPILE_UR) -DELEM_TYPE=double  -c spm_delta_mt_mul.c -o spm_delta_mt_mul_double.o
#	$(COMPILE) -c spm_delta_mt_part.c -o spm_delta_mt_part.o
#	$(LD) -i spm_delta_mt_{mul_{float,double},part}.o -o spm_delta_mt.o
#
#
#spm_delta_cv.o: spm_delta_cv_mul.c spm_delta_cv.h spm_delta.h vector.h
#	$(COMPILE_UR) -DELEM_TYPE=float  -c $< -o spm_delta_cv_mul_float.o
#	$(COMPILE_UR) -DELEM_TYPE=double -c $< -o spm_delta_cv_mul_double.o
#	$(LD) -i spm_delta_cv_mul_{float,double}.o -o spm_delta_cv.o
#
libspmv.o: $(libspmv_deps)
	$(LD) -i --allow-multiple-definition $(libspmv_deps) $(deps) -o libspmv.o

mt_lib.o: mt_lib.c mt_lib.h
	$(COMPILE) -c $< -o $@

ext_prog.o: ext_prog.c ext_prog.h
	$(COMPILE) -c $< -o $@

spmv.o: spmv.c
	$(COMPILE) -c $< -o $@

spmv: spmv.o libspmv.o
	$(COMPILE) $(LIBS) $^ -o $@

spmv_crsvi: libspmv.o spmv_crsvi.o
	$(COMPILE) $(LIBS) libspmv.o spmv_crsvi.o -o $@

spmv_crsvi_mt: libspmv.o spmv_crsvi_mt.o
	$(COMPILE) $(LIBS) libspmv.o spmv_crsvi_mt.o -o $@

spmv_crsvi_check: libspmv.o spmv_crsvi_check.o
	$(COMPILE) $(LIBS) libspmv.o spmv_crsvi_check.o -o $@

spmv_crsvh_check: libspmv.o spmv_crsvh_check.o
	$(COMPILE) -Xlinker --allow-multiple-definition $(LIBS) spmv_crsvh_check.o $^  -o $@

spmv_crsvh: libspmv.o spmv_crsvh.o
	$(COMPILE) -Xlinker --allow-multiple-definition $(LIBS) spmv_crsvh.o $^  -o $@

spmv_crsvh_mt: libspmv.o spmv_crsvh_mt.o
	$(COMPILE) -Xlinker --allow-multiple-definition $(LIBS) $^  -o $@

spmv_crs: libspmv.o spmv_crs.o
	$(COMPILE) $(LIBS) libspmv.o spmv_crs.o -o $@

spmv_crs_mt: libspmv.o spmv_crs_mt.o
	$(COMPILE) $(LIBS) libspmv.o spmv_crs_mt.o -o $@

spmv_crs_mt_check: libspmv.o spmv_crs_mt_check.o
	$(COMPILE) $(LIBS) libspmv.o spmv_crs_mt_check.o -o $@

spmv_crs64: libspmv.o spmv_crs64.o
	$(COMPILE) $(LIBS) libspmv.o spmv_crs64.o -o $@

mmf_drle: mmf_drle.o mmf.o $(dynarray_dep)
	$(COMPILE) $^ -o $@

cv_simple: cv_simple.o libspmv.o
	$(COMPILE_UR)  cv_simple.o libspmv.o -o $@

cv1: cv1.o libspmv.o
	$(COMPILE_UR)  cv1.o libspmv.o -o $@

cv1_d: cv1_d.o libspmv.o
	$(COMPILE_UR)  cv1_d.o libspmv.o -o $@

vals_idx: vals_idx.c
	$(COMPILE) $(PYCFLAGS) $(PYLIBS) $< -o $@

scripts/calc_bits: calc_bits.c
	$(COMPILE) $< -o $@

%.s: %.c
	$(COMPILE) -S -fverbose-asm $<
%.o: %.c
	$(COMPILE) -c $<
%.i: %.c
	$(COMPILE) -E $< | indent -kr > $@


#rm -rf *.s *.o *.i spmv_crs{,64,vi{,_check,_mt},_mt,_mt_check} spmv_crsvh{,_check} spmv_crsvh_mt spm_csrdu_test spmv
clean:
	rm -rf *.s *.o *.i spmv crsvi_method
