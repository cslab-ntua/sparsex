.PHONY: all clean

all: spmv_vhjit_check spmv_vhjit_mt spmv_vhjit_mt_check crsvh_spmv_template.llvm.bc

huffman_dir  = $(shell rsrc resource 'huffman')
phash_dir    = $(shell rsrc resource 'phash')
heap_dir     = $(shell rsrc resource 'heap')
dynarray_dir = $(shell rsrc resource 'dynarray')
bitutils_dir = $(shell rsrc resource 'bitutils')
prfcnt_dir   = $(shell rsrc resource 'prfcnt')

LD          ?=ld
LLVM_LD     ?=llvm-ld
GCC         ?= gcc
LLVM_GCC    ?= llvm-gcc
GXX         ?= g++
CFLAGS      ?= -Wall -Winline -O3 -Wdisabled-optimization -fPIC
CFLAGS      += -g
CXXFLAGS     =  -Wall -g
DEFS         = -D_GNU_SOURCE
INC          = -I../ -I$(huffman_dir)
INC         += -I$(phash_dir) -I$(heap_dir) -I$(dynarray_dir) -I$(bitutils_dir) -I$(prfcnt_dir)
COMPILE      = $(GCC) $(CFLAGS) $(INC) $(DEFS)

llvm_deps     = core analysis executionengine jit native bitreader ipo linker
LLVM_CFLAGS   = $(shell llvm-config --cflags)
LLVM_CXXFLAGS = $(shell llvm-config --cxxflags)
LLVM_LDFLAGS  = $(shell llvm-config --ldflags --libs $(llvm_deps))

LLVM_COMPILE = $(LLVM_GCC) -emit-llvm $(CFLAGS) $(INC) $(DEFS)

crsvh_spmv_template.llvm.bc: crsvh_spmv_template.c
	for t in double; do                                   \
		for ci in 32; do                              \
			$(LLVM_COMPILE)                       \
			 -DELEM_TYPE=$$t                      \
			 -DCI_BITS=$$ci                       \
			 -o crs$${ci}vh_spmv_$${t}_template.llvm.bc \
			 -c $<;                               \
		done                                          \
	done
	$(LLVM_LD) -link-as-library crs32vh_spmv_double_template.llvm.bc -o $@

crsvh.o: crsvh.c
	for t in double; do                                   \
		for ci in 32; do                              \
			$(COMPILE)                            \
			 -DELEM_TYPE=$$t                      \
			 -DCI_BITS=$$ci                       \
			 -o crs$${ci}_vh_$${t}.o              \
			 -c $<;                               \
		done                                          \
	done
	ld -i crs32_vh_double.o -o crsvh.o

vh_jit.o: vh_jit.cpp
	$(GXX) -DHUFF_LLVM_TEMPLATE='"$(huffman_dir)/huff_llvm_template.o"' $(INC) $(CXXFLAGS) $(LLVM_CXXFLAGS) -c $<

spmv_vhjit_check.o: spmv_vhjit_check.c
	$(GCC) $(INC) $(CFLAGS) -c $<

spmv_vhjit_mt.o: spmv_vhjit_mt.c
	$(GCC) $(INC) $(CFLAGS) -c $<

spmv_vhjit_mt_check.o: spmv_vhjit_mt_check.c
	$(GCC) $(INC) $(CFLAGS) -c $<

spmv_vhjit_check: spmv_vhjit_check.o crsvh.o vh_jit.o
	$(GXX) $^ ../libspmv.o $(huffman_dir)/libhj_llvm.o $(LLVM_LDFLAGS) -o $@

spmv_vhjit_mt: spmv_vhjit_mt.o crsvh.o vh_jit.o
	$(GXX) $^ ../libspmv.o $(huffman_dir)/libhj_llvm.o $(LLVM_LDFLAGS) -o $@

spmv_vhjit_mt_check: spmv_vhjit_mt_check.o crsvh.o vh_jit.o
	$(GXX) $^ ../libspmv.o $(huffman_dir)/libhj_llvm.o $(LLVM_LDFLAGS) -o $@

%.s: %.c
	$(COMPILE) -S -fverbose-asm $<
%.o: %.c
	$(COMPILE) -c $<
%.i: %.c
	$(COMPILE) -E $< | indent -kr > $@
%.i: %.cpp
	$(COMPILE) $(CXXFLAGS) $(LLVM_CXXFLAGS) -E $<  > $@


clean:
	rm -f *.llvm.bc *.o *.i spmv_vhjit_check spmv_vhjit_mt{,_check}
