.PHONY: all clean

all: spmv csx_llvm_tmpl.llvm.bc libcsx.so

## Configuration
llvmconf     = llvm-config-2.7

dynarray_dir = ../lib/dynarray
spm_dir      = ../lib/spm
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts
csx_template = $(shell pwd)/csx_llvm_tmpl.llvm.bc

GXX            ?= g++
CXXFLAGS       ?= -Wall -O3 -Wdisabled-optimization -DCSX_TEMPLATE=\"$(csx_template)\"
#CXXFLAGS       ?= -Wall -g -O0 -DCSX_TEMPLATE=\"$(csx_template)\"
INC            ?=

#CXXFLAGS       += -Winline
LDFLAGS         =
CXXFLAGS       += -ggdb -rdynamic -fPIC
DEFS            =
INC            += -I$(dynarray_dir) -I$(spm_dir) -I$(prfcnt_dir)

## cairomm
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS    = $(shell pkg-config cairomm-1.0 --libs)

#llvm
llvm_deps       = core analysis executionengine jit native bitreader ipo linker bitwriter
LLVM_GCC        = llvm-gcc-4.2
LLVM_CFLAGS     = $(shell $(llvmconf) --cflags)
LLVM_CXXFLAGS   = $(shell $(llvmconf) --cxxflags)
LLVM_LDFLAGS    = $(shell $(llvmconf) --ldflags --libs $(llvm_deps))

dynarray_dep    = $(dynarray_dir)/dynarray.o

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	DEFS += -DSPM_NUMA
	override LDFLAGS += -lnuma
endif

CXX_COMPILE     = $(GXX) $(CXXFLAGS) $(INC) $(DEFS)

# define __STDC_FORMAT_MACROS, so that we can use PRIu64
mmf.o: mmf.cc spm.h mmf.h
	$(CXX_COMPILE) -D__STDC_FORMAT_MACROS $< -c -o $@

spm.o: spm.cc spm.h spm_bits.h mmf.h
	$(CXX_COMPILE) $< -c -o $@

draw.o: draw.cc spm.h draw.h
	$(CXX_COMPILE) $(CAIROMM_FLAGS) $< -c -o $@

draw: draw.cc spm.h draw.h spm.o mmf.o
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $(CAIROMM_FLAGS) $(CAIROMM_LIBS) $(dynarray_dep) spm.o mmf.o $<  -o $@

csx.o: csx.cc csx.h spm.h delta.h
	$(CXX_COMPILE) -I$(dynarray_dir) $< -c -o $@

drle.o: drle.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

test.o: test.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

#drle_shell: mmf.o spm.o drle.o draw.o drle_shell.cc
#	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $(CAIROMM_LIBS) $^ -o $@

drle_test: mmf.o spm.o drle.o test.o drle_test.cc
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $^ -o $@

# define __STDC_FORMAT_MACROS, so that we can use PRIu64
csx_llvm_tmpl.llvm.bc: csx_llvm_tmpl.c ctl_ll.h
	$(LLVM_GCC) -D__STDC_FORMAT_MACROS -I../ $(INC) -Wall -emit-llvm -c csx_llvm_tmpl.c -o $@

llvm_jit_help.o: llvm_jit_help.cc llvm_jit_help.h
	$(GXX) $(CXXFLAGS) $(LLVM_CXXFLAGS) -c $<

jit.o: jit.cc jit.h spm.h drle.h llvm_jit_help.h
	$(CXX_COMPILE)  $(LLVM_CXXFLAGS) -c $<

spmv.o: spmv.cc spmv.h
	$(CXX_COMPILE) $(LLVM_CXXFLAGS) -c $<

elmerif.o: elmerif.cc elmerif.h spmv.h csr.h
	$(CXX_COMPILE) -c $<

main.o: main.cc spmv.h
	$(CXX_COMPILE) -c $<

# -Xlinker --allow-multiple-definition is set because
# dynarray code is also included in libspmv.o.
spmv: main.o spmv.o jit.o llvm_jit_help.o mmf.o spm.o csx.o drle.o $(spm_dir)/libspmv.o
	$(GXX) $(LDFLAGS) -ggdb -rdynamic -Xlinker --allow-multiple-definition $^ $(LLVM_LDFLAGS) $(dynarray_dep) -lpthread -o $@

libcsx.so: main.o spmv.o jit.o llvm_jit_help.o mmf.o spm.o csx.o drle.o elmerif.o $(spm_dir)/libspmv.o
	$(GXX) $(LDFLAGS) -ggdb -rdynamic -Xlinker --allow-multiple-definition -shared $^ $(LLVM_LDFLAGS) $(dynarray_dep) -lpthread -o $@

clean:
	rm -f *.o spmv
