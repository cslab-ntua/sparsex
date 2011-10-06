.PHONY: all clean

all: spmv csx_llvm_tmpl.llvm.bc

## Configuration
llvmconf     = llvm-config

dynarray_dir = ../lib/dynarray
spm_dir      = ../lib/spm
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts
csx_template = $(shell pwd)/csx_llvm_tmpl.llvm.bc
llvm_dir     = $(shell $(llvmconf) --includedir)

CXX         ?= g++
LD			 = $(CXX)
CXXFLAGS    ?= -g -rdynamic -Wall -O3 -Wdisabled-optimization
LIBS		 = -lnuma -lpthread

# Always add the following flags
override CPPFLAGS += -DCSX_TEMPLATE=\"$(csx_template)\" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS
override CPPFLAGS += -I$(dynarray_dir) -I$(spm_dir) -I$(prfcnt_dir) -I$(llvm_dir)
override LDFLAGS  += -Wl,--allow-multiple-definition

## cairomm
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS    = $(shell pkg-config cairomm-1.0 --libs)

#llvm
llvm_deps       = core analysis executionengine jit native bitreader ipo linker bitwriter
LLVM_CC         = llvm-gcc
LLVM_CFLAGS     = $(shell $(llvmconf) --cflags)
LLVM_CPPFLAGS	=
LLVM_CXXFLAGS   = 
LLVM_LDFLAGS    = $(shell $(llvmconf) --ldflags --libs $(llvm_deps))
LLVM_COMPILE	= $(LLVM_CC) $(LLVM_CPPFLAGS) $(CPPFLAGS) $(CXXFLAGS)
dynarray_dep    = $(dynarray_dir)/dynarray.o

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

CXX_COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS)
LINK        = $(LD) $(LDFLAGS)

# define __STDC_FORMAT_MACROS, so that we can use PRIu64
mmf.o: mmf.cc spm.h mmf.h
	$(CXX_COMPILE) -D__STDC_FORMAT_MACROS $< -c -o $@

spm.o: spm.cc spm.h mmf.h
	$(CXX_COMPILE) $< -c -o $@

draw.o: draw.cc spm.h draw.h
	$(CXX_COMPILE) $(CAIROMM_FLAGS) $< -c -o $@

draw: draw.cc spm.h draw.h spm.o mmf.o
	$(LINK) $(CAIROMM_FLAGS) $(CAIROMM_LIBS) $(dynarray_dep) spm.o mmf.o $<  -o $@

csx.o: csx.cc csx.h spm.h delta.h
	$(CXX_COMPILE) -I$(dynarray_dir) $< -c -o $@

drle.o: drle.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

test.o: test.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

#drle_shell: mmf.o spm.o drle.o draw.o drle_shell.cc
#	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $(CAIROMM_LIBS) $^ -o $@

drle_test: mmf.o spm.o drle.o test.o drle_test.cc
	$(LINK) $^ -o $@

# define __STDC_FORMAT_MACROS, so that we can use PRIu64
csx_llvm_tmpl.llvm.bc: csx_llvm_tmpl.c ctl_ll.h
	$(LLVM_COMPILE) -emit-llvm -c csx_llvm_tmpl.c -o $@

llvm_jit_help.o: llvm_jit_help.cc llvm_jit_help.h
	$(CXX_COMPILE) -c $<

jit.o: jit.cc jit.h spm.h drle.h llvm_jit_help.h
	$(CXX_COMPILE) -c $<

spmv.o: spmv.cc $(spm_dir)/spm_mt.h
	$(CXX_COMPILE) -c $<

# -Xlinker --allow-multiple-definition is set because
# dynarray code is also included in libspmv.o.
spmv: spmv.o jit.o llvm_jit_help.o mmf.o spm.o csx.o drle.o $(spm_dir)/libspmv.o
	$(LINK) $(LIBS) $^ $(LLVM_LDFLAGS) $(dynarray_dep) -o $@

clean:
	rm -f *.o spmv
