.PHONY: all clean

all: cg

## Configuration
llvmconf     = llvm-config

dynarray_dir = ../lib/dynarray
spm_dir      = ../lib/spm
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts
csx_dir      = ../csx

CC          ?= gcc
CXX         ?= g++
CFLAGS       = -Wall
LD           = $(CXX)
CXXFLAGS    ?= -g -rdynamic -O3 -Wall -Woverloaded-virtual -Wdisabled-optimization -fPIC
LIBS	     = -lm -lnuma -lpthread -L~/spmv.new/spmv/csx/libcsx.so
OBJS	     = $(csx_dir)/spmv.o $(csx_dir)/jit.o $(csx_dir)/llvm_jit_help.o $(csx_dir)/mmf.o $(csx_dir)/spm.o $(csx_dir)/csx.o $(csx_dir)/drle.o $(csx_dir)/elmerif.o $(spm_dir)/libspmv.o

#llvm
llvm_deps       = core analysis executionengine jit native bitreader ipo linker bitwriter
LLVM_CC         = llvm-gcc
LLVM_CPPFLAGS	=
LLVM_CFLAGS     = -O0
LLVM_LDFLAGS    = $(shell $(llvmconf) --ldflags --libs $(llvm_deps))
LLVM_COMPILE	= $(LLVM_CC) $(CPPFLAGS) $(LLVM_CPPFLAGS) $(CFLAGS) $(LLVM_CFLAGS)
dynarray_dep    = $(dynarray_dir)/dynarray.o

# Always add the following flags
#override CPPFLAGS += -DCSX_TEMPLATE=\"$(csx_template)\" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_GNU_SOURCE -D__STDC_FORMAT_MACROS
override CPPFLAGS += -I$(dynarray_dir) -I$(spm_dir) -I$(prfcnt_dir) -I$(csx_dir)
override LDFLAGS  += -Wl,--allow-multiple-definition

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

CXX_COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS)
LINK        = $(LD) $(CPPFLAGS) $(LDFLAGS)
	 
cg_vector.o: cg_vector.cc cg_vector.h $(spm_dir)/vector.h
	$(CXX_COMPILE) $< -c -o $@

cg.o: cg.cc cg.h cg_vector.h
	$(CXX_COMPILE) $< -c -o $@
    
cg: main.cc cg.o cg_vector.o
	$(LINK) $(LIBS) $(OBJS) $^ $(LLVM_LDFLAGS) $(dynarray_dep) -o $@
    
clean:
	rm -f *.o cg
