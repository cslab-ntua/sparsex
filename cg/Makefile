.PHONY: all clean

all: cg

## Configuration
llvmconf     = llvm-config

dynarray_dir = ../lib/dynarray
spm_dir      = ../lib/spm
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts
csx_dir      = ../csx
mult_template_dir = $(csx_dir)/mult_templates
llvm_dir     = $(shell $(llvmconf) --includedir)
llvm_prefix  = $(shell $(llvmconf) --prefix)

CC          ?= gcc
CXX         ?= g++
CFLAGS       = -Wall
LD           = $(CXX)
CXXFLAGS    ?= -g -rdynamic -O3 -Wall -Woverloaded-virtual -Wdisabled-optimization -fPIC
LIBS	     = -lm -lnuma -lpthread -lboost_regex -L~/spmv.new/spmv/csx/libcsx.so
OBJS	     = $(csx_dir)/spmv.o $(csx_dir)/compiler.o $(csx_dir)/jit_util.o $(csx_dir)/template_text.o $(csx_dir)/jit.o $(csx_dir)/llvm_jit_help.o $(csx_dir)/mmf.o $(csx_dir)/spm.o $(csx_dir)/csx.o $(csx_dir)/drle.o $(spm_dir)/libspmv.o

# Always add the following flags
override CPPFLAGS += -DCSX_TEMPLATE=\"$(mult_template_dir)\" -DCLANG_PREFIX=\"$(llvm_prefix)\" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_GNU_SOURCE -D__STDC_FORMAT_MACROS
override CPPFLAGS += -I$(dynarray_dir) -I$(spm_dir) -I$(prfcnt_dir) -I$(llvm_dir) -I$(csx_dir)
override LDFLAGS  += -Wl,--allow-multiple-definition

#llvm
llvm_deps       = core analysis executionengine jit native bitreader ipo linker bitwriter asmparser
LLVM_CC         = llvm-gcc
LLVM_CPPFLAGS	=
LLVM_CFLAGS     = -O0
LLVM_LDFLAGS    = $(shell $(llvmconf) --ldflags --libs $(llvm_deps))
LLVM_COMPILE	= $(LLVM_CC) $(CPPFLAGS) $(LLVM_CPPFLAGS) $(CFLAGS) $(LLVM_CFLAGS)
dynarray_dep    = $(dynarray_dir)/dynarray.o
CLANG_LDFLAGS   = -Wl,--start-group -lclangAnalysis -lclangAST -lclangBasic -lclangCodeGen -lclangDriver -lclangFrontend -lclangFrontendTool -lclangIndex -lclangLex -lclangParse -lclangRewrite -lclangSema -lclangSerialization -lclangStaticAnalyzerCheckers -lclangStaticAnalyzerCore -lclangStaticAnalyzerFrontend -Wl,--end-group

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

CXX_COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS)
LINK        = $(LD) $(CPPFLAGS) $(LDFLAGS)
	 
cg.o: cg.cc cg.h $(spm_dir)/vector.h
	$(CXX_COMPILE) $< -c -o $@
    
cg: main.cc cg.o $(OBJS)
	$(LINK) $(LIBS) $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(dynarray_dep) -o $@
    
clean:
	rm -f *.o cg
