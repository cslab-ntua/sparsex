.PHONY: all clean

all: spmv

## Configuration
llvmconf     = llvm-config

dynarray_dir = ../lib/dynarray
spm_dir      = ../lib/spm
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts
mult_template_dir = $(shell pwd)/mult_templates
llvm_dir     = $(shell $(llvmconf) --includedir)
llvm_prefix  = $(shell $(llvmconf) --prefix)

CC          ?= gcc
CXX         ?= g++
LD           = $(CXX)
CFLAGS       = -Wall
CXXFLAGS    ?= -g -rdynamic -O3 -Wall -Woverloaded-virtual -Wdisabled-optimization -fPIC
LIBS         = -lnuma -lpthread -lboost_regex

# Always add the following flags
override CPPFLAGS += -DMULT_TEMPLATE_DIR=\"$(mult_template_dir)\" -DCLANG_PREFIX=\"$(llvm_prefix)\" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_GNU_SOURCE -D__STDC_FORMAT_MACROS
override CPPFLAGS += -I$(dynarray_dir) -I$(spm_dir) -I$(prfcnt_dir) -I$(llvm_dir)
override LDFLAGS  += -Wl,--allow-multiple-definition

## cairomm
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS    = $(shell pkg-config cairomm-1.0 --libs)

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
LINK        = $(LD) $(LDFLAGS)

# define __STDC_FORMAT_MACROS, so that we can use PRIu64
draw.o: draw.cc spm.h draw.h
	$(CXX_COMPILE) $(CAIROMM_FLAGS) $< -c -o $@

draw: draw.cc spm.h draw.h spm.o mmf.o
	$(LINK) $(CAIROMM_FLAGS) $(CAIROMM_LIBS) $(dynarray_dep) spm.o mmf.o $<  -o $@

csx.o: csx.cc csx.h spm.h delta.h ../lib/spm/numa_util.h ../lib/dynarray/dynarray.h
	$(CXX_COMPILE) -I$(dynarray_dir) $< -c -o $@

test.o: test.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

#drle_shell: mmf.o spm.o drle.o draw.o drle_shell.cc
#	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $(CAIROMM_LIBS) $^ -o $@

drle_test: mmf.o spm.o drle.o test.o drle_test.cc
	$(LINK) $^ -o $@

%.o: %.cc %.h
	$(CXX_COMPILE) -c $<

jit.o: spm.h csx.h drle.h template_text.h jit_util.h jit_config.h
compiler.o: jit_util.h
spmv.o: spm.h mmf.h csx.h drle.h jit.h ctl_ll.h ../lib/spm/mt_lib.h ../lib/spm/spmv_method.h ../lib/spm/spm_crs.h ../lib/spm/spm_mt.h ../lib/prfcnt/timer.h
spm.o: mmf.h ../lib/dynarray/dynarray.h
drle.o: spm.h

# -Xlinker --allow-multiple-definition is set because
# dynarray code is also included in libspmv.o.
spmv: spmv.o jit.o llvm_jit_help.o mmf.o spm.o csx.o drle.o template_text.o jit_util.o compiler.o $(spm_dir)/libspmv.o
	$(LINK) $(LIBS) $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(dynarray_dep) -o $@

newjit_test: newjit.o llvm_jit_help.o mmf.o spm.o csx.o drle.o $(dynarray_dep) $(spm_dir)/libspmv.o template_text.o jit_util.o compiler.o
	$(LINK) $(LIBS) $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) -o $@

clean:
	rm -f *.o spmv csx_llvm_tmpl.llvm.bc
