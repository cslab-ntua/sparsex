.PHONY: all clean

all: libcsx_in.a libcsx_in.so sparsematrix_test
# libcsx_test set_test spmv

# Configuration
llvmconf     = llvm-config

API_dir      = ../C-API
log_dir		 = logger
dynarray_dir = ../lib/dynarray
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts
mult_template_dir = $(shell pwd)/mult_templates
csx_prefix   = $(shell pwd)/..
llvm_dir     = $(shell $(llvmconf) --includedir)
llvm_prefix  = $(shell $(llvmconf) --prefix)
clang_prefix = $(llvm_prefix)

# This is ugly and non-portable
# (we need a script to return CPU MhZs)
MHZ_SH      = $(script_dir)/cpu_mhz.sh
CPU_TYPE    = $(shell $(script_dir)/cpu_info.sh)

CC          ?= gcc
CXX         ?= g++
LD           = $(CXX)
CFLAGS      ?= -Wall -O3 -fPIC
CXXFLAGS    ?= -g -rdynamic -O3 -Wall -Woverloaded-virtual -Wdisabled-optimization -fPIC
LIBS         = -lnuma -lboost_regex -lpthread -lboost_thread -lboost_serialization
AR           = ar
AR_FLAGS     = cru

# Always add the following flags
override CPPFLAGS += -DCPU_MHZ_SH=\"$(MHZ_SH)\" -DCPU_$(CPU_TYPE)
override CPPFLAGS += -DMULT_TEMPLATE_DIR=\"$(mult_template_dir)\" -DCLANG_PREFIX=\"$(clang_prefix)\" -DCSX_PREFIX=\"$(csx_prefix)\" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_GNU_SOURCE -D__STDC_FORMAT_MACROS
override CPPFLAGS += -I$(dynarray_dir) -I$(prfcnt_dir) -I$(llvm_dir) -I$(log_dir)
override LDFLAGS  += -Wl,--allow-multiple-definition -Wl,--rpath=$(shell pwd)

# cairomm
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS    = $(shell pkg-config cairomm-1.0 --libs)

# LLVM and Clang
llvm_deps       = core analysis executionengine jit native bitreader ipo linker bitwriter asmparser instrumentation
LLVM_CPPFLAGS	=
LLVM_LDFLAGS    = $(shell $(llvmconf) --ldflags --libs $(llvm_deps))
CLANG_LDFLAGS   = -Wl,--start-group -lclangAnalysis -lclangAST -lclangBasic -lclangCodeGen -lclangDriver -lclangFrontend -lclangFrontendTool -lclangIndex -lclangLex -lclangParse -lclangRewrite -lclangSema -lclangSerialization -lclangStaticAnalyzerCheckers -lclangStaticAnalyzerCore -lclangStaticAnalyzerFrontend -Wl,--end-group
dynarray_dep    = $(dynarray_dir)/dynarray.o

CSX_OBJECTS = $(log_dir)/logger.o $(log_dir)/output_policy.o timer.o csx_util.o affinity.o csx_bench.o csx_spmv_mt.o csxsym_spmv_mt.o csx_matvec.o vector.o numa_util.o sparse_matrix_wrapper.o timer.o node.o runtime.o csx_build.o jit.o template_text.o jit_util.o compiler.o dynamic_array.o

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

CXX_COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS)
CC_COMPILE  = $(CC) $(CPPFLAGS) $(CFLAGS)
LINK        = $(LD) $(LDFLAGS)
ARCHIVE     = $(AR) $(AR_FLAGS)

# precompiled.h.gch: precompiled.h
# 	$(CXX_COMPILE) $< -o $@
# CXXFLAGS += -include precompiled.h

#draw.o: draw.cc spm.h draw.h
#	$(CXX_COMPILE) $(CAIROMM_FLAGS) $< -c -o $@

#draw: draw.cc sparsepartition.h sparsepartition_impl.hpp draw.h
#	$(LINK) $(dynarray_dep) $<  -o $@ $(CAIROMM_FLAGS) $(CAIROMM_LIBS) 

sparse_matrix_wrapper.o: sparse_matrix_wrapper.cc sparse_matrix_wrapper.h sparse_matrix.h sparse_internal.h sparse_partition.h $(API_dir)/mattype.h csx_save_restore.h timer.h
	$(CXX_COMPILE) -I$(API_dir) $< -c -o $@

test.o: test.cc drle.h sparse_parition.h
	$(CXX_COMPILE) $< -c -o $@

drle_test: mmf.o drle.o test.o drle_test.cc
	$(LINK) $^ -o $@

%.o: %.cc %.h
	$(CXX_COMPILE) -c $< $(LIBS)

%.o: %.c %.h
	$(CC_COMPILE) -c $<

#elmerif.o: elmerif.cc elmerif.h csx_build.h csr.h sparse_internal.h sparse_partition.h
#	$(CXX_COMPILE) -DSPM_CRS_BITS=32 -c $<
csx_util.o: $(log_dir)/logger.hpp
runtime.o: csx.h sparse_partition.h jit.h spm_mt.h affinity.h $(log_dir)/logger.hpp
jit.o: sparse_partition.h csx.h drle.h template_text.h jit_util.h jit_config.h $(log_dir)/logger.hpp
compiler.o: jit_util.h
csx_build.o: sparse_partition.h runtime.h csx.h drle.h jit.h ctl_ll.h spm_mt.h affinity.h spmv_method.h $(log_dir)/logger.hpp timer.h
elmerif.o: csr.h csx_build.h
main.o: csx_build.h matrix_loading.h runtime.h sparse_internal.h sparse_partition.h mmf.h csx.h csx_get_set.h $(log_dir)/logger.hpp
#libcsx_test.o: elmerif.h $(spm_dir)/spm_crs.h
sparsematrix_test.o: csx_build.h rcm.h runtime.h sparse_internal.h sparse_partition.h sparse_matrix.h mmf.h csr.h csx.h $(log_dir)/logger.hpp
dynamic_array.o: allocators.h $(dynarray_dir)/dynarray.h

# elmerif.o
libcsx_in.so: $(CSX_OBJECTS)
	$(LINK) -shared -o $@ $^ $(dynarray_dep) $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(LIBS)

libcsx_in.a: $(CSX_OBJECTS)
	$(ARCHIVE) $@ $^

spmv: main.o libcsx_in.so
	$(LINK) -o $@ $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(dynarray_dep) $(log_dep) $(LIBS)

libcsx_test: libcsx_test.o libcsx.so $(spm_dir)/libspmv.o
	$(LINK) -o $@ $^ $(LIBS)

sparsematrix_test: sparsematrix_test.o libcsx_in.so
	$(LINK) -o $@ $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(dynarray_dep) $(LIBS)

set_test: set_test.o libcsx_in.so
	$(LINK) -o $@ $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(dynarray_dep) $(LIBS)

dynamic_array_test: dynamic_array.o timer.o
	$(LINK) -o $@ $^ $(dynarray_dep) -lnuma

clean:
	rm -f *.o spmv libcsx_in.a libcsx_in.so libcsx_test dynamic_array_test sparsematrix_test
