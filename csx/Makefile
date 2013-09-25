.PHONY: all clean

all: spmv libcsx.a libcsx.so dynamic_array_test rcm_test sparsematrix_test

# Configuration
llvmconf     = llvm-config

dynarray_dir = ../lib/dynarray
spm_dir      = ../lib/spm
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts
mult_template_dir = $(shell pwd)/mult_templates
csx_prefix   = $(shell pwd)/..
llvm_dir     = $(shell $(llvmconf) --includedir)
llvm_prefix  = $(shell $(llvmconf) --prefix)
clang_prefix = $(llvm_prefix)

CC          ?= gcc
CXX         ?= g++
LD           = $(CXX)
CFLAGS      ?= -Wall -O3 -fPIC
CXXFLAGS    ?= -g -rdynamic -O3 -Wall -Woverloaded-virtual -Wdisabled-optimization -fPIC
LIBS         = -lnuma -lboost_regex -lboost_thread -lboost_serialization
AR           = ar
AR_FLAGS     = cru

# Always add the following flags
override CPPFLAGS += -DMULT_TEMPLATE_DIR=\"$(mult_template_dir)\" -DCLANG_PREFIX=\"$(clang_prefix)\" -DCSX_PREFIX=\"$(csx_prefix)\" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_GNU_SOURCE -D__STDC_FORMAT_MACROS
override CPPFLAGS += -I$(dynarray_dir) -I$(spm_dir) -I$(prfcnt_dir) -I$(llvm_dir)
override LDFLAGS  += -Wl,--allow-multiple-definition -Wl,--rpath=$(shell pwd) -pthread

# cairomm
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS    = $(shell pkg-config cairomm-1.0 --libs)

# LLVM and Clang
llvm_deps       = core analysis executionengine jit native bitreader ipo linker bitwriter asmparser instrumentation
LLVM_CPPFLAGS	=
LLVM_LDFLAGS    = $(shell $(llvmconf) --ldflags --libs $(llvm_deps))
CLANG_LDFLAGS   = -Wl,--start-group -lclangAnalysis -lclangAST -lclangBasic -lclangCodeGen -lclangDriver -lclangFrontend -lclangFrontendTool -lclangIndex -lclangLex -lclangParse -lclangRewrite -lclangSema -lclangSerialization -lclangStaticAnalyzerCheckers -lclangStaticAnalyzerCore -lclangStaticAnalyzerFrontend -Wl,--end-group
dynarray_dep    = $(dynarray_dir)/dynarray.o
libadd = $(dynarray_dir)/dynarray.o $(spm_dir)/libspmv.o

DEPDIR = .deps
df = $(DEPDIR)/$(*F)
CSX_SRC = compiler.cc csx.cc drle.cc dynamic_array.cc encodings.cc \
jit.cc jit_util.cc main.cc runtime.cc spmv.cc template_text.cc timer.cc
CSX_OBJECTS = $(CSX_SRC:%.cc=%.o)

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

CXX_COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS)
CC_COMPILE  = $(CC) $(CPPFLAGS) $(CFLAGS)
LINK        = $(LD) $(LDFLAGS)
ARCHIVE     = $(AR) $(AR_FLAGS)

%.o: %.cc
	$(CXX_COMPILE) -c -MD -o $@ $<
	@cp $*.d $(df).P; \
    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	rm -f $*.d

-include $(CSX_SRC:%.cc=$(DEPDIR)/%.P)

# elmerif.o
libcsx.so: $(CSX_OBJECTS)
	$(LINK) -shared -o $@ $^ $(libadd) $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(LIBS)

libcsx.a: $(CSX_OBJECTS)
	$(ARCHIVE) $@ $^

spmv: main.o libcsx.so
	$(LINK) -o $@ $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(libadd) $(LIBS)

#libcsx_test: libcsx_test.o libcsx.so $(spm_dir)/libspmv.o
#	$(LINK) -o $@ $^ $(LIBS)

rcm_test: rcm_test.o libcsx.so
	$(LINK) -o $@ $^ $(LIBS)

sparsematrix_test: sparsematrix_test.o libcsx.so
	$(LINK) -o $@ $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(libadd) $(LIBS)

dynamic_array_test: dynamic_array.o timer.o
	$(LINK) -o $@ $^ $(libadd) -lnuma

clean:
	rm -f *.o *.d spmv libcsx.a libcsx.so libcsx_test dynamic_array_test rcm_test sparsematrix_test
