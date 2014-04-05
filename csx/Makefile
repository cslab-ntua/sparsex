.PHONY: all clean

all: libcsximpl.a libcsximpl.so spmv element_test
#dynamic_array_test

# Configuration
llvmconf     = llvm-config

api_dir      = ../api
log_dir		 = logger
dynarray_dir = ../lib/dynarray
prfcnt_dir   = ../lib/prfcnt
script_dir   = ../scripts
mult_template_dir = $(shell pwd)/mult_templates
csx_prefix   = $(shell pwd)/..
llvm_dir     = $(shell $(llvmconf) --includedir)
llvm_prefix  = $(shell $(llvmconf) --prefix)
clang_prefix = $(llvm_prefix)
icu_dir		= /usr/lib
boost_dir	= /various/common_tools/boost-1.54

# This is ugly and non-portable
# (we need a script to return CPU MhZs)
MHZ_SH      = $(script_dir)/cpu_mhz.sh
CPU_TYPE    = $(shell $(script_dir)/cpu_info.sh)

CC          ?= gcc
CXX         ?= g++
LD           = $(CXX)
CFLAGS      ?= -Wall -O3 -fPIC
CXXFLAGS    ?= -g -rdynamic -O3 -Wall -Woverloaded-virtual -Wdisabled-optimization -fPIC -std=c++0x
LIBS         = -lnuma -lboost_thread -lboost_regex -lboost_serialization -lboost_system
AR           = ar
AR_FLAGS     = cru
MKDIR_P	     = mkdir -p

# Always add the following flags
override CPPFLAGS += -DCPU_MHZ_SH=\"$(MHZ_SH)\" -DCPU_$(CPU_TYPE)
override CPPFLAGS += -DMULT_TEMPLATE_DIR=\"$(mult_template_dir)\" -DCLANG_PREFIX=\"$(clang_prefix)\" -DCSX_PREFIX=\"$(csx_prefix)\" -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_GNU_SOURCE -D__STDC_FORMAT_MACROS
override CPPFLAGS += -I$(dynarray_dir) -I$(prfcnt_dir) -I$(llvm_dir) -I$(log_dir) -I$(api_dir) $(BOOST_INC_PATH)
override LDFLAGS  += -Wl,--allow-multiple-definition -Wl,--rpath=$(shell pwd) -pthread $(BOOST_LIB_PATH)

# cairomm
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS    = $(shell pkg-config cairomm-1.0 --libs)

# LLVM and Clang
llvm_deps       = core analysis executionengine jit native bitreader ipo linker bitwriter asmparser instrumentation
LLVM_CPPFLAGS	=
LLVM_LDFLAGS    = $(shell $(llvmconf) --ldflags --libs $(llvm_deps))
CLANG_LDFLAGS   = -Wl,--start-group -lclangAnalysis -lclangAST -lclangBasic -lclangCodeGen -lclangDriver -lclangFrontend -lclangFrontendTool -lclangIndex -lclangLex -lclangParse -lclangRewrite -lclangSema -lclangSerialization -lclangStaticAnalyzerCheckers -lclangStaticAnalyzerCore -lclangStaticAnalyzerFrontend -Wl,--end-group
dynarray_dep    = $(dynarray_dir)/dynarray.o
libadd = $(dynarray_dir)/dynarray.o numa_util.o
BOOST_INC_PATH	= -I$(boost_dir)/include
BOOST_LIB_PATH	= -L$(boost_dir)/lib 

DEPDIR = .deps
df = $(DEPDIR)/$(*F)
CSX_SRC = \
	Affinity.cpp \
	Compiler.cpp \
	CsxBench.cpp \
	CsxBuild.cpp \
	CsxMatvec.cpp \
	CsxUtil.cpp \
	CtlBuilder.cpp \
	DynamicArrayTest.cpp \
	ElementTest.cpp \
	Encodings.cpp \
	Facade.cpp \
	Jit.cpp \
	JitUtil.cpp \
	Node.cpp \
	Runtime.cpp \
	TemplateText.cpp \
	ThreadPool.cpp \
	Timer.cpp \
	Statistics.cpp \
	Vector.cpp \
	$(log_dir)/Logger.cpp \
	$(log_dir)/OutputPolicy.cpp \

CSX_OBJECTS = $(CSX_SRC:%.cpp=%.o) numa_util.o

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

CXX_COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS)
CC_COMPILE  = $(CC) $(CPPFLAGS) $(CFLAGS)
LINK        = $(LD) $(LDFLAGS)
ARCHIVE     = $(AR) $(AR_FLAGS)

%.o: %.c
	$(MKDIR_P) $(DEPDIR)
	$(CC_COMPILE) -c -MD -o $@ $<
	@cp $*.d $(df).P; \
    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	rm -f $*.d

%.o: %.cpp
	$(MKDIR_P) $(DEPDIR)
	$(CXX_COMPILE) -c -MD -o $@ $<
	@cp $*.d $(df).P; \
    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	rm -f $*.d

-include $(CSX_SRC:%.cpp=$(DEPDIR)/%.P)

libcsximpl.so: $(CSX_OBJECTS) $(libadd)
	$(LINK) -shared -o $@ $^ $(libadd) $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(LIBS)

libcsximpl.a: $(CSX_OBJECTS) $(libadd)
	$(ARCHIVE) $@ $^

spmv: main.o libcsximpl.so
	$(LINK) -o $@ $^ $(libadd) $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(LIBS)

set_test: set_test.o libcsximpl.so
	$(LINK) -o $@ $^ $(CLANG_LDFLAGS) $(LLVM_LDFLAGS) $(libadd) $(LIBS)

dynamic_array_test: DynamicArrayTest.o Timer.o $(libadd)
	$(LINK) -o $@ $^ $(libadd) -lnuma

element_test: ElementTest.o Encodings.o Statistics.o numa_util.o
	$(LINK) -o $@ $^ $(LIBS)

clean:
	rm -rf $(DEPDIR)
	rm -f *.o *.d $(log_dir)/*.o libcsximpl.a libcsximpl.so spmv element_test dynamic_array_test
