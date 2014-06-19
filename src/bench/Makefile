.PHONY: all clean

all: bench
script_dir	= ../../scripts
sparsex_dir	= /usr/local
poski_dir	= $(poski_prefix)
mkl_dir		= $(mkl_prefix)

SPARSEXINC	= $(sparsex_dir)/include/sparsex
SPARSEXLIB	= $(sparsex_dir)/lib
POSKIINC	= $(poski_dir)/include/poski
POSKILIB	= $(poski_dir)/lib
OSKIINC		= $(poski_dir)/build_oski/include
OSKILIB		= $(poski_dir)/build_oski/lib/oski
MKLINC		= $(mkl_dir)/include
MKLLIB		= $(mkl_dir)/lib/intel64

CXX         	= g++
LD          	= g++
CXXFLAGS    	= -g -Wall -O3 -Wdisabled-optimization -Winline -fPIC
POSKI_CXXFLAGS	= -fopenmp
MKL_CXXFLAGS	= -fopenmp -m64
LIBS        	= -fopenmp -lpthread -ldl -lm
$(mkl_dir)/lib/intel64
SPARSEX_LDFLAGS = -Wl,-rpath,$(SPARSEXLIB) -L$(SPARSEXLIB) -lsparseximpl
POSKI_LDFLAGS   = -Wl,-rpath,$(POSKILIB) -L$(POSKILIB) -lposki -Wl,-rpath,$(OSKILIB) -L$(OSKILIB)
MKL_LDFLAGS     = -Wl,-rpath,$(MKLLIB) -L$(MKLLIB) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core

override CPPFLAGS += -I$(SPARSEXINC)
override LDFLAGS  += $(SPARSEX_LDFLAGS)

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS)
LINK    = $(LD)

BENCH_SRC = \
	Bench.cpp \
	Timer.cpp \
	SparsexModule.cpp \

BENCH_OBJECTS = $(BENCH_SRC:%.cpp=%.o)

ifneq ($(strip $(poski_prefix)),)
	override CPPFLAGS += -DPOSKI
	override CPPFLAGS += -I$(POSKIINC) -I$(OSKIINC)
	override LDFLAGS  += $(POSKI_LDFLAGS)
	BENCH_OBJECTS += poski_module.o
endif

ifneq ($(strip $(mkl_prefix)),)
	override CPPFLAGS += -DMKL
	override CPPFLAGS += -I$(MKLINC)
	override LDFLAGS  += $(MKL_LDFLAGS)
	BENCH_OBJECTS += mkl_module.o
endif

%.o: %.cpp %.hpp
	$(COMPILE) -c $< -o $@

PoskiModule.o: PoskiModule.cpp
ifneq ($(strip $(poski_prefix)),)
	$(COMPILE) $(POSKI_CXXFLAGS) -c $< -o $@
endif

MklModule.o: MklModule.cpp
ifneq ($(strip $(mkl_prefix)),)
	$(COMPILE) $(MKL_CXXFLAGS) -c $< -o $@
endif

bench: main.o $(BENCH_OBJECTS)
	$(LINK) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o bench
