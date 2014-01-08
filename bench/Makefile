.PHONY: all clean

all: bench
poski_dir	= $(poski_prefix)
mkl_dir		= $(mkl_prefix)

script_dir	= ../scripts
log_dir		= ../csx/logger
csx_dir     = ../csx
sparsex_dir	= ../api
POSKIINC	= $(poski_dir)/include/poski
POSKILIB	= $(poski_dir)/lib
OSKIINC		= $(poski_dir)/build_oski/include
OSKILIB		= $(poski_dir)/build_oski/lib/oski

CXX         = g++
LD          = g++
CXXFLAGS    = -g -Wall -O3 -Wdisabled-optimization -Winline -fPIC
MKL_CXXFLAGS	= -fopenmp -m64
POSKI_CXXFLAGS	= -fopenmp
LIBS        	= -fopenmp -lpthread -ldl -lm

SPARSEX_LDFLAGS = -L$(sparsex_dir) -Wl,--rpath=$(sparsex_dir) -lsparsex -L$(csx_dir) -Wl,--rpath=$(csx_dir) -lcsximpl
MKL_LDFLAGS     = -Wl,--start-group $(mkl_dir)/lib/intel64/libmkl_intel_lp64.a $(mkl_dir)/lib/intel64/libmkl_gnu_thread.a $(mkl_dir)/lib/intel64/libmkl_core.a -Wl,--end-group
POSKI_LDFLAGS   = -Wl,--rpath -Wl,$(POSKILIB) -L$(POSKILIB) -lposki -Wl,--rpath -Wl,$(OSKILIB) -L$(OSKILIB)

override CPPFLAGS += -I$(sparsex_dir) -I$(csx_dir) -I$(log_dir) -I$(mkl_dir)/include -I$(POSKIINC) -I$(OSKIINC)
override LDFLAGS  += $(SPARSEX_LDFLAGS)
#$(MKL_LDFLAGS) $(POSKI_LDFLAGS)
override LDFLAGS  += -Wl,--allow-multiple-definition

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS)
LINK    = $(LD)

%.o: %.cc %.h
	$(COMPILE) -c $< -o $@

poski_module.o: poski_module.cc
	$(COMPILE) $(POSKI_CXXFLAGS) -c $< -o $@

mkl_module.o: mkl_module.cc
	$(COMPILE) $(MKL_CXXFLAGS) -c $< -o $@

bench: main.o bench.o timer.o sparsex_module.o
	$(LINK) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o bench
