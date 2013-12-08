.PHONY: all clean

all: bench

script_dir	= ../scripts
log_dir		= ../csx/logger
csx_dir     = ../csx
libcsx_dir	= ../api
mkl_dir		= /opt/intel/composer_xe_2013_sp1.1.106/mkl
POSKIINC	= /home/athena/Desktop/Products/pOSKI/build-1.0.0/include/poski
POSKILIB	= /home/athena/Desktop/Products/pOSKI/build-1.0.0/lib
OSKIINC		= /home/athena/Desktop/Products/pOSKI/build-1.0.0/build_oski/include
OSKILIB		= /home/athena/Desktop/Products/pOSKI/build-1.0.0/build_oski/lib/oski

CXX         = g++
LD          = g++
CXXFLAGS    = -g -Wall -O3 -Wdisabled-optimization -Winline -fPIC
MKL_CXXFLAGS	= -fopenmp -m64
POSKI_CXXFLAGS	= -fopenmp
LIBS        	= -fopenmp -lpthread -ldl -lm

LIBCSX_LDFLAGS  = -L$(libcsx_dir) -Wl,--rpath=$(libcsx_dir) -lcsx -L$(csx_dir) -Wl,--rpath=$(csx_dir) -lcsximpl
MKL_LDFLAGS     = -Wl,--start-group $(mkl_dir)/lib/intel64/libmkl_intel_lp64.a $(mkl_dir)/lib/intel64/libmkl_gnu_thread.a $(mkl_dir)/lib/intel64/libmkl_core.a -Wl,--end-group
POSKI_LDFLAGS   = -Wl,--rpath -Wl,$(POSKILIB) -L$(POSKILIB) -lposki -Wl,--rpath -Wl,$(OSKILIB) -L$(OSKILIB)

override CPPFLAGS += -I$(libcsx_dir) -I$(csx_dir) -I$(log_dir) -I$(mkl_dir)/include -I$(POSKIINC) -I$(OSKIINC)
override LDFLAGS  += $(LIBCSX_LDFLAGS) 
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

bench: main.o bench.o timer.o libcsx_module.o
	$(LINK) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o bench
