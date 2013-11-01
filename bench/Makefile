.PHONY: all clean

all: bench

# ifndef interface
# interface=lp64
# endif

# ifndef MKLROOT
# MKLROOT = ../..
# endif
# MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
# CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

# CC=gcc
# IFACE_COMP_PART=intel
# IFACE_THREADING_PART=intel

log_dir		 = ../csx/logger
csx_dir     = ../csx
libcsx_dir	= ../C-API
mkl_dir		= /opt/intel/composer_xe_2013.5.192/mkl/include
POSKIINC	= /home/athena/Desktop/Products/pOSKI/build-1.0.0/include/poski
POSKILIB	= /home/athena/Desktop/Products/pOSKI/build-1.0.0/lib
OSKIINC		= /home/athena/Desktop/Products/pOSKI/build-1.0.0/build_oski/include
OSKILIB		= /home/athena/Desktop/Products/pOSKI/build-1.0.0/build_oski/lib/oski

CXX         = g++
LD          = g++
CXXFLAGS    = -g -Wall -O3 -Wdisabled-optimization -Winline -fPIC
CPPFLAGS    = -I$(libcsx_dir) -I$(csx_dir) -I$(log_dir) -I$(mkl_dir) 
#-I$(POSKIINC) -I$(OSKIINC)
MKL_CXXFLAGS	= -fopenmp -m64
POSKI_CXXFLAGS	= -fopenmp
LIBCSX_LDFLAGS  = -L$(libcsx_dir) -Wl,--rpath=$(libcsx_dir) -lcsx -L$(csx_dir) -Wl,--rpath=$(csx_dir) -lcsximpl
POSKI_LDFLAGS   = -Wl,--rpath -Wl,$(POSKILIB) -L$(POSKILIB) -lposki -Wl,--rpath -Wl,$(OSKILIB) -L$(OSKILIB)
MKL_LDFLAGS     = -Wl,--start-group /opt/intel/composer_xe_2013.5.192/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2013.5.192/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/composer_xe_2013.5.192/mkl/lib/intel64/libmkl_core.a -Wl,--end-group
LIBS        	= -lpthread -ldl -lm

COMPILE = $(CXX) $(CXXFLAGS) $(CPPFLAGS)
LINK    = $(LD) 

%.o: %.cc %.h
	$(COMPILE) -c $< -o $@

poski_module.o: poski_module.cc
	$(COMPILE) $(POSKI_CXXFLAGS) -c $< -o $@

mkl_module.o: mkl_module.cc
	$(COMPILE) $(MKL_CXXFLAGS) -c $< -o $@

# bench: main.o bench.o timer.o poski_module.o libcsx_module.o mkl_module.o
# 	$(LINK) -o $@ $^ $(LIBCSX_LDFLAGS) $(MKL_LDFLAGS) $(POSKI_LDFLAGS) $(LIBS)
bench: main.o bench.o timer.o libcsx_module.o mkl_module.o
	$(LINK) -o $@ $^ $(LIBCSX_LDFLAGS) $(MKL_LDFLAGS) $(MKL_CXXFLAGS) $(LIBS)

clean:
	rm -f *.o bench
