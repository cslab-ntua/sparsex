.PHONY: all clean

all: libsparsex.so libsparsex.a mmf_example csr_example reordering_example

## Configuration
csx_dir     = ../csx
logger_dir	= ../csx/logger
script_dir	= ../scripts
CC          = gcc
LD          = gcc
CFLAGS      = -g -Wall -O3 -Wno-unused-function -fPIC
# -Wpadded
LDFLAGS     = -L$(csx_dir) -Wl,--rpath=$(csx_dir) -L$(logger_dir) -Wl,--rpath=$(logger_dir) -D_GNU_SOURCE
COMPILE     = $(CC) $(CPPFLAGS) $(CFLAGS)
LINK        = $(LD) $(LDFLAGS)
LIBS        = -lcsximpl

override CFLAGS += -I$(csx_dir) -I$(logger_dir)
override LDFLAGS  += -Wl,--allow-multiple-definition -Wl,--rpath=$(shell pwd)

ifeq ($(shell $(script_dir)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

LIB_OBJECTS = error.o common.o mat_vec.o

%.o: %.c %.h
	$(COMPILE) -D_GNU_SOURCE -c $< -o $@

libsparsex.so: $(LIB_OBJECTS)
	$(LINK) -shared -o $@ $^ $(LIBS)

libsparsex.a: $(LIB_OBJECTS)
	ar cru $@ $^

mmf_example: mmf_example.o libsparsex.so
	$(LINK) -o $@ $^ $(LIBS)

csr_example: csr_example.o libsparsex.so
	$(LINK) -o $@ $^ $(LIBS)

reordering_example: reordering_example.o libsparsex.so
	$(LINK) -o $@ $^ $(LIBS)

clean:
	rm -f *.o libsparsex.so libsparsex.a mmf_example csr_example reordering_example
