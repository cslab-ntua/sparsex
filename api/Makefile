.PHONY: all clean

all: libcsx.so libcsx.a libtest

## Configuration
csx_dir     = ../csx
logger_dir	= ../csx/logger
script_dir	= ../scripts
CC          = gcc
LD          = gcc
CFLAGS      = -g -Wall -O3 -Wdisabled-optimization -Winline -fPIC
# -Wpadded
LDFLAGS     = -L$(csx_dir) -Wl,--rpath=$(csx_dir) -L$(logger_dir) -Wl,--rpath=$(logger_dir)
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
	$(COMPILE) -c $< -o $@

libcsx.so: $(LIB_OBJECTS)
	$(LINK) -shared -o $@ $^ $(LIBS)

libcsx.a: $(LIB_OBJECTS)
	ar cru $@ $^

libtest: main.o libcsx.so
	$(LINK) -o $@ $^ $(LIBS)

clean:
	rm -f *.o libcsx.so libcsx.a libtest
