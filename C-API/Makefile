.PHONY: all clean

all: libcsx.so libcsx.a libtest

## Configuration
csx_dir     = ../csx

CC          = gcc
LD          = gcc
CFLAGS      = -g -Wall -O3 -Wdisabled-optimization -Winline -fPIC -DLIBCSX_LOG_FILE
# -Wpadded
LDFLAGS     = -L$(csx_dir) -Wl,--rpath=$(csx_dir)
COMPILE     = $(CC) $(CFLAGS)
LINK        = $(LD) $(LDFLAGS)
LIBS        = -lcsx_in

override CFLAGS += -I$(csx_dir)

LIB_OBJECTS = error.o common.o mat.o

%.o: %.c %.h
	$(COMPILE) -c $< -o $@

mat.o: $(csx_dir)/SparseMatrixWrapper.h

libcsx.so: $(LIB_OBJECTS)
	$(LINK) -shared -o $@ $^ $(LIBS)

libcsx.a: $(LIB_OBJECTS)
	ar cru $@ $^

libtest: main.o libcsx.a
	$(LINK) -o $@ $^ $(LIBS)

clean:
	rm -f *.o libcsx.so libcsx.a libtest
