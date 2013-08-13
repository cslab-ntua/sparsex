.PHONY: all clean

all: libtest

## Configuration
csx_dir     = ../csx
spm_dir     = ../lib/spm
prfcnt_dir   = ../lib/prfcnt

CC          = gcc
LD          = gcc
CFLAGS      = -g -Wall -O3 -Wdisabled-optimization -Winline -fPIC
# -Wpadded
LDFLAGS     = -L$(csx_dir) -Wl,--rpath=$(csx_dir)
COMPILE     = $(CC) $(CFLAGS)
LINK        = $(LD) $(LDFLAGS)
LIBS        = -lcsx

override CFLAGS += -I$(csx_dir) -I$(spm_dir)

LIB_OBJECTS = error.o common.o mat.o vec.o

%.o: %.c %.h
	$(COMPILE) -c $< -o $@

mat.o: $(spm_dir)/vector.h $(spm_dir)/spmv_method.h $(spm_dir)/spm_mt.h $(spm_dir)/spmv_loops_mt.h $(prfcnt_dir)/tsc.h $(csx_dir)/SparseMatrixWrapper.h

libtest: main.o $(LIB_OBJECTS)
	$(LINK) -o $@ $^ $(LIBS)

clean:
	rm -f *.o libtest
