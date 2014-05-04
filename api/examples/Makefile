.PHONY: all clean

all: mmf_example csr_example reordering_example

SPARSEX_DIR = ..
CSX_DIR     = ../../csx
LOGGER_DIR	= ../../csx/logger
CC          = gcc
LD          = gcc
CFLAGS      = -g -Wall -O3 -fPIC
LDFLAGS     = 
COMPILE     = $(CC) $(CPPFLAGS) $(CFLAGS)
LINK        = $(LD) $(LDFLAGS)
LIBS        = -lsparsex -lcsximpl

override CFLAGS += -I$(SPARSEX_DIR) -I$(CSX_DIR) -I$(LOGGER_DIR)
override LDFLAGS += -L$(SPARSEX_DIR) -Wl,--rpath=$(SPARSEX_DIR) -L$(CSX_DIR) -Wl,--rpath=$(CSX_DIR)

mmf_example: mmf_example.o
	$(LINK) -o $@ $^ $(LIBS)

csr_example: csr_example.o
	$(LINK) -o $@ $^ $(LIBS)

reordering_example: reordering_example.o
	$(LINK) -o $@ $^ $(LIBS)

clean:
	rm -f *.o mmf_example csr_example reordering_example
