.PHONY: all clean

all: libsparsex.so libsparsex.a

## Configuration
SRC_DIR		= .
CSX_DIR     = ../csx
LOGGER_DIR	= ../csx/logger
SCRIPT_DIR	= ../scripts

CC          = gcc
LD          = gcc
CFLAGS      = -g -Wall -O3 -fPIC
LDFLAGS     = -L$(CSX_DIR) -Wl,--rpath=$(CSX_DIR) -L$(LOGGER_DIR) -Wl,--rpath=$(LOGGER_DIR)
COMPILE     = $(CC) $(CPPFLAGS) $(CFLAGS)
LINK        = $(LD) $(LDFLAGS)
LIBS        = -lcsximpl

override CFLAGS += -I$(CSX_DIR) -I$(LOGGER_DIR)
override LDFLAGS  += -Wl,--allow-multiple-definition -Wl,--rpath=$(shell pwd)

ifeq ($(shell $(SCRIPT_DIR)/numa_lib.sh FOO), FOO)
	override CPPFLAGS += -DSPM_NUMA
endif

SPARSEX_SOURCES := $(wildcard $(SRC_DIR)/*.c)
SPARSEX_OBJECTS := $(SPARSEX_SOURCES:%.c=%.o)

$(SPARSEX_OBJECTS): %.o : %.c
	$(COMPILE) -c $< -o $@

libsparsex.so: $(SPARSEX_OBJECTS)
	$(LINK) -shared -o $@ $^ $(LIBS)

libsparsex.a: $(SPARSEX_OBJECTS)
	ar cru $@ $^

clean:
	rm -f *.o libsparsex.so libsparsex.a
