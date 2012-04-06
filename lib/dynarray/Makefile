.PHONY: all clean

all: dynarray.o

common_dir = ../

## Configuration
CC          ?= gcc
CFLAGS      ?= -g -Wall -O3 -Wdisabled-optimization -fPIC
override CPPFLAGS += -I$(common_dir)
COMPILE      = $(CC) $(CPPFLAGS) $(CFLAGS)

dynarray.o: dynarray.c dynarray.h
	$(COMPILE) -c $< -o $@

dynarray.llvm.bc: dynarray.c dynarray.h
	llvm-gcc -Wall $(INC) -emit-llvm -c $< 

%.s: %.c
	$(COMPILE) -S -fverbose-asm $<
%.o: %.c
	$(COMPILE) -c $<
%.i: %.c
	$(COMPILE) -E $< | indent -kr > $@

clean:
	rm -rf *.s *.o *.i
