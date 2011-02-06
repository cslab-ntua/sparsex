.PHONY: all clean

all: dynarray.o

common_dir = ../

## Configuration
GCC         ?= gcc
CFLAGS      ?= -Wall -Winline -O3 -Wdisabled-optimization -fPIC
CFLAGS      += -g
COMPILE      = $(GCC) $(CFLAGS) $(INC) $(DEFS)
INC          = -I$(common_dir)

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
