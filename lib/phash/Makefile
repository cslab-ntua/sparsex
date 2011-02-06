.PHONY: all clean

all: phash.o phash_i phash_test

GCC         ?= gcc
CFLAGS      ?= -Wall -Winline -O3 -Wdisabled-optimization -fPIC
CFLAGS      += -g
DEFS         = -DVAL_OVERLOAD
COMPILE      = $(GCC) $(CFLAGS) $(INC) $(DEFS)

phash.o: phash.c phash.h
	$(COMPILE) -c $< -o $@

phash_stats.o: phash.c phash.h
	$(COMPILE) -DPHASH_STATS -c $< -o $@

phash_i: phash.c
	$(COMPILE) -DPHASH_MAIN -DPHASH_STATS $< -o $@

phash_test: phash_test.o phash_stats.o
	$(COMPILE) $^ -o $@

%.s: %.c
	$(COMPILE) -S -fverbose-asm $<
%.o: %.c
	$(COMPILE) -c $<
%.i: %.c
	$(COMPILE) -E $< | indent -kr > $@

clean:
	rm -rf *.s *.o *.i phash_test phash_i
