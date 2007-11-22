.PHONY: all
.PHONY: clean

CC=gcc

CFLAGS= -Wall -O2 -g
progs := $(patsubst %.c,%,$(wildcard *.c))

all: $(progs) 

% : %.c
	$(CC) $(CFLAGS) -o $* $< 

%.s: %.c
	$(CC) $(CFLAGS) $(INC) $(DEFS) -S $< 

%.i: %.c
	$(CC) $(CFLAGS) $(INC) $(DEFS) -E $< | indent -kr > $@

clean:
	rm -rf $(progs) 
