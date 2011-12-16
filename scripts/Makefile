.PHONY: all
.PHONY: clean

CC = gcc

CFLAGS = -Wall -O2 -g
progs := $(patsubst %.c,%,$(wildcard *.c))

all: $(progs) 

% : %.c
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $* $< 

%.s: %.c
	$(CC) $(CPPFLAGS) -S $< 

%.i: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -E $< | indent -kr > $@

clean:
	rm -rf $(progs) 
