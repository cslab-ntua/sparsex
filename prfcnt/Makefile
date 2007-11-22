.phony: all clean

CC      = gcc
CFLAGS  = -Wall -Winline -O2
CPU    ?= $(shell ../cpu/cpu_info.sh)
MHZ    ?= $(shell ../cpu/cpu_mhz.sh)
DEFS   += -DCPU_$(CPU) -DCPU_MHZ=$(MHZ)

cc_def  = $(CC) $(CFLAGS) $(INC) $(DEFS)

progs := $(patsubst %.c,%,$(wildcard *.c))

all: $(progs)

% : %.c
	$(cc_def) -o $* $< 

%.s: %.c
	$(cc_def) -S $<

%.i: %.c
	$(cc_def) -E $< | indent -kr > $@

clean: 
	rm -rf $(progs) *.s *.i
