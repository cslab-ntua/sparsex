.phony: all clean

## Configuration
cpu_dir = $(shell pwd)/../cpu

CC      = gcc
CFLAGS  = -Wall -Winline -O2
CPU    ?= $(shell $(cpu_dir)/cpu_info.sh)
MHZ_SH ?= $(cpu_dir)/cpu_mhz.sh
DEFS   += -DCPU_$(CPU) -DCPU_MHZ_SH=\"$(MHZ_SH)\"

ifeq ($(CPU),NIAGARA)
	CFLAGS += -m64
endif

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
