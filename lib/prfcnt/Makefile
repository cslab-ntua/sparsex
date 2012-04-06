.PHONY: all clean

## Configuration
scripts      = $(shell pwd)/../../scripts

CC      = gcc
CFLAGS  = -Wall -O2
# run scripts to get information
ifneq ($(MAKECMDGOALS), clean)
	CPU ?= $(shell $(scripts)/cpu_info.sh)
endif

MHZ_SH  ?= $(scripts)/cpu_mhz.sh
override CPPFLAGS += -DCPU_$(CPU) -DCPU_MHZ_SH=\"$(MHZ_SH)\"
COMPILE  = $(CC) $(CPPFLAGS) $(CFLAGS)

ifeq ($(CPU),NIAGARA)
	CFLAGS += -m64
endif

# compile all .c as executables
progs := $(patsubst %.c,%,$(wildcard *.c))

all: $(progs)

% : %.c
	$(COMPILE) $(LDFLAGS) -o $* $<

%.s: %.c
	$(COMPILE) -S $<

%.i: %.c
	$(COMPILE) -E $< | indent -kr > $@

clean: 
	rm -rf $(progs) *.s *.i
