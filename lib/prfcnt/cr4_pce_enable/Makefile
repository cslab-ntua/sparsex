ifneq ($(KERNELRELEASE),)

obj-m    := cr4_pce_enable.o

else

KDIR     := /lib/modules/$(shell uname -r)/build
PWD      := $(shell pwd)

default: 
	$(MAKE) -C $(KDIR) SUBDIRS=$(PWD) modules

clean:
	rm -rf *.o *.ko .cr4* *.mod.c *.tmp*
endif
