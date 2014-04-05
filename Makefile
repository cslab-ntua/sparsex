SUBDIRS = scripts csx api bench
# lib/spm cg

.PHONY: all clean $(SUBDIRS)
.NOTPARALLEL: $(SUBDIRS)

all clean: $(SUBDIRS)

$(SUBDIRS):
	 $(MAKE) -C $@ $(MAKECMDGOALS)
