
SUBDIRS = lib/dynarray lib/prfcnt lib/spm csx

.PHONY: all clean $(SUBDIRS)

all clean: $(SUBDIRS)

$(SUBDIRS):
	 $(MAKE) -C $@ $(MAKECMDGOALS)
