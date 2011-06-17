SUBDIRS = scripts lib/prfcnt lib/dynarray lib/spm csx

.PHONY: all clean $(SUBDIRS)

all clean: $(SUBDIRS)

$(SUBDIRS):
	 $(MAKE) -C $@ $(MAKECMDGOALS)
