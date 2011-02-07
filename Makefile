
SUBDIRS = lib/dynarray lib/spm lib/cpu lib/prfcnt csx

.PHONY: all clean $(SUBDIRS)

all clean: $(SUBDIRS) 

$(SUBDIRS):
	 $(MAKE) -C $@ $(MAKECMDGOALS)
