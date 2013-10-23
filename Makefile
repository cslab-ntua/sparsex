SUBDIRS = scripts lib/prfcnt lib/dynarray csx C-API
# bench lib/spm cg

.PHONY: all clean $(SUBDIRS)
.NOTPARALLEL: $(SUBDIRS)

all clean: $(SUBDIRS)

$(SUBDIRS):
	 $(MAKE) -C $@ $(MAKECMDGOALS)
