SUBDIRS = scripts lib/prfcnt lib/dynarray lib/spm csx cg

.PHONY: all clean $(SUBDIRS)
.NOTPARALLEL: $(SUBDIRS)

all clean: $(SUBDIRS)

$(SUBDIRS):
	 $(MAKE) -C $@ $(MAKECMDGOALS)
