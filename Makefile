SUBDIRS = lib/dynarray lib/spm csx cg

.PHONY: all clean $(SUBDIRS)

all clean: $(SUBDIRS)

$(SUBDIRS):
	 $(MAKE) -C $@ $(MAKECMDGOALS)
