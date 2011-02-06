
SUBDIRS = lib/dynarray lib/phash lib/spm csx

all: $(SUBDIRS)

.PHONY: all $(SUBDIRS)

$(SUBDIRS):
	 $(MAKE) -C $@
