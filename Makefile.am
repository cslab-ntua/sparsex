#
# Makefile.am -- Top level Makefile.am
#
# Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
# Copyright (C) 2014, Vasileios Karakasis
# All rights reserved.
#
# This file is distributed under the BSD License. See LICENSE.txt for details.
#

AM_DEFAULT_VERBOSITY = 0
ACLOCAL_AMFLAGS = -I m4check
SUBDIRS = include scripts src test

if SPX_GENERATE_DOC
SUBDIRS += doc
endif

silent_flag   = $(silent_flag_$(V))
silent_flag_  = $(silent_flag_$(AM_DEFAULT_VERBOSITY))
silent_flag_0 = --silent
AM_MAKEFLAGS  = $(silent_flag)

dist_noinst_DATA = LICENSE.txt doc/SparseX_UG.pdf

# Recursively invoke our new and/or redefined targets
doc doc-api doc-devel clean-doc mostlyclean-local:
	$(MAKE) $(AM_MAKEFLAGS) -C @dx_docdir@ $@

# Always install user guide, whether doxygen documentation is built or not
install-data-local:
	@echo -n "Installing user guide... "
	@(MKDIR_P) @docdir@
	@INSTALL_DATA@ $(top_srcdir)/doc/SparseX_UG.pdf @docdir@ && \
echo "done" || echo "failed"

uninstall-local:
	$(RM) @docdir@/SparseX_UG.pdf

.PHONY: doc doc-api doc-devel
