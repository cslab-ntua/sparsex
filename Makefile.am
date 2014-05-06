#
# Makefile.am -- Top level Makefile.am
#
# Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
# Copyright (C) 2014, Vasileios Karakasis
# All righrs reserved.
#
# This file is distributed under the BSD License. See LICENSE.txt for details.
#

AM_DEFAULT_VERBOSITY = 0
ACLOCAL_AMFLAGS = -I m4check
SUBDIRS = include src test

silent_flag   = $(silent_flag_$(V))
silent_flag_  = $(silent_flag_$(AM_DEFAULT_VERBOSITY))
silent_flag_0 = --silent
AM_MAKEFLAGS  = $(silent_flag)

# Recursively invoke our new and/or redefined targets
doc doc-api doc-devel clean-doc mostlyclean-local:
	$(MAKE) $(AM_MAKEFLAGS) -C @dx_docdir@ $@

.PHONY: doc doc-api doc-devel
