#
# Makefile.am -- SparseX C api
#
# Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
# Copyright (C) 2014, Vasileios Karakasis
# All rights reserved.
#
# This file is distributed under the BSD License. See LICENSE.txt for details.
#

AM_DEFAULT_VERBOSITY = 0
ACLOCAL_AMFLAGS = -I m4check

AM_CPPFLAGS = \
	-I$(top_builddir)/include \
	-I$(top_srcdir)/include \
	$(AX_CPPFLAGS)

AM_CFLAGS = $(AX_CFLAGS)
AM_CXXFLAGS = $(AX_CXXFLAGS)
AM_LDFLAGS	= \
	-L$(top_builddir)/src/internals/.libs \
	$(AX_LDFLAGS)

lib_LTLIBRARIES = libsparsex.la
libsparsex_la_SOURCES = \
	common.c \
	error.c \
	matvec.c

libsparsex_la_LDFLAGS = $(AM_LDFLAGS) -version-info 1:0:0
libsparsex_la_LIBADD = -lsparseximpl
