#
# Makefile.am -- Makefile for library benchmarks
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
	-I$(top_srcdir)/src/examples/examples \
	$(AX_CPPFLAGS)

AM_CFLAGS = $(AX_CFLAGS)
AM_CXXFLAGS = $(AX_CXXFLAGS)
AM_LDFLAGS	= \
	-L$(top_builddir)/src/internals/.libs \
	$(AX_LDFLAGS)

noinst_PROGRAMS = mmf_example csr_example reordering_example advanced_example matrix_caching_example_p1 matrix_caching_example_p2

mmf_example_SOURCES = \
	mmf_example.c

csr_example_SOURCES = \
	csr_example.c

reordering_example_SOURCES = \
	reordering_example.c

advanced_example_SOURCES = \
	advanced_example.c

matrix_caching_example_p1_SOURCES = \
	matrix_caching_example_p1.c

matrix_caching_example_p2_SOURCES = \
	matrix_caching_example_p2.c

mmf_example_LDADD = \
	$(top_builddir)/src/api/.libs/libsparsex.la \
	-ldl -lm

csr_example_LDADD = \
	$(top_builddir)/src/api/.libs/libsparsex.la \
	-ldl -lm

reordering_example_LDADD = \
	$(top_builddir)/src/api/.libs/libsparsex.la \
	-ldl -lm

advanced_example_LDADD = \
	$(top_builddir)/src/api/.libs/libsparsex.la \
	-ldl -lm

matrix_caching_example_p1_LDADD = \
	$(top_builddir)/src/api/.libs/libsparsex.la \
	-ldl -lm

matrix_caching_example_p2_LDADD = \
	$(top_builddir)/src/api/.libs/libsparsex.la \
	-ldl -lm
