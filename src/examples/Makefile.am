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
	-I$(top_srcdir)/src/examples \
	$(AX_CPPFLAGS)

AM_CFLAGS = $(AX_CFLAGS)
AM_CXXFLAGS = $(AX_CXXFLAGS)
AM_LDFLAGS = \
	-L$(top_builddir)/src/api/.libs \
	$(AX_LDFLAGS)

SPARSEX_LDADD = \
	$(top_builddir)/src/api/.libs/libsparsex.la

noinst_PROGRAMS = mmf_example csr_example reordering_example advanced_example matrix_caching_example_p1 matrix_caching_example_p2

# Dummy C++ source to cause C++ linking. This is required becaused the
# sublibraries include C++ source code, so we must force automake to select
# the appropriate linker.
nodist_EXTRA_mmf_example_SOURCES = dummy.cxx
mmf_example_SOURCES = \
	mmf_example.c
mmf_example_LDFLAGS = -rpath $(LLVM_LIBDIR)
mmf_example_LDADD = $(SPARSEX_LDADD)

nodist_EXTRA_csr_example_SOURCES = dummy.cxx
csr_example_SOURCES = \
	csr_example.c
csr_example_LDFLAGS = -rpath $(LLVM_LIBDIR)
csr_example_LDADD = $(SPARSEX_LDADD)

nodist_EXTRA_reordering_example_SOURCES = dummy.cxx
reordering_example_SOURCES = \
	reordering_example.c
reordering_example_LDFLAGS = -rpath $(LLVM_LIBDIR)
reordering_example_LDADD = $(SPARSEX_LDADD)

nodist_EXTRA_advanced_example_SOURCES = dummy.cxx
advanced_example_SOURCES = \
	advanced_example.c
advanced_example_LDFLAGS = -rpath $(LLVM_LIBDIR)
advanced_example_LDADD = $(SPARSEX_LDADD)

nodist_EXTRA_matrix_caching_example_p1_SOURCES = dummy.cxx
matrix_caching_example_p1_SOURCES = \
	matrix_caching_example_p1.c
matrix_caching_example_p1_LDFLAGS = -rpath $(LLVM_LIBDIR)
matrix_caching_example_p1_LDADD = $(SPARSEX_LDADD)

nodist_EXTRA_matrix_caching_example_p2_SOURCES = dummy.cxx
matrix_caching_example_p2_SOURCES = \
	matrix_caching_example_p2.c
matrix_caching_example_p2_LDFLAGS = -rpath $(LLVM_LIBDIR)
matrix_caching_example_p2_LDADD = $(SPARSEX_LDADD)
