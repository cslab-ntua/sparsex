#
# Makefile.am -- SparseX tests
#
# Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
# Copyright (C) 2014, Vasileios Karakasis
# All righrs reserved.
#
# This file is distributed under the BSD License. See LICENSE.txt for details.
#

AM_DEFAULT_VERBOSITY = 0
ACLOCAL_AMFLAGS = -I m4check

AM_CPPFLAGS = \
	-I$(top_builddir)/include \
	-I$(top_srcdir)/include \
	$(AX_CPPFLAGS)
AM_CXXFLAGS = $(AX_CXXFLAGS)
AM_LDFLAGS	= \
	-L$(top_builddir)/src/internals/.libs \
	$(AX_LDFLAGS)

SPARSEXIMPL_LDADD = \
	$(top_builddir)/src/internals/.libs/libsparseximpl.la \
	$(CLANG_LIBS) \
	$(LLVM_LIBS) \
	$(BOOST_REGEX_LIB) \
	$(BOOST_SERIALIZATION_LIB) \
	$(BOOST_SYSTEM_LIB) \
	$(BOOST_THREAD_LIB)

check_PROGRAMS = test_spmv test_sparsex test_element test_dynarray test_timer

test_sparsex_SOURCES = \
	CsxBench.cpp \
	sparsex_test.c
test_sparsex_LDADD = \
	$(SPARSEXIMPL_LDADD) \
	$(top_builddir)/src/api/.libs/libsparsex.la

test_spmv_SOURCES = \
	CsxBench.cpp \
	SpmvTest.cpp
test_spmv_LDADD = $(SPARSEXIMPL_LDADD)

test_element_SOURCES = ElementTest.cpp
test_element_LDADD = $(SPARSEXIMPL_LDADD)

test_dynarray_SOURCES = DynamicArrayTest.cpp
test_dynarray_LDADD = $(SPARSEXIMPL_LDADD)

test_timer_SOURCES = TimerTest.cpp
test_timer_LDADD = $(SPARSEXIMPL_LDADD)
