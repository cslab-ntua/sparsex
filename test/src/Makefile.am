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

AM_CFLAGS = $(AX_CFLAGS)
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
	$(BOOST_THREAD_LIB) -ldl -lm

check_PROGRAMS = test_sparsex test_element test_dynarray test_timer \
	test_binary_p1 test_binary_p2

noinst_HEADERS = CsxCheck.hpp

test_sparsex_SOURCES = \
	CsxCheck.cpp \
	sparsex_test.c
test_sparsex_LDADD = \
	$(top_builddir)/src/api/.libs/libsparsex.la \
	$(SPARSEXIMPL_LDADD)

test_element_SOURCES = ElementTest.cpp
test_element_LDADD = $(SPARSEXIMPL_LDADD)

test_dynarray_SOURCES = DynamicArrayTest.cpp
test_dynarray_LDADD = $(SPARSEXIMPL_LDADD)

test_timer_SOURCES = TimerTest.cpp
test_timer_LDADD = $(SPARSEXIMPL_LDADD)

test_binary_p1_SOURCES = BinaryTest_p1.cpp
test_binary_p1_LDADD = $(SPARSEXIMPL_LDADD)

test_binary_p2_SOURCES = BinaryTest_p2.cpp
test_binary_p2_LDADD = $(SPARSEXIMPL_LDADD)
