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

AM_CPPFLAGS = \
	-I$(top_builddir)/include \
	-I$(top_srcdir)/include \
	$(AX_CPPFLAGS)
AM_CXXFLAGS = $(AX_CXXFLAGS)
AM_LDFLAGS	= \
	-L$(top_builddir)/src/internals/.libs \
	$(AX_LDFLAGS)

check_PROGRAMS = \
	test_spmv

test_spmv_SOURCES = SpmvTest.cpp
test_spmv_LDADD = \
	$(top_builddir)/src/internals/.libs/libsparseximpl.la \
	$(CLANG_LIBS) \
	$(LLVM_LIBS) \
	$(BOOST_REGEX_LIB) \
	$(BOOST_SERIALIZATION_LIB) \
	$(BOOST_SYSTEM_LIB) \
	$(BOOST_THREAD_LIB)
