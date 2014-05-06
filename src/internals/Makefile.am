#
# Makefile.am -- SparseX library
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
AM_LDFLAGS  = $(AX_LDFLAGS)

lib_LTLIBRARIES = libsparseximpl.la
libsparseximpl_la_SOURCES = \
	numa_util.c \
	Affinity.cpp \
	Barrier.cpp \
	Compiler.cpp \
	CsxBench.cpp \
	CsxBuild.cpp \
	CsxMatvec.cpp \
	CsxUtil.cpp \
	CtlBuilder.cpp \
	Encodings.cpp \
	Facade.cpp \
	Jit.cpp \
	JitUtil.cpp \
	Node.cpp \
	Runtime.cpp \
	Statistics.cpp \
	TemplateText.cpp \
	ThreadPool.cpp \
	Timer.cpp \
	Vector.cpp \
	logger/Logger.cpp \
	logger/OutputPolicy.cpp

libsparseximpl_la_LIBADD = \
	$(BOOST_REGEX_LIB) \
	$(BOOST_SERIALIZATION_LIB) \
	$(BOOST_SYSTEM_LIB) \
	$(BOOST_THREAD_LIB) \
	$(LLVM_LIBS)

if USE_LIB_SUFFIX
libsparseximpl_la_LDFLAGS = $(AM_LDFLAGS) -release $(ax_lib_suffix)

install-exec-hook:
	cd $(DESTDIR)$(libdir) && \
mv -f libsparseximpl.a libsparseximpl-$(ax_lib_suffix).a && \
$(LN_S) libsparseximpl-$(ax_lib_suffix).a libsparseximpl.a

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/libsparseximpl-$(ax_lib_suffix).a
else
libsparseximpl_la_LDFLAGS = $(AM_LDFLAGS)
endif
