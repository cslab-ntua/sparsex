dnl
dnl configure.ac -- System configuration checks
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl

AC_PREREQ([2.62])
AC_INIT([SparseX], [1.0.0], [], [sparsex])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([build_config.h])
AC_CONFIG_MACRO_DIR([m4check])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

dnl Initialize automake
AM_INIT_AUTOMAKE([1.11])
AM_SILENT_RULES

dnl Add options to configure script
AX_SELECT_BUILD

dnl Make `AC_PROG_CXX' macro think that `CXXFLAGS' is already defined by the
dnl user, so that it does not touch it.

if test -z "$CXXFLAGS"; then
    have_cxxflags="no"
    CXXFLAGS="-Ibogus"
else
    have_cxxflags="yes"
fi

if test -z "$CFLAGS"; then
    have_cflags="no"
    CFLAGS="-Ibogus"
else
    have_cflags="yes"
fi

AC_PROG_CC([gcc])
AC_PROG_CXX([g++])
if test x"$have_cxxflags" == x"no"; then
    # Restore `CXXFLAGS' to nothing.
    CXXFLAGS=""
else
    # Check user's flag for validity
    AX_CHECK_CXXFLAGS
fi

if test x"$have_cflags" == x"no"; then
    # Restore `CFLAGS' to nothing.
    CFLAGS=""
else
    # Check user's flag for validity
    AX_CHECK_CFLAGS
fi

# Check also user supplied CPPFLAGS for validity
AX_CHECK_CPPFLAGS

dnl Initialize libtool after compilers; its sets CFLAGS otherwise
LT_PREREQ([2.2.2])
LT_INIT

dnl --enable-bench option
AC_ARG_ENABLE([bench],
   [AS_HELP_STRING([--enable-bench@<:@={yes|no}@:>@],
       [Compile also library benchmarks. Default value is `no'.])]
   [AS_HELP_STRING([--disable-bench], [Disable library benchmarks.])],
   [
       if test x"$enableval" != x"yes" && test x"$enableval" != x"no"; then
           AC_MSG_ERROR([Invalid value for option `--enable-bench'. dnl
Type `configure --help' for more information.])
       fi
   ], [enable_bench=no])

dnl Check for programs
AC_PROG_INSTALL
AC_PROG_SED
AC_PROG_GREP
AX_PROG_DOXYGEN

dnl Check for libraries.
AX_CHECK_BOOST([1.48])
AX_CHECK_BOOST_LIB([regex])
AX_CHECK_BOOST_LIB([serialization])
AX_CHECK_BOOST_LIB([system])
AX_CHECK_BOOST_LIB([thread])
AX_CHECK_LLVM([3.0])
AC_CHECK_LIB([pthread], [pthread_create], [],
             [
                 AC_MSG_ERROR([Could not find a proper Pthread implementation.])
             ])

dnl Always check for libnuma, since it is needed for linking
AC_CHECK_LIB([numa], [numa_realloc], [],
             [
                 AC_MSG_ERROR([Could not find a proper libnuma installation dnl
Please make sure you have installed correctly libnuma >= 2.0.8.])
             ])

dnl If library benchmarks are enabled, check also for MKL and pOSKI
if test x"$enable_bench" = x"yes"; then
    SPX_BUILD_BENCH=1
else
    SPX_BUILD_BENCH=0
fi

dnl Other system checks
AX_SYSTEM_MEMORY_NODES
if test $ax_memory_nodes -eq 0; then
    AC_MSG_WARN([Could not determine the number of memory nodes of your dnl
system. Is sysfs properly mounted? Assuming a single-node system.])
    SPX_USE_NUMA=0
elif test $ax_memory_nodes -gt 1; then
    SPX_USE_NUMA=1
else
    SPX_USE_NUMA=0
fi

dnl Set necessary flags
AX_CXXFLAGS=""
AX_CPPFLAGS="$BOOST_CPPFLAGS $LLVM_CPPFLAGS"
AX_LDFLAGS="$BOOST_LDFLAGS $LLVM_LDFLAGS"

dnl Set up the compilation flags for the package.
dnl If user did not specify any custom flags, use our default
AX_PACKAGE_FLAGS

if test -z "$CPPFLAGS"; then
   AX_PKG_CPPFLAGS=$ax_pkg_cppflags
else
   AX_PKG_CPPFLAGS=""
fi

if test -z "$CXXFLAGS"; then
   AX_PKG_CXXFLAGS=$ax_pkg_cxxflags
else
   AX_PKG_CXXFLAGS=""
fi

if test -z "$CFLAGS"; then
   AX_PKG_CFLAGS=$ax_pkg_cflags
else
   AX_PKG_CFLAGS=""
fi

if test -z "$LDFLAGS"; then
   AX_PKG_LDFLAGS=$ax_pkg_ldflags
else
   AX_PKG_LDFLAGS=""
fi

AX_CPPFLAGS="$AX_CPPFLAGS $AX_PKG_CPPFLAGS"
AX_CXXFLAGS="$AX_CXXFLAGS $AX_PKG_CXXFLAGS"
AX_CFLAGS="$AX_CFLAGS $AX_PKG_CFLAGS"
AX_LDFLAGS="$AX_LDFLAGS $AX_PKG_LDFLAGS"

dnl Initialize variables for package
dnl AX_PACKAGE_INIT

SPX_MULT_TEMPLATE_DIR="`pwd`/$srcdir/src/templates"
SPX_JIT_INCLUDE="`pwd`/$srcdir/include"

AM_CONDITIONAL([USE_LIB_SUFFIX], [test ! -z $ax_lib_suffix])
AM_CONDITIONAL([SPX_BUILD_BENCH], [ test x"$enable_bench" = x"yes" ])
AC_SUBST([ax_lib_suffix])
AC_SUBST([AX_PKG_CPPFLAGS])
AC_SUBST([AX_PKG_CXXFLAGS])
AC_SUBST([AX_PKG_LDFLAGS])
AC_SUBST([AX_CPPFLAGS])
AC_SUBST([AX_CXXFLAGS])
AC_SUBST([AX_CFLAGS])
AC_SUBST([AX_LDFLAGS])
AC_SUBST([AX_TEST_LDFLAGS])
AC_SUBST([BOOST_REGEX_LIB])
AC_SUBST([BOOST_SERIALIZATION_LIB])
AC_SUBST([BOOST_SYSTEM_LIB])
AC_SUBST([BOOST_THREAD_LIB])
AC_SUBST([SPX_MULT_TEMPLATE_DIR])
AC_SUBST([SPX_JIT_INCLUDE])
AC_SUBST([SPX_USE_NUMA])

AC_CONFIG_FILES([
        Makefile
        include/Makefile
        include/sparsex/internals/Config.hpp
        src/Makefile
        src/api/Makefile
        src/bench/Makefile
        src/internals/Makefile
        src/templates/Makefile
        test/Makefile
        test/src/Makefile
        ])

dnl AC_CONFIG_FILES([scripts/tn-base-config],
dnl                 [chmod +x scripts/tn-base-config])

AC_OUTPUT

dnl Print a summary of the configuration
AX_SUMMARY_CONFIG
