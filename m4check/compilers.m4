dnl
dnl compilers.m4 -- Macros for setting up compilers
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl

dnl
dnl AX_CHECK_{CPPFLAGS,CFLAGS,CXXFLAGS}([flags])
dnl
dnl     Checks if `flags' or the user supplied CPP/CXX/CFLAGS (if `flags' is
dnl     empty) are supported by the system's compiler.
dnl

AC_DEFUN([AX_CHECK_CPPFLAGS],
[
# save user's flags first
CPPFLAGS_save=$CPPFLAGS
if test x"$1" != x; then
    CPPFLAGS="$1"
fi

AC_MSG_CHECKING([if preprocessor supports `$CPPFLAGS' option(s)])
AC_PREPROC_IFELSE(
    [AC_LANG_PROGRAM([], [return 0;])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_FAILURE([`$CPPFLAGS' option(s) are not supported by the dnl
preprocessor])])
    
CPPFLAGS="$CPPFLAGS_save"
])

AC_DEFUN([AX_CHECK_CFLAGS],
[
# save user's flags first
CPPFLAGS_save="$CPPFLAGS"
CFLAGS_save="$CFLAGS"

CPPFLAGS=""
if test x"$1" != x; then
    CFLAGS="$1"
fi

AC_MSG_CHECKING([if C compiler supports `$CFLAGS' option(s)])
AC_LANG_PUSH([C])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([], [return 0;])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_FAILURE([`$CFLAGS' option(s) are not supported by the C compiler])])
    
AC_LANG_POP

# restore user's flags now
CPPFLAGS="$CPPFLAGS_save"
CFLAGS="$CFLAGS_save"
])

AC_DEFUN([AX_CHECK_CXXFLAGS],
[
# save user's flags first
CPPFLAGS_save="$CPPFLAGS"
CXXFLAGS_save="$CXXFLAGS"

CPPFLAGS=""
if test x"$1" != x; then
    CXXFLAGS="$1"
fi

AC_MSG_CHECKING([if C++ compiler supports `$CXXFLAGS' option(s)])
AC_LANG_PUSH([C++])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([], [return 0;])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_FAILURE([`$CXXFLAGS' option(s) are not supported by the dnl
C++ compiler])])
    
AC_LANG_POP

# restore user's flags now
CPPFLAGS="$CPPFLAGS_save"
CXXFLAGS="$CXXFLAGS_save"
])

dnl
dnl AX_PACKAGE_FLAGS
dnl
dnl     Defines the default CPPFLAGS, CXXFLAGS and LDFLAGS for configured
dnl		packages and stores them in ax_pkg_cppflags, ax_pkg_cxxflags and
dnl		ax_pkg_ldflags, respectively.
dnl

AC_DEFUN([AX_PACKAGE_FLAGS],
[
AC_REQUIRE([AX_PLATFORM])
AC_REQUIRE([AX_SELECT_BUILD])

# Check compiler flags
AX_CHECK_CPPFLAGS([-DFOO])
AX_CHECK_CPPFLAGS([-UFOO])
AX_CHECK_CFLAGS([-g])
AX_CHECK_CFLAGS([-O0])
AX_CHECK_CFLAGS([-O3])
AX_CHECK_CFLAGS([-std=c99])
AX_CHECK_CFLAGS([-Wall])
AX_CHECK_CXXFLAGS([-std=c++0x])
AX_CHECK_CXXFLAGS([-Wall])
AX_CHECK_CXXFLAGS([-Woverloaded-virtual])

case "$ax_build_mode" in
    "debug")
        # undefine NDEBUG; LLVM doesn't want it
        ax_pkg_cppflags="-UNDEBUG -DSPX_DEBUG=1"
        ax_pkg_cxxflags="-g -O0 -std=c++0x -pedantic -Wall -Woverloaded-virtual"
        ax_pkg_cflags="-g -O0 -std=c99 -pedantic -Wall" ;;
    "release")
        ax_pkg_cppflags="-DNDEBUG"
        ax_pkg_cxxflags="-O3 -std=c++0x -pedantic -Wall"
        ax_pkg_cflags="-O3 -std=c99 -pedantic -Wall" ;;
    *)
        AC_MSG_ERROR(
            [@<:@BUG@:>@ should not have entered here: m4_location]) ;;
esac
])
