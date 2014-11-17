dnl
dnl compilers.m4 -- Macros for setting up compilers
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl

dnl AX_CHECK_{CC,CXX}
dnl
dnl     Checks the C/C++ compiler provided by the user and returns a
dnl     platform-independent string describing the compiler in `ax_cc' or
dnl     `ax_cxx' respectively

AC_DEFUN([AX_CHECK_CC],
[
CPPFLAGS_save="$CPPFLAGS"
CFLAGS_save="$CFLAGS"
LDFLAGS_save="$LDFLAGS"
AC_LANG_PUSH([C])
AC_LINK_IFELSE([
    AC_LANG_PROGRAM([
        #include <stdio.h>
        ],
        [
            #if defined(__ICC)
            printf("%s\n", "icc");
            #elif defined(__GNUC__)
            printf("%s\n", "gcc");
            #else
            printf("%s\n", "unknown");
            #endif
        ])],
    [
        ax_cc=`./conftest$EXEEXT`
    ],
    [
        AC_MSG_ERROR([[AX_CHECK_CC] failed: This is probably a bug! dnl
Please report this to $PACKAGE_BUGREPORT.])
    ])

AC_LANG_POP([C])
CPPFLAGS="$CPPFLAGS_save"
CFLAGS="$CFLAGS_save"
LDFLAGS="$LDFLAGS_save"
])

AC_DEFUN([AX_CHECK_CXX],
[
CPPFLAGS_save="$CPPFLAGS"
CXXFLAGS_save="$CXXFLAGS"
LDFLAGS_save="$LDFLAGS"
AC_LANG_PUSH([C++])
AC_LINK_IFELSE([
    AC_LANG_PROGRAM([
        #include <stdio.h>
        ],
        [
            #if defined(__ICC)
            printf("%s\n", "icpc");
            #elif defined(__GNUC__)
            printf("%s\n", "g++");
            #else
            printf("%s\n", "unknown");
            #endif
        ])],
    [
        ax_cxx=`./conftest$EXEEXT`
    ],
    [
        AC_MSG_ERROR([[AX_CHECK_CXX] failed: This is probably a bug! dnl
Please report this to $PACKAGE_BUGREPORT.])
    ])

AC_LANG_POP([C++])
CPPFLAGS="$CPPFLAGS_save"
CXXFLAGS="$CXXFLAGS_save"
LDFLAGS="$LDFLAGS_save"
])


dnl
dnl AX_CXX_CHECK_CXX11_FEATURES
dnl
dnl     Checks if the supplied C++ compiler supports the C++11 features
dnl     used by SparseX
dnl

AC_DEFUN([AX_CXX_CHECK_CXX11_FEATURES],
[
AC_MSG_CHECKING([whether C++ compiler supports the required C++11 features])

CPPFLAGS_save="$CPPFLAGS"
CXXFLAGS_save="$CXXFLAGS"
LDFLAGS_save="$LDFLAGS"

CPPFLAGS=""
CXXFLAGS="-std=c++0x -pedantic"
LDFLAGS=""

AC_LANG_PUSH([C++])
AC_LINK_IFELSE([
    AC_LANG_PROGRAM([
        #include <atomic>
        #include <vector>

        using namespace std;

        class C
        {
        public:
            C(int v) { v_ = v; }
                C(const C &other) { }
                C(C &&other) { }
        private:
            int v_;
        };
    ],
    [
        vector<C> v;
        v.emplace_back(1);

        C *vd = v.data();
        atomic<int> a;
        a = 0;
        atomic_fetch_sub(&a, 1);
        a.store(10);
    ])],
    [
        AC_MSG_RESULT([yes])
    ],
    [
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([dnl
The supplied compiler does not support the required C++11 features dnl
for building $PACKAGE_NAME. Consider using a more recent version.])
    ])

AC_LANG_POP([C++])

CPPFLAGS="$CPPFLAGS_save"
CXXFLAGS="$CXXFLAGS_save"
LDFLAGS="$LDFLAGS_save"
])


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
        ax_pkg_cxxflags="-g -O0 -std=c++0x -pedantic -fopenmp -Wall -Woverloaded-virtual"
        ax_pkg_cflags="-g -O0 -std=c99 -pedantic -fopenmp -Wall" ;;
    "release")
        ax_pkg_cppflags="-DNDEBUG"
        ax_pkg_cxxflags="-O3 -std=c++0x -pedantic -fopenmp -Wall"
        ax_pkg_cflags="-O3 -std=c99 -pedantic -fopenmp -Wall" ;;
    *)
        AC_MSG_ERROR(
            [@<:@BUG@:>@ should not have entered here: m4_location]) ;;
esac
])
