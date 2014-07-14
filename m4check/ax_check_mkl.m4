dnl
dnl ax_check_mkl.m4 -- Macros relatede to the MKL library installation
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl Copy
dnl

AC_DEFUN([AX_CHECK_MKL],
[
    if test -z $1; then
        mkl_required_version="0.0.0"
    else
        mkl_required_version="$1"
    fi

    system_prefixes="/usr/local /opt/local /usr /opt"
    AC_MSG_CHECKING([for Intel MKL library >= $mkl_required_version])
    AC_ARG_WITH([mkl],
        [AS_HELP_STRING([--with-mkl@<:@=DIR@:>@],
                        [use Intel MKL library installed in DIR.])],
        [], [with_mkl="check"])

    if test x"$with_mkl" = x"check"; then
        search_prefixes=$system_prefixes
    else
        search_prefixes=$with_mkl
    fi

    CPPFLAGS_save="$CPPFLAGS"
    CXXFLAGS_save="$CXXFLAGS"
    AC_LANG_PUSH([C])
    for p in $search_prefixes; do
        CPPFLAGS=""
        CXXFLAGS=""
        AC_LINK_IFELSE([
            AC_LANG_PROGRAM([
                #include <stdio.h>
                #include "$p/include/mkl.h"
                ],
                [
                    printf("%d,%d,%d\n", __INTEL_MKL__,
                           __INTEL_MKL_MINOR__, __INTEL_MKL_UPDATE__);
                ])],
            [
                mkl_found=1
                MKL_CPPFLAGS="-I$p/include"
                MKL_LDFLAGS="-L$p/lib/intel64"
                mkl_version=`./conftest$EXEEXT`
                AS_VERSION_COMPARE([$mkl_version],
                                   [$mkl_required_version],
                                   [mkl_found=0])
            ],
            [ mkl_found=0 ])

        if test $mkl_found -eq 1; then
            break;
        fi
    done

    AC_LANG_POP([C])
    if test $mkl_found -eq 0; then
        AC_MSG_RESULT([no])
    else
        AC_DEFINE([SPX_BENCH_MKL], [1], [Build SpMV benchmarks with MKL.])
        AC_MSG_RESULT([yes (found version $mkl_version)])
        MKL_LIBS="-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core"
    fi

    # restore compiler flags
    CPPFLAGS="$CPPFLAGS_save"
    CXXFLAGS="$CXXFLAGS_save"

    AC_SUBST([MKL_CPPFLAGS])
    AC_SUBST([MKL_LDFLAGS])
    AC_SUBST([MKL_LIBS])
    AM_CONDITIONAL([SPX_BUILD_BENCH_MKL], [test $mkl_found -eq 1])
])

AC_DEFUN([AX_CHECK_POSKI],
[
    AC_ARG_WITH([poski],
        [AS_HELP_STRING([--with-poski@<:@=DIR@:>@],
                        [use pOSKI library installed in DIR.])],
        [])
    
    LD_LIBRARY_PATH_save=$LD_LIBRARY_PATH
    LD_LIBRARY_PATH="$with_poski/lib:$LD_LIBRARY_PATH"
    AC_CHECK_LIB([poski], [poski_Init],
                 [poski_found=1], [poski_found=0])

    if test $poski_found -eq 1; then
        AC_MSG_RESULT([yes])
        AC_DEFINE([SPX_BENCH_POSKI], [1], [Build SpMV benchmarks with pOSKI.])
        POSKI_CPPFLAGS="-I$with_poski/include -I$with_poski/build_oski/include"
        POSKI_LDFLAGS="-L$with_poski/lib -L$with_poski/build_oski/lib/oski"
        POSKI_LIBS="-lposki"
    fi

    LD_LIBRARY_PATH=$LD_LIBRARY_PATH_save
    AC_SUBST([POSKI_CPPFLAGS])
    AC_SUBST([POSKI_LDFLAGS])
    AC_SUBST([POSKI_LIBS])
    AM_CONDITIONAL([SPX_BUILD_BENCH_POSKI], [test $poski_found -eq 1])
])
