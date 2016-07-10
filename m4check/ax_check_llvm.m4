dnl
dnl ax_check_llvm.m4 -- Check for LLVM and Clang installation
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl


dnl
dnl AX_CHECK_STRICT_LLVM(required_version)
dnl
dnl     Check for a specific set of versions of LLVM
dnl
dnl     AC_SUBST(LLVM_CPPFLAGS)
dnl     AC_SUBST(LLVM_LDFLAGS)
dnl     AC_SUBST(LLVM_LIBS)
dnl

AC_DEFUN([AX_CHECK_STRICT_LLVM],
[
    AC_REQUIRE([AC_PROG_SED])
    AC_REQUIRE([AC_PROG_GREP])

    m4_if([$#], [0],
          [m4_errprintn([Too few arguments to $0.]) m4_exit(1)])

    llvm_supported_versions="$1"

    AC_MSG_CHECKING([for supported LLVM versions dnl
@<:@$llvm_supported_versions@:>@])
    AC_ARG_WITH([llvm],
        [AS_HELP_STRING([--with-llvm=CONFIG],
                [use CONFIG as LLVM configuration script.])],
        [], [with_llvm="llvm-config"])

    if test -z $withval; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([No LLVM config script given. Please supply one.])
    fi

    dnl Default LLVM components required by SparseX
    default_llvm_components="core analysis executionengine jit native dnl
nativecodegen codegen option bitreader ipo linker bitwriter asmparser dnl
instrumentation"

    AC_ARG_WITH([llvm-components],
                [AS_HELP_STRING([--with-llvm-components@<:@=COMPONENTS@:>@],
                                [additional LLVM components to link against.])],
                [], [with_llvm_components="$default_llvm_components"])

    llvm_config_prog=$with_llvm

    dnl Run a sanity check that the user has not given us an irrelevant program
    $llvm_config_prog --help 2>&1 | grep LLVM &> /dev/null && \
        seems_llvm_config=1

    if test -z $seems_llvm_config; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([`$llvm_config_prog' doesn't seem like an actual dnl
llvm-config script. Please check your supplied argument.])
    fi

    llvm_version=`$llvm_config_prog --version 2> /dev/null | \
$SED -e 's/@<:@^0-9@:>@*$//g'`
    if test -z $llvm_version; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Could not find a proper LLVM config script. dnl
Tried `$llvm_config_prog' and failed.])
    fi

    # Strictly check for the versions supplied
    llvm_found=0
    for version in $llvm_supported_versions; do
        if test $llvm_version == $version; then
            llvm_found=1
            break
        fi
    done

    if test $llvm_found -eq 0; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Could not find any of the supported LLVM versions dnl
($llvm_supported_versions). Please check your LLVM installation or dnl
try using the `--with-llvm' option.])
    else
        AC_MSG_RESULT([yes (found version $llvm_version)])
    fi

    dnl Everything's fine with LLVM; set the variables
    LLVM_CPPFLAGS=`$llvm_config_prog --cppflags`
    LLVM_LDFLAGS=`$llvm_config_prog --ldflags`
    LLVM_LIBS=`$llvm_config_prog --libs`
    #  $with_llvm_components 2> /dev/null

    dnl The zlib and tinfo (or ncurses) libraries need to be linked for llvm-3.5
    if test $llvm_version == "3.5.0"; then
        AC_CHECK_LIB([z], [compress2], [],
            [AC_MSG_ERROR([Could not find zlib library required by dnl
LLVM $llvm_version.])])
        crc32_found=1
        AC_CHECK_LIB([ncurses], [crc32], [], [crc32_found=0]) 
        if test $crc32_found -eq 0; then
            AC_CHECK_LIB([tinfo], [crc32], [],
            [AC_MSG_ERROR([Could not find tinfo or ncurses library required by dnl
LLVM $llvm_version.])])
        fi
    fi

    dnl Check for Clang
    CLANG_PREFIX=`$llvm_config_prog --prefix`
    AX_CHECK_PROG([clang], [CLANG], [$CLANG_PREFIX/bin])
    if test -z $CLANG; then
        AC_MSG_ERROR([Could not find the Clang compiler. dnl
Please check your installation of LLVM and Clang.])
    fi

    dnl Required Clang libs
    clang_libs="-lclangFrontendTool -lclangFrontend -lclangDriver dnl
    -lclangSerialization -lclangCodeGen -lclangParse -lclangSema"
    if test $llvm_version == "3.5.0"; then
        clang_libs="$clang_libs -lclangRewriteFrontend"
    fi
    clang_libs="$clang_libs -lclangRewrite -lclangStaticAnalyzerFrontend dnl
    -lclangStaticAnalyzerCheckers -lclangStaticAnalyzerCore -lclangAnalysis"
    if test $llvm_version == "3.5.0"; then
        clang_libs="$clang_libs -lclangEdit"
    fi
    clang_libs="$clang_libs -lclangAST -lclangLex -lclangBasic"
    CLANG_LIBS=$clang_libs

    # dnl Retrieve all Clang's internal libs
    # clang_libfiles=`$llvm_config_prog --libdir`/libclang*
    # clang_libs=""
    # for libf in $clang_libfiles; do
    #     soname=`basename $libf`
    #     libname=`echo $soname | sed -e 's/^lib//' | sed -e 's/\..*//'`
    #     clang_libs="$clang_libs -l$libname"
    # done
    
    # dnl Surround libraries with appropriate linker flags to allow random search
    # dnl FIXME: these flags work for gcc only
    # clang_libs="-Wl,--start-group $clang_libs -Wl,--end-group"
    # CLANG_LIBS="$clang_libs"

    dnl Retrieve Clang's system header search path
    clang_inc_search_dirs=`${CLANG_PREFIX}/bin/clang -E -v -xc -std=c99 - dnl
< /dev/null 2>&1 | grep '^@<:@@<:@:space:@:>@@:>@\/.*'`
    for d in $clang_inc_search_dirs; do
        clang_inc_search_path="$clang_inc_search_path:$d"
    done

    CLANG_INC_SEARCH_PATH=$clang_inc_search_path

    dnl Substitute the variables and return
    AC_SUBST([LLVM_CPPFLAGS])
    AC_SUBST([LLVM_LDFLAGS])
    AC_SUBST([LLVM_LIBS])
    AC_SUBST([CLANG_LIBS])
    AC_SUBST([CLANG_INC_SEARCH_PATH])

    dnl Add to LIBS
dnl    LIBS="$CLANG_LIBS $LLVM_LIBS $LIBS"
])
