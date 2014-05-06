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
dnl AX_CHECK_LLVM(required_version)
dnl
dnl     Check for a specific version of LLVM
dnl
dnl     AC_SUBST(LLVM_CPPFLAGS)
dnl     AC_SUBST(LLVM_LDFLAGS)
dnl     AC_SUBST(LLVM_LIBS)
dnl

AC_DEFUN([AX_CHECK_LLVM],
[
    AC_REQUIRE([AC_PROG_SED])
    AC_REQUIRE([AC_PROG_GREP])
    m4_if([$#], [0],
          [m4_errprintn([Too few arguments to $0.]) m4_exit(1)])

    llvm_required_version=$1
    AC_MSG_CHECKING([for LLVM >= $llvm_required_version])
    AC_ARG_WITH([llvm],
                [AS_HELP_STRING([--with-llvm@<:@=CONFIG@:>@],
                                [use CONFIG as LLVM configuration script.])],
                [], [with_llvm="llvm-config"])

    dnl Default LLVM libraries required by SparseX
    default_llvm_libs="core analysis executionengine jit native bitreader ipo \
linker bitwriter asmparser instrumentation"

    AC_ARG_WITH([llvm-libs],
                [AS_HELP_STRING([--with-llvm-libs@<:@=LIBS@:>@],
                                [additional LLVM libraries to link against.])],
                [], [with_llvm_libs="$default_llvm_libs"])

    llvm_config_prog=$with_llvm
    llvm_version=`$llvm_config_prog --version 2> /dev/null`
    if test -z $llvm_version; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Could not find a proper LLVM config script. dnl 
Tried `$llvm_config_prog' and failed.])
    fi

    llvm_found=1
    AS_VERSION_COMPARE([$llvm_version], [$llvm_required_version],
                       [llvm_found=0])

    if test $llvm_found -eq 0; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Could not find the required LLVM version dnl
($llvm_required_version). Please check your LLVM installation or dnl
try using the `--with-llvm' option.])
    else
        AC_MSG_RESULT([yes (found version $llvm_version)])
    fi

    dnl Everything's fine with LLVM; set the variables
    LLVM_CPPFLAGS=`$llvm_config_prog --cppflags`
    LLVM_LDFLAGS=`$llvm_config_prog --ldflags`
    LLVM_LIBS=`$llvm_config_prog --libs $with_llvm_libs`

    dnl Check for Clang
    AX_CHECK_PROG([clang])
    if test -z $CLANG; then
        AC_MSG_ERROR([Could not find the Clang compiler. dnl
Please check your installation of LLVM and Clang.])
    fi

    dnl Required Clang libs
    CLANG_LIBS="-lclangFrontend -lclangSerialization -lclangDriver -lclangCodeGen -lclangParse -lclangSema -lclangStaticAnalyzerFrontend -lclangStaticAnalyzerCheckers -lclangStaticAnalyzerCore -lclangAnalysis -lclangRewrite -lclangAST -lclangLex -lclangBasic"

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
    clang_inc_search_dirs=`clang -E -v -xc -std=c99 - < /dev/null 2>&1 | grep '^@<:@@<:@:space:@:>@@:>@\/.*'`
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
