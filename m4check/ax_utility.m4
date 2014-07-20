dnl
dnl ax_utility.m4 -- Miscellaneous utility macros.
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl 

dnl
dnl UPPERCASE(text)
dnl
dnl     Convert `text' to uppercase and convert each `-' to `_'.
dnl
m4_define([UPPERCASE],
          [m4_bpatsubst(m4_translit([$1], [[a-z]], [[A-Z]]), [[-\.]], [_])])

dnl
dnl MAKE_VARIABLE(base_part, suffix)
dnl
dnl     Computes the name of a variable.
dnl
dnl     This macro computes the name of variable by concatenating `base_part',
dnl     an `_', and `suffix'. Before concatenating, `base_part' is converted to
dnl     uppercase and all `.' and `-' characters are converted to `_'.
dnl
dnl     Note: This macro is not intended to be public. It is advised to
dnl           refrain from using it directly in `configure.ac'.
dnl
dnl m4_define([MAKE_VARIABLE],
dnl           [m4_bpatsubst(m4_translit([$1], [[a-z]], [[A-Z]]), [[-\.]], [_])_$2])
m4_define([MAKE_VARIABLE], [m4_join([_], UPPERCASE($1), $2)])


dnl
dnl PACKAGE_VERSION_NUMBER(version_string)
dnl
dnl     Computes the package version number from `version_string'.
dnl
dnl     This macro computes the version number of the package being
dnl     configured given a version string in dotted notation, i.e.,
dnl     <major>.<minor>.<patch_level>. This package sets, but does *NOT*
dnl     substitute, the following variables:
dnl
dnl     AX_MAJOR_VERSION : the major version of the package being configured
dnl     AX_MINOR_VERSION : the minor version of the package being configured
dnl     AX_PATCH_LEVEL   : the patch level of the package being configured.
dnl     AX_VERSION_NUMBER: the computed version number of the package being
dnl                        configured.
dnl
dnl     Note: This macro is not intended to be public. It is advised to
dnl           refrain from using it directly in `configure.ac'. Use
dnl           `AX_PACKAGE_INIT' instead.
dnl
AC_DEFUN([PACKAGE_VERSION_NUMBER],
[
m4_if([$#], [0],[
            m4_errprintn([Too few arguments to PACKAGE_VERSION_NUMBER.])
            m4_exit(1)])

ax_version=$1

AX_MAJOR_VERSION=`echo $ax_version | \
    awk '{ split($[]0, a, "."); print a[[1]] }'`

AX_MINOR_VERSION=`echo $ax_version | \
    awk '{ split($[]0, a, "."); print a[[2]] }'`
AX_PATCH_LEVEL=`echo $ax_version |
    awk '{ split($[]0, a, "."); print a[[3]] }'`

let "AX_VERSION_NUMBER = 100000*$AX_MAJOR_VERSION +
                         100*$AX_MINOR_VERSION + $AX_PATCH_LEVEL"
])

dnl
dnl AX_CHECK_PRAGMA_MACRO
dnl
AC_DEFUN([AX_CHECK_PRAGMA_MACRO],
[
AC_MSG_CHECKING([if _Pragma macro is supported])
AC_LANG_PUSH([C++])

dnl Save CXXFLAGS
CXXFLAGS_save=$CXXFLAGS
CXXFLAGS=""

AC_COMPILE_IFELSE(
        AC_LANG_PROGRAM([
                #define AX_PRAGMA _Pragma("dummy")
                AX_PRAGMA
                ], []),
        [AX_PRAGMA_MACRO=1
        AC_DEFINE([HAVE_PRAGMA_MACRO], [1], [_Pragma macro is supported])
        AC_MSG_RESULT([yes])
        ],
        [AX_PRAGMA_MACRO=0
        AC_MSG_RESULT([no])
        AC_MSG_WARN([Compiler does not support the _Pragma macro])])

dnl Restore CXXFLAGS
CXXFLAGS=$CXXFLAGS_save
AC_LANG_POP
])

dnl
dnl AX_PACKAGE_INIT
dnl
dnl     Initialization of some package variables (version numbers, system date)
dnl
dnl     This macro sets and substitutes the following variables:
dnl
dnl     AX_PACKAGE_VERSION_NUMBER: the version number of the package being
dnl                                configured
dnl
AC_DEFUN([AX_PACKAGE_INIT],
[
    AC_REQUIRE([AC_PROG_SED])
    AC_REQUIRE([AX_SYSTEM_DATE])

    # strip the revision number from the package version
    pkg_version_stripped=`echo $PACKAGE_VERSION | $SED -e 's/_.*//'`

    PACKAGE_VERSION_NUMBER([$pkg_version_stripped])
    AX_PACKAGE_VERSION_NUMBER="$AX_VERSION_NUMBER"

    AC_SUBST([AX_PACKAGE_VERSION_NUMBER])
])

dnl
dnl AX_OPTION_REVNUM
dnl
dnl Adds `--enable-revnum' option to the configure script.
dnl
dnl This option, if enabled, adds the suffix `_<revnum>' to the package version
dnl string, where `revnum' is the latest svn revision number of the package.
dnl
AC_DEFUN([AX_OPTION_REVNUM],
[
    AC_REQUIRE([AC_PROG_GREP])
    AC_REQUIRE([AC_PROG_SED])
    AC_ARG_ENABLE(
        [revnum],
        [AS_HELP_STRING([--enable-revnum@<:@={yes|no}@:>@],
                        [Adds revision number in package version.
Default value is `no'.])
AS_HELP_STRING([--disable-revnum],
		       [Disables revision number in package version.])],
        [
            if test x"$enableval" != x"yes" && test x"$enableval" != x"no"; then
                AC_MSG_ERROR([Invalid value for option `--enable-revnum'. dnl
Type `configure --help' for more information.])
            fi
        ], [enable_revnum=no])
    
    if test x"$enable_revnum" = x"yes"; then
        AX_CHECK_PROG([svn])
        if test ! -z $SVN; then
            revnum=`$SVN info $srcdir | $GREP 'Revision:' | $SED -e 's/Revision: //'`
        else
            AC_MSG_ERROR([Could not find `svn' in PATH. dnl
Check your installation of `Subversion' and retry or don't use revision dnl
numbers at all.])
        fi

        if test ! -z $revnum; then
            PACKAGE_VERSION=${PACKAGE_VERSION}_$revnum
        else
            AC_MSG_WARN([Could not retrieve the svn revision number. dnl
Check if you have checked out correctly the repository. Revision numbers dnl
will not be used.])
        fi
    fi
])

dnl
dnl AX_SELECT_BUILD
dnl
dnl		Add support for selecting the build type (debug, release, etc.) 
dnl
dnl This macro defines the command-line variable `build_type' to let the user
dnl decide on the preferred build type of the project. Supported build types
dnl are:
dnl
dnl	 * debug-static: No optimization, debug symbols, static linking
dnl  * debug-shared: No optimization, debug symbols, dynamic linking
dnl	 * release-static: Full optimization, no debug symbols, static linking
dnl  * release-shared: Full optimization, no debug symbols, dynamic linking
dnl
dnl This macro simply defines three shell variables, depending on the build and
dnl link modes selected and whether a suffix will be appended to the resulting
dnl library names. The necessary compiler flags for implementing the available
dnl build types are not set by this macro (see AX_PACKAGE_FLAGS @
dnl compilers.m4). The defined variables are the following:
dnl
dnl	Variable Name		Values
dnl --------------------------
dnl ax_build_mode		debug | release
dnl ax_link_mode		static | shared
dnl ax_lib_suffix		<string> | null (i.e., undefined)
dnl
dnl These variables are *not* AC_SUBST'ed.
dnl
dnl The `ax_lib_suffix' must be appended to the library name, if defined.
dnl

AC_DEFUN([AX_SELECT_BUILD],
[
    AC_ARG_VAR([build_type], [@<:@release|debug@:>@ Select build type for SparseX. Default is `release'])
	if test -z $build_type; then
	    build_type="release"
	fi

	case "$build_type" in
		"debug-static")
		    ax_build_mode="debug"
			ax_link_mode="static"
			ax_lib_suffix="gd" ;;
		"debug" | "debug-shared")
		    ax_build_mode="debug"
			ax_link_mode="shared"
			ax_lib_suffix="gd" ;;
		"release-static")
		    ax_build_mode="release"
			ax_link_mode="static" ;;
		"release" | "release-shared")
		    ax_build_mode="release"
			ax_link_mode="shared" ;;
		*)
			AC_MSG_ERROR([Unknown build type `$build_type']) ;;
	esac

	AC_MSG_NOTICE([selected build: $build_type])
])


dnl
dnl AX_VERSION_COMPARE(version-1, version-2,
dnl                    [action-if-less], [action-if-equal],
dnl                    [action-if-greater])
dnl
dnl     Compare two version strings.
dnl
dnl     This macro assumes version strings of the form
dnl     <major>.<minor>.<patch-level>. If `version-1' is less than `version-2',
dnl     `action-if-less' is executed, if they are equal, `action-if-equal' is
dnl     executed, and if `version-1 is greater than `version-2',
dnl     `action-if-greater' is executed.
dnl
dnl     NOTE: Use with care; not tested.
dnl

AC_DEFUN([AX_VERSION_COMPARE],
[
    m4_if(m4_cmp([$#], [2]), [-1],
          [m4_errprintn([Too few arguments to macro AX_VERSION_COMPARE.]) m4_exit([1])])

    version1=$1
    version2=$2
    action_if_less=$3
    action_if_equal=$4
    action_if_greater=$5

    version1_major=$(echo $version1 | awk -F. '{ print $[]1 }')
    version1_minor=$(echo $version1 | awk -F. '{ print $[]2 }')
    version1_patch=$(echo $version1 | awk -F. '{ print $[]3 }')

    version2_major=$(echo $version2 | awk -F. '{ print $[]1 }')
    version2_minor=$(echo $version2 | awk -F. '{ print $[]2 }')
    version2_patch=$(echo $version2 | awk -F. '{ print $[]3 }')

    version_major_cmp=$((version1_major - $version2_major))
    version_minor_cmp=$((version1_minor - $version2_minor))
    version_patch_cmp=$((version1_patch - $version2_patch))

    if ((version_major_cmp < 0)); then
        $action_if_less ;
    elif ((version_major_cmp > 0)); then
        $action_if_greater ;
    elif ((version_minor_cmp < 0)); then
        $action_if_less ;
    elif ((version_minor_cmp > 0)); then
        $action_if_greater ;
    elif ((version_patch_cmp < 0)); then
        $action_if_less ;
    elif ((version_patch_cmp > 0)); then
        $action_if_greater ;
    else
        $action_if_equal ;
    fi
])

dnl
dnl AX_SYSTEM_CPUS
dnl
dnl     Retrieve the number of cpus in the system
dnl
dnl     AC_SUBST([NR_CPUS])
dnl

AC_DEFUN([AX_SYSTEM_CPUS],
[
    AC_REQUIRE([AC_PROG_GREP])
    NR_CPUS=`cat /proc/cpuinfo | $GREP 'processor' | wc -l`
    AC_SUBST([NR_CPUS])
])

dnl
dnl AX_SYSTEM_DATE
dnl
dnl Get the current system date and time. This macro calls:
dnl
dnl AC_SUBST(CURRENT_TIME)
dnl AC_SUBST(CURRENT_MONTH)
dnl AC_SUBST(CURRENT_YEAR)
dnl

AC_DEFUN([AX_SYSTEM_DATE],
[
    AC_REQUIRE([AC_PROG_AWK])
    current_date=`date "+%R %Z,%B,%Y"`
    CURRENT_TIME=`echo $current_date | $AWK '{ split($[]0, a, ","); print a[[1]] }'`
    CURRENT_MONTH=`echo $current_date | $AWK '{ split($[]0, a, ","); print a[[2]] }'`
    CURRENT_YEAR=`echo $current_date | $AWK '{ split($[]0, a, ","); print a[[3]] }'`

    AC_SUBST(CURRENT_TIME)
    AC_SUBST(CURRENT_MONTH)
    AC_SUBST(CURRENT_YEAR)
])

AC_DEFUN([AX_SUMMARY_CONFIG],
[
    if test $ax_memory_nodes -gt 1; then
        is_numa="yes"
    else
        is_numa="no"
    fi

    AC_REQUIRE([AC_CANONICAL_BUILD])
    echo "*** CONFIGURATION SUMMARY ***"
    echo "    Target     $build"
    echo "    Build      $build_type"
    echo "    NUMA       $is_numa ($ax_memory_nodes memory node(s) detected)"
    echo "    CC         $CC"
    echo "    CXX        $CXX"
    echo "    CPPFLAGS   $CPPFLAGS $AX_CPPFLAGS"
    echo "    CFLAGS     $AX_CFLAGS $CFLAGS"
    echo "    CXXFLAGS   $AX_CXXFLAGS $CXXFLAGS"
    echo "    LDFLAGS    $AX_LDFLAGS $LDFLAGS"
    echo "    LIBS       $LIBS"
    echo "    PREFIX     $prefix"
    echo "    INDEX_TYPE $SPX_INDEX_TYPE"
    echo "    VALUE_TYPE $SPX_VALUE_TYPE"
])