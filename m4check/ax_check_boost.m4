dnl
dnl ax_check_boost.m4 -- Macros related to the Boost Library installation.
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl

dnl
dnl AX_CHECK_BOOST(required_version)
dnl
dnl     Check for a specific version of the C++ Boost library.
dnl
dnl     AC_SUBST(BOOST_CPPFLAGS)
dnl     AC_SUBST(BOOST_LDFLAGS)
dnl

AC_DEFUN([AX_CHECK_BOOST],
[
    AC_REQUIRE([AC_PROG_SED])
    m4_if([$#], [0],
          [m4_errprintn([Too few arguments to $0.]) m4_exit(1)])

    system_dirs="/usr/local /opt/local /usr /opt"
    boost_required_version=$1
    AC_MSG_CHECKING([for Boost library >= $boost_required_version])
    AC_ARG_WITH([boostdir],
                [AS_HELP_STRING([--with-boostdir=PREFIX],
                                [search for Boost library (headers and binaries) installed in PREFIX.])],
                [], [with_boostdir=check])

    AC_ARG_WITH([boostlibdir],
                [AS_HELP_STRING([--with-boostlibdir=DIR],
                                [search Boost library binaries in DIR.])],
                [], [])

     if test x"$with_boostdir" = x"check"; then
         search_dirs=$system_dirs
     else
         search_dirs=$with_boostdir
     fi

     CPPFLAGS_save="$CPPFLAGS"
     CXXFLAGS_save="$CXXFLAGS"
     AC_LANG_PUSH([C++])
     for d in $search_dirs; do
         CPPFLAGS=""
         CXXFLAGS=""
         AC_LINK_IFELSE([
             AC_LANG_PROGRAM([
                 #include <iostream>
                 #include "$d/include/boost/version.hpp"
                 ],
                 [ std::cout << BOOST_LIB_VERSION << std::endl; ])],
             [
                 boost_found=1
                 BOOST_CPPFLAGS="-I$d/include"
                 if test ! -z $with_boostlibdir; then
                     BOOST_LDFLAGS="-L$with_boostlibdir"
                 else
                     BOOST_LDFLAGS="-L$d/lib"
                 fi
                 boost_version=`./conftest$EXEEXT`
                 boost_version=`echo $boost_version | $SED -e 's/_/./g'`
                 AS_VERSION_COMPARE([$boost_version],
                                    [$boost_required_version],
                                    [boost_found=0])
             ],
             [ boost_found=0 ])
         if test $boost_found -eq 1; then
             break
         fi
     done
     AC_LANG_POP([C++])

     if test $boost_found -eq 0; then
         AC_MSG_RESULT([no])
         AC_MSG_ERROR([Could not find the required Boost version dnl
($boost_required_version). Please check your installation of Boost or dnl
try using the `--with-boostdir' option.])
     else
         AC_MSG_RESULT([yes (found version $boost_version)])
     fi

     # restore compiler flags
     CPPFLAGS="$CPPFLAGS_save"
     CXXFLAGS="$CXXFLAGS_save"

     AC_SUBST([BOOST_CPPFLAGS])
     AC_SUBST([BOOST_LDFLAGS])
])

dnl
dnl AX_CHECK_BOOST_LIB(boost_library)
dnl
dnl     Check for a specific Boost library
dnl
dnl     AC_SUBST(BOOST_<boost_library>_LIB)
dnl     AC_SUBST(BOOST_<boost_library>_LIB_STATIC)
dnl

AC_DEFUN([AX_CHECK_BOOST_LIB],
[
    AC_REQUIRE([AC_PROG_SED])
    AC_REQUIRE([AC_CANONICAL_BUILD])
    AC_REQUIRE([AC_CANONICAL_HOST])
    AC_REQUIRE([AX_SELECT_BUILD])

    m4_if([$#], [0],
          [m4_errprintn([Too few arguments to $0.]) m4_exit(1)])

    AC_ARG_WITH([boost_$1],
                [AS_HELP_STRING([--with-boost_$1@<:@=VARIANT@:>@],
                                [use the specified variant for this library.
                                 See naming convention for Boost libraries.])],
                [user_variant=$withval], [])

    check_lib=$1

    m4_if([$1], [program_options], [boost_release=1.32.0],
          [$1], [regex], [boost_release=1.18.0],
          [$1], [serialization], [boost_release=1.32.0],
          [$1], [system], [boost_release=1.35.0],
          [$1], [thread], [boost_release=1.25.0],
          [$1], [unit_test_framework], [boost_release=1.21.0],
          [boost_release=1.11])

    if test -z $BOOST_LDFLAGS; then
        # Check for Boost
        AX_CHECK_BOOST([$boost_release])
    fi

    AC_MSG_CHECKING([for Boost.$check_lib])
    
    boost_libdir=`echo $BOOST_LDFLAGS | $SED -e 's/-L//'`

    # set the correct extensions for library names
    is_darwin=`echo "$build_os" | $GREP -oe "darwin"`
    if test ! -z $is_darwin; then
        shared_ext=".dylib"
        static_ext=".a"
    else
        shared_ext=".so"
        static_ext=".a"
    fi

    if test ! -z $user_variant; then
        # user specified a variant
        libfiles_shared="$boost_libdir/libboost_${check_lib}${user_variant}$shared_ext"
        libfiles_static="$boost_libdir/libboost_${check_lib}${user_variant}$static_ext"
    else
        # inspect the system and the build type and prepend preferrable variants
        variants="-mt"       # prefer the thread-safe libraries by default
        if test x"$ax_build_mode" = x"debug"; then
            variants="-mt-d -d $variants"
        fi

        libfiles_shared=`/bin/ls $boost_libdir/libboost_${check_lib}*$shared_ext 2> /dev/null`
        libfiles_static=`/bin/ls $boost_libdir/libboost_${check_lib}*$static_ext 2> /dev/null`
        for v in $variants; do
            shared_=`/bin/ls $boost_libdir/libboost_${check_lib}*$v*$shared_ext 2> /dev/null`
            static_=`/bin/ls $boost_libdir/libboost_${check_lib}*$v*$static_ext 2> /dev/null`
            libfiles_shared="$shared_ $libfiles_shared"
            libfiles_static="$static_ $libfiles_static"
        done
    fi

    # save compiler flags
    CPPFLAGS_save=$CPPFLAGS
    CXXFLAGS_save=$CXXFLAGS
    LDLFAGS_save=$LDFLAGS
    LIBS_save=$LIBS

    AC_LANG_PUSH([C++])

    CPPFLAGS=$BOOST_CPPFLAGS
    CXXFLAGS=""
    LIBS=""
    lib_shared_found=0
    for libf in $libfiles_shared; do
        soname=`basename $libf`
        libname=`echo $soname | sed -e 's/^lib//' | sed -e 's/\..*//'`
        CXXFLAGS="-fPIC"
        LIBS="$libf"
        LDFLAGS="-shared"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([],[])],
                       [ lib_shared_found=1 ], [])
        if test $lib_shared_found -eq 1; then
            lib_shared="-l$libname"
            break
        fi
    done

    lib_static_found=0
    for libf in $libfiles_static; do
        CXXFLAGS=""
        LIBS="$libf"
        LDFLAGS=""
        AC_LINK_IFELSE([AC_LANG_PROGRAM([],[])],
                       [ lib_static_found=1 ], [])
        if test $lib_static_found -eq 1; then
            lib_static="$libf"
            break
        fi
    done

    AC_LANG_POP([C++])

    # restore compiler flags
    CPPFLAGS=$CPPFLAGS_save
    CXXFLAGS=$CXXFLAGS_save
    LDFLAGS=$LDFLAGS_save
    LIBS=$LIBS_save

    # check search results
    if test $lib_shared_found -eq 0 && test $lib_static_found -eq 0; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Boost.$check_lib was not found. dnl
Check your installation of the Boost library or, if you specified a variant dnl
with the `--with-boost_$1' option, either omit it or try another variant.])
    fi

    if test $lib_shared_found -eq 0 && test x"$ax_link_mode" = x"shared"; then
        AC_MSG_WARN([No shared library for $check_lib was found, dnl
but a shared library build for this package is requested.])
    fi

    if test $lib_static_found -eq 0 && test x"$ax_link_mode" = x"static"; then
        AC_MSG_WARN([No static library for $check_lib was found, dnl
but a static library build for this package is requested.])
    fi

    # everything is fine; set the variables and finish
    MAKE_VARIABLE([BOOST_$1], [LIB])=$lib_shared
    MAKE_VARIABLE([BOOST_$1], [LIB_STATIC])=$lib_static

    AC_SUBST(MAKE_VARIABLE([BOOST_$1], [LIB]))
    AC_SUBST(MAKE_VARIABLE([BOOST_$1], [LIB_STATIC]))

    # prepend libraries to LIBS
    dnl LIBS="$[]MAKE_VARIABLE([BOOST_$1], [LIB]) $LIBS"
    AC_MSG_RESULT([yes])
])
