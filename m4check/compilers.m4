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

ax_cxx=$1
if test -z $ax_cxx; then
   ax_cxx=$CXX
fi

if test x"$ax_cxx" == x"g++"; then
    # g++ flags
	case "$ax_build_mode" in
	    "debug")
            # undefine DEBUG; LLVM doesn't want it
            ax_pkg_cppflags="-UNDEBUG -DSPX_DEBUG=1"
            ax_pkg_cxxflags="-g -O0 -std=c++0x dnl
-fstrict-aliasing -Wall -Woverloaded-virtual"
            ax_pkg_cflags="-g -O0 -std=c99 -Wall -fstrict-aliasing" ;;
        "release")
            ax_pkg_cppflags="-DNDEBUG"
            ax_pkg_cxxflags="-O3 -std=c++0x"
            ax_pkg_cflags="-O3 -std=c99" ;;
		*)
			AC_MSG_ERROR(
				[@<:@BUG@:>@ should not have entered here: m4_location]) ;;
    esac

	if test x"$ax_link_mode" = x"static"; then
	    ax_pkg_ldflags="-all-static -pthread -Wl,--allow-multiple-definition"
    else
        ax_pkg_ldflags="-Wl,--allow-multiple-definition"
	fi
fi
])
