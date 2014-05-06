dnl
dnl platforms.m4 -- Macros related to system platforms.
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl

dnl
dnl AX_PLATFORM
dnl
dnl     Returns in `ax_platform' a simple string describing the underlying
dnl     platform. The returned value may be one of `linux' or `darwin'. If the
dnl     underlying platform is not any of these, it returns the value of the
dnl     `build_os' variable.
dnl

AC_DEFUN([AX_PLATFORM],
[
AC_REQUIRE([AC_CANONICAL_BUILD])
AC_REQUIRE([AC_CANONICAL_HOST])

case $build_os in
    linux*)
        ax_platform="linux" ;;
    darwin*)
        ax_platform="darwin" ;;
    *)
        ax_platform=$build_os ;;
esac
])
