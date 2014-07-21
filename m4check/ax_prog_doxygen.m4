dnl
dnl ax_prog_doxygen.m4 -- Check for and set up doxygen
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl 

AC_DEFUN([AX_PROG_DOXYGEN],
[
AX_CHECK_PROG([doxygen])
if test -z $DOXYGEN; then
   AC_MSG_WARN([Could not find `doxygen' in your path. dnl
You may not be able to generate the documentation of the package.])
else
    dx_docdir=doc
    dx_api_subdir=api
    dx_devel_subdir=devel

    DX_CONF_API=Doxyfile.api
    DX_CONF_DEVEL=Doxyfile.devel
    DX_CONF_TAGFILE=$PACKAGE_TARNAME.tag
    DX_CONF_HTML_DIR=html

    AC_SUBST(dx_docdir)
    AC_SUBST(dx_api_subdir)
    AC_SUBST(dx_devel_subdir)
    AC_SUBST(DX_CONF_API)
    AC_SUBST(DX_CONF_DEVEL)
    AC_SUBST(DX_CONF_HTML_DIR)
    AC_SUBST(DX_CONF_TAGFILE)
fi

AM_CONDITIONAL([SPX_GENERATE_DOC], [ test ! -z $DOXYGEN ])
])

