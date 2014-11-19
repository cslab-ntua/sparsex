dnl
dnl ax_check_prog.m4 -- Check for a specific program
dnl
dnl Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
dnl Copyright (C) 2014, Vasileios Karakasis
dnl All rights reserved.
dnl
dnl This file is distributed under the BSD License. See LICENSE.txt for details.
dnl

dnl
dnl AX_CHECK_PROG(prog, [variable], [path = $PATH])
dnl
dnl     Check for `prog' in `path' and set `variable' to the absolute filename
dnl     of the program, if found. Path directories must not contain spaces.
dnl
dnl     AC_SUBST([variable])
dnl
AC_DEFUN([AX_CHECK_PROG],
[
    AC_REQUIRE([AC_PROG_SED])
    m4_if([$#], [0],
          [m4_errprintn([Too few arguments to $0.]) m4_exit(1)])
       
    user_prog=$1
    user_path=$3

    AC_MSG_CHECKING([$user_prog])
    if test -z $user_path; then
        user_path=$PATH
    fi

    search_dirs=`echo $user_path | $SED -e 's/:/ /g'`
    for d in $search_dirs; do
        if test -e "$d/$user_prog" && test -x "$d/$user_prog"; then
            prog_found="$d/$user_prog"
            break
        fi
    done

    if test -z $prog_found; then
        AC_MSG_RESULT([no])
    else
        m4_ifblank($2, UPPERCASE($1), $2)=$prog_found
        AC_SUBST(m4_ifblank($2, UPPERCASE($1), $2))
        AC_MSG_RESULT([yes])
    fi
])
