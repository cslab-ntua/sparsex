/*
 * Copyright (C) 2011-2016, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2016, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file futex.h
 * \brief Linux futex wrapper
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2016
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/futex.h>
#include <errno.h>
#include <linux/futex.h>
#include <sys/syscall.h>
#include <unistd.h>

void futex_wait(int *addr, int val)
{
    int err = syscall(SYS_futex, addr, FUTEX_WAIT, val, NULL);
    if (err < 0 && errno == ENOSYS)
        syscall(SYS_futex, addr, FUTEX_WAIT, val, NULL);
}

void futex_wake(int *addr, int count)
{
    // Wakes at most count processes waiting on the addr
    int err = syscall(SYS_futex, addr, FUTEX_WAKE, count);
    if (err < 0 && errno == ENOSYS)
        syscall(SYS_futex, addr, FUTEX_WAKE, count);
}
