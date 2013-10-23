/**
 * libcsx/common.c -- Common utilities.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "common.h"

void libcsx_init()
{
    /* log_disable_all(); */
    log_enable_all_file("logfile");
//    log_enable_all_console();
    log_disable_debug();
}

void *malloc_internal(size_t x, const char *sourcefile, unsigned long lineno,
                      const char *function)
{
    void *ret;
    ret = malloc(x);
    if (!ret) {
        err_handle(LIBCSX_ERR_MEM_ALLOC, sourcefile, lineno, function, NULL);
        exit(1);
    }
    return ret;
}

void free_internal(void *ptr, const char *sourcefile, unsigned long lineno,
                   const char *function)
{
    if (!ptr) {
        err_handle(LIBCSX_ERR_MEM_FREE, sourcefile, lineno, function, NULL);
        exit(1);
    }
    free(ptr);
    ptr = NULL;
}
