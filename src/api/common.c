/**
 * \file common.c
 *
 * \brief Common utilities
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/common.h>

void spx_log_disable_all()
{
    DisableLogging();
}

void spx_log_error_console()
{
    EnableErrorConsole();
}

void spx_log_warning_console()
{
    EnableWarningConsole();
}

void spx_log_info_console()
{
    EnableInfoConsole();
}

void spx_log_debug_console()
{
    EnableDebugConsole();
}

void spx_log_error_file()
{
    EnableErrorFile();
}

void spx_log_warning_file()
{
    EnableWarningFile();
}

void spx_log_info_file()
{
    EnableInfoFile();
}

void spx_log_debug_file()
{
    EnableDebugFile();
}

void spx_log_all_console()
{
    EnableAllConsole();
}

void spx_log_all_file(const char *file)
{
    EnableAllFile(file);
}

void spx_log_set_file(const char *file)
{
    UseFile(file);
}

void spx_init()
{
    EnableWarningConsole();
}

void spx_finalize()
{
    /* Cleanups */
}

void *malloc_internal(size_t x, const char *sourcefile, unsigned long lineno,
                      const char *function)
{
    void *ret;
    ret = malloc(x);
    if (!ret) {
        err_handle(SPX_ERR_MEM_ALLOC, sourcefile, lineno, function, NULL);
        exit(1);
    }

    return ret;
}

void free_internal(void *ptr, const char *sourcefile, unsigned long lineno,
                   const char *function)
{
    if (!ptr) {
        err_handle(SPX_ERR_MEM_FREE, sourcefile, lineno, function, NULL);
        exit(1);
    }

    free(ptr);
    ptr = NULL;
}
