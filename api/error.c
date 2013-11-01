/**
 * libcsx/error.c -- Error handling interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "error.h"
#include "LoggerUtil.hpp"

/**
 *  The first global variable points to the current error handler. The default
 *  error handler #err_handle() can be overidden with a user defined error
 *  handling routine by setting this variable through #err_set_handler(). The new
 *  handler must conform to the given signature (\see libcsx/error.h).
 */
static libcsx_errhandler_t err_handlerptr = err_handle;

#define MAXLINE 1024

/**
 *  Returns the default error message corresponding to supplied error code
 *  as a string.
 *
 *  @param[in]  code         nonzero error code, see the list above.
 *  @param[out] message      default message as defined in libcsxErrorStrings[].
 */
static const char *err_get_message(libcsx_error_t code);

/**
 *  Helper function that prints a message and returns to caller.
 *  Printing modes:
 *      1. errnoflag set to 0: for library related errors/warnings.
 *         The corresponding default or user-defined error message
 *         is simply printed to the chosen output (logfile/console).
 *      2. errnoflag set to 1: for OS related errors.
 *         The OS internal error message corresponding to the errno
 *         value is added to the default or user-defined error message.
 *         The errno value must be given as the "error" parameter.
 *
 *  @param[in] levelflag    flag that defines if the message is an error or
 *                          warning (0 for warning, 1 for error).
 *  @param[in] errnoflag    flag that defines the error handling mode.
 *  @param[in] error        nonzero error code.
 *  @param[in] sourcefile   the name of the file where the error occured.
 *  @param[in] lineno       the line in the file where the error occured.
 *  @param[in] function     the function in which the error occured.
 *  @param[in] fmt          default or user-defined error message.
 *  @param[in] ap           list f addintional arguments.          
 */
static void err_print(int levelflag, int errnoflag, int error,
                      const char *sourcefile, unsigned long lineno,
                      const char *function, const char *fmt, va_list ap);

void err_handle(libcsx_error_t code, const char *sourcefile,
                unsigned long lineno, const char *function, const char *fmt,
                ...)
{
    va_list ap;
    int errno_saved = errno;

    /* In case a valid error code is given find the corresponding default
       error/warning message */
    if (!fmt) {
        fmt = err_get_message(code);
    }
    assert(fmt && "Non-existent error code!");
    va_start(ap, fmt);
    if (code > LIBCSX_ERR_MIN_VALUE && code < LIBCSX_ERR_SYSTEM) {
        err_print(1, 0, 0, sourcefile, lineno, function, fmt, ap);
    } else if (code > LIBCSX_ERR_SYSTEM && code < LIBCSX_ERR_MAX_VALUE) {
        err_print(1, 1, errno_saved, sourcefile, lineno, function, fmt, ap);
        exit(1);
    } else if (code > LIBCSX_ERR_MAX_VALUE && code < LIBCSX_WARN_MAX_VALUE) {
        err_print(0, 0, 0, sourcefile, lineno, function, fmt, ap);
    }
    va_end(ap);
}

libcsx_errhandler_t err_get_handler()
{
    return err_handlerptr;
}

void err_set_handler(libcsx_errhandler_t new_handler)
{
    /* If NULL handler is set to default */
    if (new_handler) {
        err_handlerptr = new_handler;
    } else {
        err_handlerptr = err_handle;
    }
}

/**
 *  \brief Strings corresponding to the defined generic error codes.
 */
static const char *libcsx_errors[] = {
    "invalid argument",
    "file doesn't exist",
    "loading of input matrix failed",
    "conversion to CSX failed",
    "incompatible vector dimension",
    "reordering failed",
    "matrix entry doesn't exist",
    "unexpected data in file",
    "index out of bounds",
    "dummy",
    "failed to open file",
    "failed to read from file",
    "failed to write to file",
    "memory allocation failed",
    "memory deallocation failed"
};

/**
 *  \brief Strings corresponding to the defined generic warning codes.
 */
static const char *libcsx_warnings[] = {
    "no specific file given to save CSX, using default: \"csx_file\"",
    "invalid tuning option",
    "invalid runtime option",
    "reordering wasn't feasible on this matrix"
};

static const char *err_get_message(libcsx_error_t code)
{
    if (code > LIBCSX_ERR_MIN_VALUE && code < LIBCSX_ERR_MAX_VALUE) {
        return libcsx_errors[code - LIBCSX_ERR_MIN_VALUE - 1];
    } else if (code > LIBCSX_ERR_MAX_VALUE && code < LIBCSX_WARN_MAX_VALUE) {
        return libcsx_warnings[code - LIBCSX_ERR_MAX_VALUE - 1];
    } else {
        return NULL;
    }
    return NULL;
}

static void err_print(int levelflag, int errnoflag, int error,
                      const char *sourcefile, unsigned long lineno,
                      const char *function, const char *fmt, va_list ap)
{
    char buf[MAXLINE] = "";
    char *postfix = " [\"%s\":%ld:%s()]\n";

    /* Add default or user defined error/warning message */
    vsnprintf(buf+strlen(buf), MAXLINE-strlen(buf), fmt, ap);
    /* Add os error message */
    if (errnoflag)
        snprintf(buf+strlen(buf), MAXLINE-strlen(buf), ": %s", strerror(error));
    /* Add postfix */
    snprintf(buf+strlen(buf), MAXLINE-strlen(buf), postfix, sourcefile, lineno,
             function);
    if (levelflag) {
        log_error(buf);
    } else {
        log_warning(buf);
    }
}
