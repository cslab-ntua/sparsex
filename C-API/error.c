/**
 *
 * libcsx/error.c -- Error handling interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "error.h"

/**
 *  The first global variable points to the current error handler. The default
 *  error handler #err_handle() can be overidden with a user defined error
 *  handling routine by setting this variable through #err_set_handler(). The new
 *  handler must conform to the given signature (\see libcsx/error.h).
 *
 *  The second global variable defines the current error log file i.e. where
 *  error messages are dumped to. The default option is stderr. Can be
 *  overridden with #err_set_logfile().
 */
static libcsx_errhandler_t err_handlerptr = err_handle;
static libcsx_logfile_t err_logfile;

#define MAXLINE 1024

/**
 *  Helper function that prints a message and returns to caller.
 *  Printing modes:
 *      1. errnoflag set to 0: for library related errors/warnings.
 *         The corresponding default or user-defined error message
 *         is simply printed to the logfile.
 *      2. errnoflag set to 1: for OS related errors.
 *         The OS internal error message corresponding to the errno
 *         value is added to the default or user-defined error message.
 *         The errno value must be given as the "error" parameter.
 *
 *  @param[in] errnoflag    flag that defines the error handling mode.
 *  @param[in] error        nonzero error code.
 *  @param[in] sourcefile   the name of the file where the error occured.
 *  @param[in] lineno       the line in the file where the error occured.
 *  @param[in] function     the function in which the error occured.
 *  @param[in] prefix       [ERROR] or [WARNING].
 *  @param[in] fmt          default or user-defined error message.
 *  @param[in] ap           list f addintional arguments.          
 */
static void err_print(int errnoflag, int error, const char* sourcefile,
                      unsigned long lineno, const char *function,
                      const char *prefix, const char *fmt, va_list ap)
{
    char buf[MAXLINE] = "";
    char *postfix = " [\"%s\":%ld:%s()]\n";

    /* Add prefix */
    strcat(buf, prefix);
    /* Add default or user defined error/warning message */
    vsnprintf(buf+strlen(buf), MAXLINE-strlen(buf), fmt, ap);
    /* Add os error message */
    if (errnoflag)
        snprintf(buf+strlen(buf), MAXLINE-strlen(buf), ": %s", strerror(error));
    /* Add postfix */
    snprintf(buf+strlen(buf), MAXLINE-strlen(buf), postfix, sourcefile, lineno,
             function);
    fputs(buf, err_logfile);
    fflush(err_logfile);
}

void err_handle(libcsx_error_t code, const char *sourcefile,
                unsigned long lineno, const char *function, const char *fmt,
                ...)
{
    va_list ap;
    int errno_saved = errno;
    char *err_prefix = "[ERROR]: ";
    char *warn_prefix = "[WARNING]: ";

    /* In case a valid error code is given find the corresponding default
       error/warning message */
    if (!fmt) {
        err_get_message(code, &fmt);
    }
    assert(fmt && "Non-existent error code!");
    va_start(ap, fmt);                                                      //FIXTHIS ap
    if (code > LIBCSX_ERR_MIN_VALUE && code < LIBCSX_ERR_SYSTEM) {
        err_print(0, 0, sourcefile, lineno, function, err_prefix, fmt, ap);
    } else if (code > LIBCSX_ERR_SYSTEM && code < LIBCSX_ERR_MAX_VALUE) {
        err_print(1, errno_saved, sourcefile, lineno, function, err_prefix,
                  fmt, ap);
        exit(1);
    } else if (code > LIBCSX_ERR_MAX_VALUE && code < LIBCSX_WARN_MAX_VALUE) {
        err_print(0, 0, sourcefile, lineno, function, warn_prefix, fmt, ap);
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

void err_set_logfile(const char *filename)
{
    if (!filename) {
        err_logfile = stderr;
    } else {
        err_logfile = fopen(filename, "a+");
        if (err_logfile == NULL) {
            fprintf(stderr, "[ERROR]: failed to open logfile\n");
            exit(1);
        }
    }
}

void err_close_logfile()
{
    int error = fclose(err_logfile);
    if (error != 0) printf("ok\n");
    err_logfile = NULL;
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
    "reoredering wasn't feasible on this matrix"
};

libcsx_error_t err_get_message(libcsx_error_t code, const char *message[])
{
    if (message && code > LIBCSX_ERR_MIN_VALUE && code < LIBCSX_ERR_MAX_VALUE) {
        *message = libcsx_errors[code - LIBCSX_ERR_MIN_VALUE - 1];
    } else if (message && code > LIBCSX_ERR_MAX_VALUE && code < LIBCSX_WARN_MAX_VALUE) {
        *message = libcsx_warnings[code - LIBCSX_ERR_MAX_VALUE - 1];
    } else if (message) {
        *message = NULL;
    }
    return 0;
}

#if 0
/* Error testing */
int main()
{
    error_logfile= stderr;
    FILE *p;
    char text[] = "blabla";
    p = fopen(text, "r");

    if (p == NULL) {
        SETERROR_0(LIBCSX_ERR_FILE_OPEN);
    }
    return 0;
}
#endif
