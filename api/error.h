/**
 * libcsx/error.h -- Error handling interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_ERROR_H
#define LIBCSX_ERROR_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

/**
 *  \brief List of available generic error codes.
 */
#define LIBCSX_FAILURE             -1
#define LIBCSX_SUCCESS              0
#define LIBCSX_ERR_MIN_VALUE        1

#define LIBCSX_ERR_ARG_INVALID      2   /* invalid argument */
#define LIBCSX_ERR_FILE             3   /* generic file error */
#define LIBCSX_ERR_INPUT_MAT        4   /* input matrix wasn't properly created */
#define LIBCSX_ERR_TUNED_MAT        5   /* tuned matrix wasn't properly created */
#define LIBCSX_ERR_VEC_DIM          6   /* incompatible vector dimension */
#define LIBCSX_ERR_REORDER          7   /* reordering failed */
#define LIBCSX_ERR_ENTRY_NOT_FOUND  8   /* matrix entry not found */
#define LIBCSX_ERR_FILE_DATA        9   /* unexpected data in file */
#define LIBCSX_OUT_OF_BOUNDS        10  /* index out of bounds */
/* OS related errors */
#define LIBCSX_ERR_SYSTEM           15
#define LIBCSX_ERR_FILE_OPEN        16  /* unable to open file */
#define LIBCSX_ERR_FILE_READ        17  /* unable to read from file */
#define LIBCSX_ERR_FILE_WRITE       18  /* unable to write to file */
#define LIBCSX_ERR_MEM_ALLOC        19  /* memory allocation failed */
#define LIBCSX_ERR_MEM_FREE         20  /* memory deallocation failed */

#define LIBCSX_ERR_MAX_VALUE        21

#define LIBCSX_WARN_CSXFILE         22  /* no specific csxfile given */
#define LIBCSX_WARN_TUNING_OPT      23  /* invalid tuning option */
#define LIBCSX_WARN_RUNTIME_OPT     24  /* invalid runtime option */
#define LIBCSX_WARN_REORDER         25  /* reordering failed */
#define LIBCSX_WARN_MAX_VALUE       26

typedef int libcsx_error_t;

/**
 *  \brief Pointer to an error handling routine. Explain signature here!
 */
typedef void (*libcsx_errhandler_t) (libcsx_error_t, const char *, unsigned long,
                                     const char *, const char *, ...);

/**
 *  Macros that are called when an error has been detected. Both versions call
 *  the current error handling routine (either default or user-defined).
 *
 *  Note: The _0 version uses the default error message, as defined in
 *        libcsx_ErrorStrings[] (\see libcsx/error.c).                
 *
 *  @param[in] code     nonzero error code, see the list above
 *  @param[in] message  error message, if NULL a default error message is used
 */
#define SETERROR_0(code) \
    err_get_handler()(code, __FILE__, __LINE__, __FUNCTION__, NULL)
#define SETERROR_1(code, message) \
    err_get_handler()(code, __FILE__, __LINE__, __FUNCTION__, message)
#define SETWARNING(code) \
    err_get_handler()(code, __FILE__, __LINE__, __FUNCTION__, NULL)

/**
 *  Default error handler that is called through the macros #SETERROR_0/1
 *  when an error has been detected. Prints the error messages to stderr
 *  by default.                  
 *
 *  @param[in] code             nonzero error code, see the list above.
 *  @param[in] sourcefile       the name of the file where the error occured.
 *  @param[in] lineno           the line in the file where the error occured.
 *  @param[in] function         the function in which the error occured.
 *  @param[in] fmt              error message.
 */
void err_handle(libcsx_error_t code, const char *sourcefile, unsigned long lineno,
                const char *function, const char *fmt, ...);

/**
 *  Returns a pointer to the current error handler (either default or
 *  user-defined).
 *
 *  @return     current error handling routine.
 */
libcsx_errhandler_t err_get_handler();

/**
 *  This function allows the user to change the default error handling
 *  policy with a new one, which must conform to the signature provided above.
 *
 *  @param[in] new_handler   user-defined routine.
 */
void err_set_handler(libcsx_errhandler_t new_handler);

#endif // LIBCSX_ERROR_H
