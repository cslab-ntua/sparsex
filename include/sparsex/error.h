/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * Copyright (C) 2014,      Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */ 

/**
 * \file error.h
 * \brief Error handling interface.
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_ERROR_H
#define SPARSEX_ERROR_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

/**
 *  @defgroup errorcode_group List of available generic error codes.
 *  @{
 */
#define SPX_FAILURE             -1  /** failure */
#define SPX_SUCCESS             0   /** success */
#define SPX_ERR_MIN_VALUE       1
#define SPX_ERR_ARG_INVALID     2   /** invalid argument */
#define SPX_ERR_FILE            3   /** generic file error */
#define SPX_ERR_INPUT_MAT       4   /** input matrix wasn't properly created */
#define SPX_ERR_TUNED_MAT       5   /** tuned matrix wasn't properly created */
#define SPX_ERR_VEC             6   /** vector creation failed */
#define SPX_ERR_PART            7   /** partitioning object wasn't properly
                                        created */
#define SPX_ERR_PERM            8   /** error in permutation */
#define SPX_ERR_DIM             9   /** incompatible matrix and vectors
                                         dimensions */
#define SPX_ERR_VEC_DIM         10   /** incompatible vector dimension */
#define SPX_ERR_ENTRY_NOT_FOUND 11  /** matrix entry not found */
#define SPX_OUT_OF_BOUNDS       12  /** index out of bounds */
#define SPX_ERR_SYSTEM          15
#define SPX_ERR_FILE_OPEN       16  /** unable to open file */
#define SPX_ERR_FILE_READ       17  /** unable to read from file */
#define SPX_ERR_FILE_WRITE      18  /** unable to write to file */
#define SPX_ERR_MEM_ALLOC       19  /** memory allocation failed */
#define SPX_ERR_MEM_FREE        20  /** memory deallocation failed */
#define SPX_ERR_MAX_VALUE       21
/**
 *  @}
 */

/**
 *  @defgroup warncode_group List of available warning codes.
 *  @{
 */
#define SPX_WARN_CSXFILE         22  /** no specific csx binary file given */
#define SPX_WARN_TUNING_OPT      23  /** invalid tuning option */
#define SPX_WARN_RUNTIME_OPT     24  /** invalid runtime option */
#define SPX_WARN_REORDER         25  /** reordering failed */
#define SPX_WARN_ENTRY_NOT_SET   26  /** entry not set in matrix or vector */
#define SPX_WARN_MAX_VALUE       27
/**
 *  @}
 */

/**
 *  Error type.
 */
typedef int spx_error_t;

/**
 *  Pointer to an error handling routine.
 */
typedef void (*spx_errhandler_t) (spx_error_t, const char *, unsigned long,
                                  const char *, const char *, ...);

/**
 *  @defgroup seterror_group Convenient macros for logging errors and warnings
 *  Macros that are called when an error has been detected. Both versions call
 *  the current error handling routine (either default or user-defined).
 *
 *  Note: The _0 version uses the default error message, as defined in
 *        spx_ErrorStrings[] (\see SparseX/error.c).                
 */

/**
 *  @ingroup seterror_group 
 *  @param[in] code     nonzero error code, see the list above
 */
#define SETERROR_0(code) \
    spx_err_get_handler()(code, __FILE__, __LINE__, __func__, NULL)

/**
 *  @ingroup seterror_group 
 *  @param[in] code     nonzero error code, see the list above
 *  @param[in] message  error message, if NULL a default error message is used
 */
#define SETERROR_1(code, message) \
    spx_err_get_handler()(code, __FILE__, __LINE__, __func__, message)

/**
 *  @ingroup seterror_group 
 *  @param[in] code     nonzero error code, see the list above
 */
#define SETWARNING(code) \
    spx_err_get_handler()(code, __FILE__, __LINE__, __func__, NULL)

/**
 *  Default error handler that is called through the macros #SETERROR_0/1
 *  when an error has been detected. Prints the error messages to stderr
 *  by default.                  
 *
 *  @param[in] code             nonzero error code, see the list above.
 *  @param[in] sourcefile       the name of the file where the error occurred.
 *  @param[in] lineno           the line in the file where the error occurred.
 *  @param[in] function         the function in which the error occurred.
 *  @param[in] fmt              error message.
 */
void err_handle(spx_error_t code, const char *sourcefile, unsigned long lineno,
                const char *function, const char *fmt, ...);

/**
 *  Returns a pointer to the current error handler (either default or
 *  user-defined).
 *
 *  @return     current error handling routine.
 */
spx_errhandler_t spx_err_get_handler();

/**
 *  This function allows the user to change the default error handling
 *  policy with a new one, which must conform to the signature provided above.
 *
 *  @param[in] new_handler   user-defined routine.
 */
void spx_err_set_handler(spx_errhandler_t new_handler);

#endif /* SPARSEX_ERROR_H */
