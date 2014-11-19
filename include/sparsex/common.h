/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * Copyright (C) 2014,      Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */ 

/**
 * \file common.h
 * \brief Common library utilities and definitions.
 * 
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_COMMON_H
#define SPARSEX_COMMON_H

#include <sparsex/error.h>
#include <sparsex/types.h>
#include <sparsex/internals/logger/LoggerUtil.hpp>
#include <sparsex/internals/Vector.hpp>
#include <stdlib.h>

/**
 *  The sparse matrix handle type.
 */
typedef struct matrix spx_matrix_t;

/**
 *  The input matrix handle type.
 */
typedef struct input spx_input_t;

/**
 *  The vector handle type.
 */
typedef struct vector_struct spx_vector_t;

/**
 *  The partitioning handle type.
 */
typedef struct partition spx_partition_t;

/**
 *  Dense array object that represents a permutation.
 */
typedef spx_index_t spx_perm_t;

/**
 *  Generic option type.
 */
typedef int spx_option_t;

/**
 *  Vector copy-mode type.
 */
typedef int spx_copymode_t;

/**
 *  @defgroup invalid_handles_group Invalid handle types
 *  @{
 */

/**
 *  Invalid input handle.
 */
#define SPX_INVALID_INPUT   ((spx_input_t *) NULL)
/**
 *  Invalid matrix handle.
 */
#define SPX_INVALID_MAT     ((spx_matrix_t *) NULL)
/**
 *  Invalid vector handle.
 */
#define SPX_INVALID_VEC     ((spx_vector_t *) NULL)
/**
 *  Invalid partitioning handle.
 */
#define SPX_INVALID_PART    ((spx_partition_t *) NULL)
/**
 *  Invalid permutation handle.
 */
#define SPX_INVALID_PERM    ((spx_perm_t *) NULL)
/**
 *  @}
 */

/**
 *  @defgroup options_group Available options
 *  @{
 */

/**
 *  Reorder the input matrix. @sa spx_mat_tune()
 */
#define SPX_MAT_REORDER          42

/**
 *  The user and the library agree to share a vector. The user promises not to 
 *  free, reallocate or modify the corresponding array (except through the
 *  interface's set-value routine, @sa spx_vec_set_entry()), while the library
 *  promises to directly update the input array values.
 */
#define SPX_VEC_SHARE            43

/**
 *  The library makes a copy of the input array. The user is henceforth allowed 
 *  to free, reallocate or modify the corresponding array, without affecting
 *  the previously created vector object, while the library will only modify
 *  its internal copy.
 */
#define SPX_VEC_COPY             44

/**
 *  Array indices start at 0.
 */
#define SPX_INDEX_ZERO_BASED     0

/**
 *  Array indices start at 1.
 */
#define SPX_INDEX_ONE_BASED      1
/**
 *  @}
 */

/**
 *  Checks indexing validity.
 */
static inline
int check_indexing(spx_option_t base)
{
    return (base == 0 || base == 1);
}

/**
 *  Checks copy mode validity.
 */
static inline
int check_copymode(spx_copymode_t mode)
{
    return (mode == SPX_VEC_SHARE || mode == SPX_VEC_COPY);
}


/**
 *  Checks validity of a matrix dimension.
 */
static inline
int check_mat_dim(spx_index_t dim)
{
    return (dim >= 0);
}

/**
 *  Checks compatibility of a vector's size.
 */
static inline
int check_vec_dim(const spx_vector_t *x, unsigned long dim)
{
    return (x->size == dim);
}

/**
 *  Disables logging in SparseX.
 */
void spx_log_disable_all();

/**
 *  Activates logging of the Error level on stderr.
 */
void spx_log_error_console();

/**
 *  Activates logging of the Warning level on stderr.
 */
void spx_log_warning_console();

/**
 *  Activates logging of the Info level on stderr.
 */
void spx_log_info_console();

/**
 *  Activates logging of the Verbose level on stderr.
 */
void spx_log_verbose_console();

/**
 *  Activates logging of the Debug level on stderr.
 */
void spx_log_debug_console();

/**
 *  Activates logging of the Error level on a file.
 *  The file name must be previously provided through
 *  the spx_log_set_file() routine. If not, a default
 *  "sparsex.log" file will be created. 
 */
void spx_log_error_file();

/**
 *  Activates logging of the Warning level on a file.
 *  The file name must be previously provided through
 *  the spx_log_set_file() routine. If not, a default
 *  "sparsex.log" file will be created. 
 */
void spx_log_warning_file();

/**
 *  Activates logging of the Info level on a file.
 *  The file name must be previously provided through
 *  the spx_log_set_file() routine. If not, a default
 *  "sparsex.log" file will be created. 
 */
void spx_log_info_file();

/**
 *  Activates logging of the Verbose level on a file.
 *  The file name must be previously provided through
 *  the spx_log_set_file() routine. If not, a default
 *  "sparsex.log" file will be created. 
 */
void spx_log_verbose_file();

/**
 *  Activates logging of the Debug level on a file.
 *  The file name must be previously provided through
 *  the spx_log_set_file() routine. If not, a default
 *  "sparsex.log" file will be created. 
 */
void spx_log_debug_file();

/**
 *  Activates all logging levels and redirects output to stderr.
 */
void spx_log_all_console();

/**
 *  Activates all logging levels and redirects output to a file.
 *  The file name must be previously provided through the
 *  spx_log_set_file() routine. If not, a default
 *  "sparsex.log" file will be created. If the file already
 *  exists it will be overwritten.
 *
 *  @param[in] file            a filename.
 */
void spx_log_all_file(const char *file);

/**
 *  Sets the file that will be used when logging is redirected
 *  to a file. If the file already exists it will be overwritten.
 *
 *  @param[in] file            a filename.
 */
void spx_log_set_file(const char *file);

/**
 *  Library initialization routine.
 */
void spx_init();

/**
 *  Library shutdown routine.
 */
void spx_finalize();

/**
 *  A malloc() wrapper.
 */
#define spx_malloc(type, size) \
	(type *) malloc_internal(size, __FILE__, __LINE__, __FUNCTION__)
void *malloc_internal(size_t x, const char *sourcefile, unsigned long lineno,
                      const char *function);

/**
 *  A free() wrapper.
 */
#define spx_free(object) \
    free_internal(object, __FILE__, __LINE__, __FUNCTION__)
void free_internal(void *ptr, const char *sourcefile, unsigned long lineno,
                   const char *function);

#endif /* SPARSEX_COMMON_H */
