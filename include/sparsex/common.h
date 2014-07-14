/**
 * \file common.h
 *
 * \brief Common library utilities and definitions
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_COMMON_H
#define SPARSEX_COMMON_H

#include <sparsex/error.h>
#include <sparsex/types.h>
#include <sparsex/internals/logger/LoggerUtil.hpp>
#include <sparsex/internals/Vector.hpp>

#include <stdlib.h>

/**
 *  \brief The sparse matrix handle type.
 */
typedef struct matrix spx_matrix_t;

/**
 *  \brief The input matrix handle type.
 */
typedef struct input spx_input_t;

/**
 *  \brief The vector handle type.
 */
typedef struct vector_struct spx_vector_t;

/**
 *  \brief Partitioning object.
 */
typedef struct partition spx_partition_t;

/**
 *  \brief Dense array object that represents a permutation.
 */
typedef spx_index_t spx_perm_t;

/**
 *  \brief Matrix property type.
 */
typedef int spx_property_t;

/**
 *  \brief A vector copy-mode type.
 */
typedef int spx_copymode_t;

#define INVALID_INPUT   ((spx_input_t *) NULL)
#define INVALID_MAT     ((spx_matrix_t *) NULL)
#define INVALID_VEC     ((spx_vector_t *) NULL)
#define INVALID_PART    ((spx_partition_t *) NULL)
#define INVALID_PERM    ((spx_perm_t *) NULL)

#define OP_REORDER          42
#define OP_SHARE            43
#define OP_COPY             44
#define INDEXING_ZERO_BASED 0
#define INDEXING_ONE_BASED  1

static inline
int check_indexing(int base)
{
    return (base == 0 || base == 1);
}

static inline
int check_mat_dim(spx_index_t dim)
{
    return (dim >= 0);
}

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
 *  \brief Library initialization routine.
 */
void spx_init();

/**
 *  \brief Library shutdown routine.
 */
void spx_finalize();

/**
 *  \brief malloc() wrapper.
 */
#define spx_malloc(type, size) \
	(type *) malloc_internal(size, __FILE__, __LINE__, __FUNCTION__)
void *malloc_internal(size_t x, const char *sourcefile, unsigned long lineno,
                      const char *function);

/**
 *  \brief free() wrapper.
 */
#define spx_free(object) \
    free_internal(object, __FILE__, __LINE__, __FUNCTION__)
void free_internal(void *ptr, const char *sourcefile, unsigned long lineno,
                   const char *function);

#endif /* SPARSEX_COMMON_H */
