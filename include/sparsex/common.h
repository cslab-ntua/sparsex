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
 *  Vector data-mode type.
 */
typedef unsigned int spx_vecmode_t;

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
#define SPX_MAT_REORDER         42

/**
 *  @sa spx_vec_create_from_buff(). When this option is set,
 *  the buffer provided by the user will be used as is in the newly
 *  created vector handle.
 */
#define SPX_VEC_AS_IS           43

/**
 *  @sa spx_vec_create_from_buff(). When this option is set,
 *  the buffer provided by the user will be copied into an optimally
 *  allocated buffer, that will be returned to the user for further use.
 */
#define SPX_VEC_TUNE            44

/**
 *  Array indices start at 0.
 */
#define SPX_INDEX_ZERO_BASED    45

/**
 *  Array indices start at 1.
 */
#define SPX_INDEX_ONE_BASED     46
/**
 *  @}
 */

/**
 *  Checks indexing validity.
 */
static inline
int check_indexing(spx_option_t base)
{
    return (base == SPX_INDEX_ZERO_BASED || base == SPX_INDEX_ONE_BASED);
}

/**
 *  Checks data mode validity for vectors created with 
 *  spx_vec_create_from_buff().
 */
static inline
int check_vecmode(spx_vecmode_t mode)
{
    return (mode == SPX_VEC_AS_IS || mode == SPX_VEC_TUNE);
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
	(type *) malloc_internal(size, __FILE__, __LINE__, __func__)
void *malloc_internal(size_t x, const char *sourcefile, unsigned long lineno,
                      const char *function);

/**
 *  A free() wrapper.
 */
#define spx_free(object) \
    free_internal(object, __FILE__, __LINE__, __func__)
void free_internal(void *ptr, const char *sourcefile, unsigned long lineno,
                   const char *function);

#endif /* SPARSEX_COMMON_H */
