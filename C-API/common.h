/**
 * libcsx/common.h -- Common utilities.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_COMMON_H__
#define LIBCSX_COMMON_H__

#include <stdlib.h>
#include "error.h"
#include "mattype.h"

#define INVALID_INPUT ((input_t *) NULL)
#define INVALID_MAT ((matrix_t *) NULL)
#define INVALID_VEC ((vector_t *) NULL)
#define INVALID_PERM ((perm_t *) NULL)

enum TuningOptions {
    OPT_SYMMETRIC               = 0x01,
    OPT_REORDER                 = 0x02,
    OPT_SPLIT_BLOCKS            = 0x04,
    OPT_ONE_DIM_BLOCKS          = 0x08,
    OPT_NR_THREADS              = 0x10,
    OPT_AFFINITY                = 0x20,
    OPT_XFORMS                  = 0x40,
    OPT_WINDOW_SIZE             = 0x80,
    OPT_NR_SAMPLES              = 0x100,
    OPT_SAMPLING_PORTION        = 0x200,
/* Timing related */
    OPT_ENABLE_PREPROC_TIMERS   = 0x400,
    OPT_ENABLE_SPMV_TIMERS      = 0x800,
/* Verbosity related */
    OPT_PRINT_ENCODING_STATS    = 0x1000
};
typedef unsigned int option_t;

static inline int
check_indexing(int base)
{
    return (base == 0 || base == 1);
}

static inline int
check_dim(index_t dim)
{
    return (dim >= 0);
}

/**
 *  \brief Library initialization routine.
 */
void libcsx_init();

/**
 *  \brief Library shutdown routine.
 */
void libcsx_close();

/**
 *  \brief malloc() wrapper.
 */
#define libcsx_malloc(type, size) \
	(type *)malloc_internal(size, __FILE__, __LINE__, __FUNCTION__)
void * malloc_internal(size_t x, const char *sourcefile, unsigned long lineno,
                       const char *function);

/**
 *  \brief free() wrapper.
 */
#define libcsx_free(object) \
    free_internal(object, __FILE__, __LINE__, __FUNCTION__)
void free_internal(void *ptr, const char *sourcefile, unsigned long lineno,
                   const char *function);

#endif // LIBCSX_COMMON_H__
