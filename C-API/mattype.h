/**
 * libcsx/mattype.h -- Available matrix types.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_MATTYPE_H__
#define LIBCSX_MATTYPE_H__

#include <inttypes.h>

/* Default index and value types of a sparse matrix */
#define index_t int
#define value_t double

#if defined(USE_UNSIGNED_INDICES)
#define index_t uint32_t
#endif

#if defined(USE_64BIT_INDICES)
#define index_t uint64_t
#endif

#if defined(USE_SINGLE_PRECISION)
#define value_t float
#elif defined(USE_DOUBLE_PRECISION)
#define value_t double
#endif

#endif // LIBCSX_MATTYPE_H__
