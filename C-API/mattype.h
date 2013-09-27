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

#ifndef SCALAR
#define SCALAR double
#endif

/* Default index and value types of a sparse matrix */
#define index_t uint64_t
#define value_t double

#ifdef USE_UNSIGNED_INDICES
#define index_t uint32_t
#endif

#ifdef USE_64BIT_INDICES
#define index_t uint64_t
#endif

#ifdef USE_SINGLE_PRECISION
#define value_t float
#else
#define value_t double
#endif

#endif // LIBCSX_MATTYPE_H__
