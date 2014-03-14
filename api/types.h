/**
 * SparseX/types.h -- Available indexing and value types.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSEX_TYPES_H
#define SPARSEX_TYPES_H

#include <inttypes.h>

#ifdef USE_UNSIGNED_INDICES
#   undef index_t
#   define index_t uint32_t
#endif

#ifdef USE_64BIT_INDICES
#   undef index_t
#   define index_t uint64_t
#endif

#ifndef index_t
#   define index_t int
#endif

#ifdef USE_SINGLE_PRECISION
#   define value_t float
#else
#   define value_t double
#endif

#ifndef scalar_t
#   define scalar_t double
#endif

#endif // SPARSEX_MATTYPE_H
