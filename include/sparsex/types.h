/**
 * \file types.h -- Available indexing and value types.
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
#   undef spx_index_t
#   define spx_index_t uint32_t
#endif

#ifdef USE_64BIT_INDICES
#   undef spx_index_t
#   define spx_index_t uint64_t
#endif

#ifndef spx_index_t
#   define spx_index_t int
#endif

#ifdef USE_SINGLE_PRECISION
#   define spx_value_t float
#else
#   define spx_value_t double
#endif

#ifndef spx_scalar_t
#   define spx_scalar_t double
#endif

#endif // SPARSEX_MATTYPE_H
