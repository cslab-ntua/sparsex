/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file types.h
 * \brief Available indexing and value types.
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \author Athena Elafrou
 * \author Vasileios Karakasis
 * \date 2013&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_TYPES_H
#define SPARSEX_TYPES_H

#include <sparsex/config.h>

#ifdef SPX_INDEX_TYPE
typedef SPX_INDEX_TYPE spx_index_t;
#else
typedef int spx_index_t;
#endif

#ifdef SPX_VALUE_TYPE
typedef SPX_VALUE_TYPE spx_value_t;
#else
typedef double spx_value_t;
#endif

#endif /* SPARSEX_TYPES_H */
