/*
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file SpmvMethod.hpp
 * \brief SpMV function types
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_SPMV_METHOD_H
#define SPARSEX_INTERNALS_SPMV_METHOD_H

#include <sparsex/internals/Vector.hpp>

typedef void (*spmv_fn_t)(const void *matrix, const vector_t *in, vector_t *out,
                          spx_value_t scale_f, vector_t *local_buff);

#endif  // SPARSEX_INTERNALS_SPMV_METHOD_H
