/*
 * SpmvMethod.hpp -- SpMV function types.
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPMV_METHOD_H
#define SPMV_METHOD_H

#include "sparsex/internals/Vector.hpp"

typedef void (*spmv_fn_t)(void *matrix, vector_t *in, vector_t *out,
                          spx_scalar_t scale_f, vector_t *local_buff);
// typedef void (*spmv_sym_fn_t)(void *matrix, vector_t *in,
//                               vector_t *out, vector_t *temp, 
//                               spx_scalar_t scale_f);

#endif  // SPMV_METHOD_H
