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

#ifndef SPARSEX_INTERNALS_SPMV_METHOD_H
#define SPARSEX_INTERNALS_SPMV_METHOD_H

#include "sparsex/internals/Vector.hpp"

typedef void spmv_double_fn_t(void *matrix, vector_t *in, vector_t *out,
                              spx_scalar_t scale_f);
typedef void spmv_double_sym_fn_t(void *matrix, vector_t *in,
                                  vector_t *out, vector_t *temp, 
                                  spx_scalar_t scale_f);

#endif  // SPARSEX_INTERNALS_SPMV_METHOD_H
