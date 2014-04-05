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

#include "Vector.hpp"

typedef void spmv_double_fn_t(void *matrix, spx_vector_t *in, spx_vector_t *out,
                              spx_scalar_t scale_f);
typedef void spmv_double_sym_fn_t(void *matrix, spx_vector_t *in,
                                  spx_vector_t *out, spx_vector_t *temp, 
                                  spx_scalar_t scale_f);

#endif  // SPMV_METHOD_H
