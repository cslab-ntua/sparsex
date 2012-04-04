/*
 * cg_vector.h -- Interface of Expansion of Vector Operations for CG.
 *
 * Copyright (C) 2011,      Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
 
#ifndef CG_VECTOR_H_
#define CG_VECTOR_H_

#include <stdint.h>
#include <assert.h>
#include <cstring>

extern "C" {
#include "vector.h"
}

void vector_double_sub(vector_double_t *in1, vector_double_t *in2,
                       vector_double_t *out);
void vector_double_sub_part(vector_double_t *in1, vector_double_t *in2,
                            vector_double_t *out, uint64_t start, uint64_t end);
double vector_double_mul(vector_double_t *in1, vector_double_t *in2);
double vector_double_mul_part(vector_double_t *in1, vector_double_t *in2,
                              uint64_t start, uint64_t end);
void vector_double_scale(vector_double_t *in, double scale,
                         vector_double_t *out);
void vector_double_scale_add(vector_double_t *in1, vector_double_t *in2, 
                             vector_double_t *out, double scale);
void vector_double_scale_add_part(vector_double_t *in1, vector_double_t *in2, 
                                  vector_double_t *out, double scale,
                                  uint64_t start, uint64_t end);
void vector_double_copy(vector_double_t *in, vector_double_t *out);
  
#endif /* CG_VECTOR_H_ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
