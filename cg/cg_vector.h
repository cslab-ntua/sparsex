/* -*- C -*-
 *
 * cg_vector.h -- Interface of Expansion of Vector Operations for CG.
 *
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
 
#ifndef CG_VECTOR_H_
#define CG_VECTOR_H_

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <numa.h>
#include <numaif.h>
#include <sys/mman.h>

#include "vector.h"

void vector_double_add(vector_double_t *in1, vector_double_t *in2, vector_double_t *out);
void vector_double_sub(vector_double_t *in1, vector_double_t *in2, vector_double_t *out);
void vector_double_mul(vector_double_t *in, double mul, vector_double_t *out);
void vector_double_div(vector_double_t *in1, vector_double_t *in2, vector_double_t *out);
double vector_double_max(vector_double_t *in);
double vector_double_min(vector_double_t *in);
double vector_double_avg(vector_double_t *in);
void vector_double_copy(vector_double_t *in, vector_double_t *out);
void vector_double_abs(vector_double_t *in);
  
#endif /* CG_VECTOR_H_ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
