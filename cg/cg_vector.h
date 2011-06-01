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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "vector.h"

void vector_double_add(vector_double_t *in1, vector_double_t *in2, vector_double_t *out);
void vector_double_sub(vector_double_t *in1, vector_double_t *in2, vector_double_t *out);
void vector_double_copy(vector_double_t *in, vector_double_t *out);
  
#endif /* CG_VECTOR_H_ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
