/*
 * cg_vector.h -- Expansion of Vector Operations for CG.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
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

/**
 *  Subtraction of two vectors (C = A - B).
 *
 *  @param in1   minuend vector (A).
 *  @param in2   subtrahend vector (B).
 *  @param out   difference vector (C).
 *  @param start element from which the subtraction starts when its executed
 *               partially.
 *  @param end   element to which the subtraction ends when its executed
 *               partially.
 */
void vector_double_sub(vector_double_t *in1, vector_double_t *in2,
                       vector_double_t *out);
void vector_double_sub_part(vector_double_t *in1, vector_double_t *in2,
                            vector_double_t *out, uint64_t start, uint64_t end);
                            
/**
 *  Multiplication of two vectors (C = A * B).
 *
 *  @param in1   first input vector (A).
 *  @param in2   second input vector (B).
 *  @param out   product vector (C).
 *  @param start element from which the multiplication starts when its executed
 *               partially.
 *  @param end   element to which the multiplication ends when its executed
 *               partially.
 */
double vector_double_mul(vector_double_t *in1, vector_double_t *in2);
double vector_double_mul_part(vector_double_t *in1, vector_double_t *in2,
                              uint64_t start, uint64_t end);
                              
/**
 *  Multiply one vector with a constant (B = cA).
 *
 *  @param in    input vector (A).
 *  @param scale constant (c).
 *  @param out   output vector (B).
 */
void vector_double_scale(vector_double_t *in, double scale,
                         vector_double_t *out);
                         
/**
 *  Multiply one vector with a constant and add to another vector (C = A + cB).
 *
 *  @param in1   first input vector (A).
 *  @param in2   second input vector (B).
 *  @param scale constant (c).
 *  @param out   product vector (C).
 *  @param start element from which the operation starts when its executed
 *               partially.
 *  @param end   element to which the operation ends when its executed
 *               partially.
 */
void vector_double_scale_add(vector_double_t *in1, vector_double_t *in2, 
                             vector_double_t *out, double scale);
void vector_double_scale_add_part(vector_double_t *in1, vector_double_t *in2, 
                                  vector_double_t *out, double scale,
                                  uint64_t start, uint64_t end);
                                  
/**
 *  Copy one vector to another (B = A).
 *
 *  @param in    input vector (A).
 *  @param out   output vector (B).
 */
void vector_double_copy(vector_double_t *in, vector_double_t *out);
  
#endif /* CG_VECTOR_H_ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
