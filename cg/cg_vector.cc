/*
 * cg_vector.cc -- Implementation of Expansion of Vector Operations for CG.
 *
 * Copyright (C) 2011,      Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "cg_vector.h"

void vector_double_add(vector_double_t *in1, vector_double_t *in2, vector_double_t *out)
{
    unsigned long i;
    
    assert(in1->size == in2->size && in1->size == out->size && "vectors for add have different size");
    
    for (i = 0; i < in1->size; i++)
        out->elements[i] = in1->elements[i] + in2->elements[i];
}

void vector_double_sub(vector_double_t *in1, vector_double_t *in2, vector_double_t *out)
{
    unsigned long i;
    
    assert(in1->size == in2->size && in1->size == out->size && "vectors for sub have different size");
    
    for (i = 0; i < in1->size; i++)
        out->elements[i] = in1->elements[i] - in2->elements[i];
}

double vector_double_mul(vector_double_t *in1, vector_double_t *in2)
{
    unsigned long i;
    double out;

    assert(in1->size == in2->size &&  "vectors for mul have different size");

    out = 0;
    for (i = 0; i < in1->size; i++)
        out += in1->elements[i] * in2->elements[i];
    return out;
}

void vector_double_scale(vector_double_t *in, double scale, vector_double_t *out)
{
    unsigned long i;

    assert(in->size == out->size &&  "vectors for scale have different size");

    for (i = 0; i < in->size; i++)
        out->elements[i] = in->elements[i] * scale;
}

void vector_double_copy(vector_double_t *in, vector_double_t *out)
{
    unsigned long i;
    
    assert(in->size == out->size && "vectors for copy have different size");
    
    for (i = 0; i < in->size; i++)
        out->elements[i] = in->elements[i];
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
