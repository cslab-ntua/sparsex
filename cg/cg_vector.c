/* -*- C -*-
 *
 * cg_vector.c -- Implementation of Expansion of Vector Operations for CG.
 *
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "cg_vector.h"

void vector_double_add(vector_double_t *in1, vector_double_t *in2, vector_double_t *out)
{
    unsigned long int i;
    
    assert(in1->size == in2->size && in1->size == out->size && "vectors for add have different size");
    
    for (i = 0; i < in1->size; i++)
        out->elements[i] = in1->elements[i] + in2->elements[i];
}

void vector_double_sub(vector_double_t *in1, vector_double_t *in2, vector_double_t *out)
{
    unsigned long int i;
    
    assert(in1->size == in2->size && in1->size == out->size && "vectors for sub have different size");
    
    for (i = 0; i < in1->size; i++)
        out->elements[i] = in1->elements[i] - in2->elements[i];
}

void vector_double_mul(vector_double_t *in, double mul, vector_double_t *out)
{
    unsigned long int i;
    
    assert(in->size == out->size && "vectors for mul have different size");
    
    for (i = 0; i < in->size; i++)
        out->elements[i] = mul*in->elements[i];
}

void vector_double_div(vector_double_t *in1, vector_double_t *in2, vector_double_t *out)
{
    unsigned long int i;
    
    assert(in1->size == in2->size && in1->size == out->size && "vectors for div have different size");
    
    for (i = 0; i < in1->size; i++)
        out->elements[i] = in1->elements[i]/in2->elements[i];
}

double vector_double_max(vector_double_t *in)
{
    unsigned long int i;
    double max = in->elements[0];
    
    for (i=1; i < in->size; i++)
        if (in->elements[i] > max)
            max = in->elements[i];
    return max;
}

double vector_double_min(vector_double_t *in)
{
    unsigned long int i;
    double min = in->elements[0];
    
    for (i=1; i < in->size; i++)
        if (in->elements[i] < min)
            min = in->elements[i];
    return min;
}

double vector_double_avg(vector_double_t *in)
{
    unsigned long int i;
    double sum = 0.0;
    
    for (i=0; i < in->size; i++)
        sum += in->elements[i];
    return sum/in->size;
}

void vector_double_copy(vector_double_t *in, vector_double_t *out)
{
    unsigned long int i;
    
    assert(in->size == out->size && "vectors for copy have different size");
    
    for (i = 0; i < in->size; i++)
        out->elements[i] = in->elements[i];
}


void vector_double_abs(vector_double_t *in)
{
    unsigned long int i;
    
    for (i=0; i < in->size; i++)
        if (in->elements[i] < 0)
            in->elements[i] = 0-in->elements[i];
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
