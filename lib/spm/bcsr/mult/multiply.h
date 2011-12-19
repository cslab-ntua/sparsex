/*
 * bcsr/multiply.h -- BCSR block-specific multiplication routines
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __BCSR_MULTIPLY_H__
#define __BCSR_MULTIPLY_H__

#include "spm_bcsr.h"
#include "spm_bcsr_mt.h"

#define DECLARE_BCSR_METH(r,c)                  \
    void SPM_BCSR_NAME(_multiply_ ## r ## x ## c)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out); \
    void SPM_BCSR_MT_NAME(_multiply_ ## r ## x ## c)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out);

/* generic multiplication methods */
void SPM_BCSR_NAME(_multiply_generic)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out);
void SPM_BCSR_MT_NAME(_multiply_generic)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out);

DECLARE_BCSR_METH(1,2);
DECLARE_BCSR_METH(1,3);
DECLARE_BCSR_METH(1,4);
DECLARE_BCSR_METH(1,5);
DECLARE_BCSR_METH(1,6);
DECLARE_BCSR_METH(1,7);
DECLARE_BCSR_METH(1,8);
DECLARE_BCSR_METH(2,1);
DECLARE_BCSR_METH(2,2);
DECLARE_BCSR_METH(2,3);
DECLARE_BCSR_METH(2,4);
DECLARE_BCSR_METH(3,1);
DECLARE_BCSR_METH(3,2);
DECLARE_BCSR_METH(3,3);
DECLARE_BCSR_METH(4,1);
DECLARE_BCSR_METH(4,2);
DECLARE_BCSR_METH(5,1);
DECLARE_BCSR_METH(6,1);
DECLARE_BCSR_METH(7,1);
DECLARE_BCSR_METH(8,1);

#define BCSR_MULT_ROW_MAX 8
#define BCSR_MULT_COL_MAX 8

// table of single-thread functions
static spmv_double_fn_t *SPM_BCSR_NAME(_mult_table)[][BCSR_MULT_COL_MAX] = {
    /* 1x? */
    {
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_1x2),
        SPM_BCSR_NAME(_multiply_1x3),
        SPM_BCSR_NAME(_multiply_1x4),
        SPM_BCSR_NAME(_multiply_1x5),
        SPM_BCSR_NAME(_multiply_1x6),
        SPM_BCSR_NAME(_multiply_1x7),
        SPM_BCSR_NAME(_multiply_1x8),
    },

    /* 2x? */
    {
        SPM_BCSR_NAME(_multiply_2x1),
        SPM_BCSR_NAME(_multiply_2x2),
        SPM_BCSR_NAME(_multiply_2x3),
        SPM_BCSR_NAME(_multiply_2x4),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
    },

    /* 3x? */
    {
        SPM_BCSR_NAME(_multiply_3x1),
        SPM_BCSR_NAME(_multiply_3x2),
        SPM_BCSR_NAME(_multiply_3x3),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
    },

    /* 4x? */
    {
        SPM_BCSR_NAME(_multiply_4x1),
        SPM_BCSR_NAME(_multiply_4x2),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
    },

    /* 5x? */
    {
        SPM_BCSR_NAME(_multiply_5x1),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
    },

    /* 6x? */
    {
        SPM_BCSR_NAME(_multiply_6x1),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
    },

    /* 7x? */
    {
        SPM_BCSR_NAME(_multiply_7x1),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
    },

    /* 8x? */
    {
        SPM_BCSR_NAME(_multiply_8x1),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
        SPM_BCSR_NAME(_multiply_generic),
    },
};

// Table of multithreaded functions
static spmv_double_fn_t *SPM_BCSR_MT_NAME(_mult_table)[][BCSR_MULT_COL_MAX] = {
    /* 1x? */
    {
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_1x2),
        SPM_BCSR_MT_NAME(_multiply_1x3),
        SPM_BCSR_MT_NAME(_multiply_1x4),
        SPM_BCSR_MT_NAME(_multiply_1x5),
        SPM_BCSR_MT_NAME(_multiply_1x6),
        SPM_BCSR_MT_NAME(_multiply_1x7),
        SPM_BCSR_MT_NAME(_multiply_1x8),
    },

    /* 2x? */
    {
        SPM_BCSR_MT_NAME(_multiply_2x1),
        SPM_BCSR_MT_NAME(_multiply_2x2),
        SPM_BCSR_MT_NAME(_multiply_2x3),
        SPM_BCSR_MT_NAME(_multiply_2x4),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
    },

    /* 3x? */
    {
        SPM_BCSR_MT_NAME(_multiply_3x1),
        SPM_BCSR_MT_NAME(_multiply_3x2),
        SPM_BCSR_MT_NAME(_multiply_3x3),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
    },

    /* 4x? */
    {
        SPM_BCSR_MT_NAME(_multiply_4x1),
        SPM_BCSR_MT_NAME(_multiply_4x2),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
    },

    /* 5x? */
    {
        SPM_BCSR_MT_NAME(_multiply_5x1),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
    },

    /* 6x? */
    {
        SPM_BCSR_MT_NAME(_multiply_6x1),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
    },

    /* 7x? */
    {
        SPM_BCSR_MT_NAME(_multiply_7x1),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
    },

    /* 8x? */
    {
        SPM_BCSR_MT_NAME(_multiply_8x1),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
        SPM_BCSR_MT_NAME(_multiply_generic),
    },
};

static inline spmv_double_fn_t *SPM_BCSR_NAME(_mult_method)(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c)
{
    if (r > BCSR_MULT_ROW_MAX || c > BCSR_MULT_COL_MAX)
        return SPM_BCSR_NAME(_multiply_generic);
    else
        return SPM_BCSR_NAME(_mult_table)[r-1][c-1];
}


static inline spmv_double_fn_t *SPM_BCSR_MT_NAME(_mult_method)(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c)
{
    if (r > BCSR_MULT_ROW_MAX || c > BCSR_MULT_COL_MAX)
        return SPM_BCSR_MT_NAME(_multiply_generic);
    else
        return SPM_BCSR_MT_NAME(_mult_table)[r-1][c-1];
}

#endif /* __BCSR_MULTIPLY_H__ */
