/*
 * csx_spmv_tmpl.c -- The CSX multiplication template
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>

#include "ctl_ll.h"

#define ELEM_TYPE double
#include "vector.h"
#include "csx.h"

#define CSX_SPMV_FN_MAX CTL_PATTERNS_MAX

typedef void (csx_spmv_fn_t)(uint8_t *ctl, uint8_t size, double *values,
                             double *x, double *y);

${spmv_func_definitions}

static csx_spmv_fn_t *mult_table[] = {
    ${spmv_func_entries}
};

void csx_spmv_template(void *spm, vector_double_t *in, vector_double_t *out)
{
    csx_double_t *csx = (csx_double_t *) spm;
    double *x;
    double *y;
    double *v = csx->values;
    double *myx;
    register double yr = 0;
    uint8_t *ctl = csx->ctl;
    uint8_t *ctl_end = ctl + csx->ctl_size;
    uint64_t y_indx = csx->row_start;
    uint8_t size, flags;
    uint64_t i;
    uint8_t patt_id;
    //printf("csx->ctl: %p\n", csx->ctl);
    x = in  ? in->elements : NULL;
    y = out ? out->elements: NULL;
    myx = x;
    for (i = csx->row_start; i < csx->row_start + csx->nrows; i++)
        y[i] = 0;

    do {
        //printf("ctl:%p\n", ctl);
        flags = *ctl++;
        size = *ctl++;
        //printf("size=%d\n", size);
        if (test_bit(&flags, CTL_NR_BIT)){
            ${new_row_hook}
            myx = x;
            yr = 0;
            //y[y_indx] = yr;
        }
        //printf("x_indx before jmp: %lu\n", myx - x);
        myx += ul_get(&ctl);
        //printf("x_indx after jmp: %lu\n", myx - x);
        patt_id = flags & CTL_PATTERN_MASK;
        mult_table[patt_id](ctl, size, v, x, y);
        //printf("x_indx at end: %lu\n", myx - x);
    } while (ctl < ctl_end);

    y[y_indx] += yr;
}
