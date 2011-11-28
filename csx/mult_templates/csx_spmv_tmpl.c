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

static void ctl_print(uint8_t *ctl, uint64_t start, uint64_t end,
                      const char *descr)
{
    for (uint8_t i = start; i < end; i++)
        printf("%s: ctl[%d] = %d\n", descr, i, ctl[i]);
}

typedef double (csx_spmv_fn_t)(uint8_t **ctl, uint8_t size, double **values,
                               double **x, double **y);

${spmv_func_definitions}

static csx_spmv_fn_t *mult_table[] = {
    ${spmv_func_entries}
};

void spm_csx32_double_multiply(void *spm, vector_double_t *in, vector_double_t *out)
{
    csx_double_t *csx = (csx_double_t *) spm;
    double *x = in->elements;
    double *y = out->elements;
    double *v = csx->values;
    double *x_curr = x;
    double *y_curr = y + csx->row_start;
    register double yr = 0;
    uint8_t *ctl = csx->ctl;
    uint8_t *ctl_end = ctl + csx->ctl_size;
    uint8_t size, flags;
    uint64_t i;
    uint8_t patt_id;

    for (i = 0; i < csx->nrows; i++) {
        y_curr[i] = 0;
    }

//    ctl_print(ctl, 0, csx->ctl_size, "all");

    do {
        flags = *ctl++;
        size = *ctl++;
        if (test_bit(&flags, CTL_NR_BIT)){
            *y_curr += yr;
            yr = 0;
            ${new_row_hook}
            x_curr = x;
        }
        
//        ctl_print(ctl, 0, 1, "before patt_id");
        patt_id = flags & CTL_PATTERN_MASK;
        yr += mult_table[patt_id](&ctl, size, &v, &x_curr, &y_curr);
//        ctl_print(ctl, 0, 1, "after patt_id");
    } while (ctl < ctl_end);

    *y_curr += yr;
}
