/*
 * csx_spmv_tmpl.c -- The CSX multiplication template
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

#include  "${header_prefix}/csx/ctl_ll.h"

#define ELEM_TYPE double
#include "${header_prefix}/lib/spm/vector.h"
#include "${header_prefix}/csx/csx.h"

#define CSX_SPMV_FN_MAX CTL_PATTERNS_MAX

#define ALIGN(buf,a) (void *) (((unsigned long) (buf) + (a-1)) & ~(a-1))

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
// the following function will be unused when only delta8 units
// are present
static void align_ptr(uint8_t **ctl, int align)
{
	*ctl = ALIGN(*ctl, align);
}
#pragma GCC diagnostic pop

#ifdef CSX_DEBUG
static void ctl_print(uint8_t *ctl, uint64_t start, uint64_t end,
                      const char *descr)
{
	for (uint64_t i = start; i < end; i++)
		printf("%s: ctl[%ld] = %d\n", descr, i, ctl[i]);
}

static void deref(void *ptr)
{
	volatile unsigned long val = *((unsigned long *) ptr);
	val++;
}
#endif

typedef double (csx_spmv_fn_t)(uint8_t **ctl, uint8_t size, double **values,
                               double **x, double **y);

${spmv_func_definitions}

void spm_csx32_double_multiply(void *spm, vector_double_t *in,
                               vector_double_t *out)
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
	uint8_t patt_id;

	for (uint64_t i = 0; i < csx->nrows; i++)
		y_curr[i] = 0;

	do {
		flags = *ctl++;
		size = *ctl++;
		if (test_bit(&flags, CTL_NR_BIT)) {
			*y_curr += yr;
			yr = 0;
			${new_row_hook}
			x_curr = x;
		}

		${next_x}
		patt_id = flags & CTL_PATTERN_MASK;
		${body_hook}
	} while (ctl < ctl_end);

	*y_curr += yr;
}
