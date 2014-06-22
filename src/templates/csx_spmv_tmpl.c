/*
 * csx_spmv_tmpl.c -- The CSX multiplication template
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

#include <sparsex/internals/CtlUtil.hpp>
/* #include  "${header_prefix}/csx/CtlUtil.hpp" */

#define ELEM_TYPE double

#include <sparsex/internals/Vector.hpp>
#include <sparsex/internals/Csx.hpp>

/* #include "${header_prefix}/csx/Vector.hpp" */
/* #include "${header_prefix}/csx/Csx.hpp" */

#define CSX_SPMV_FN_MAX CTL_PATTERNS_MAX

#define ALIGN(buf, a) (void *) (((unsigned long) (buf) + (a-1)) & ~(a-1))

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
	for (uint64_t i = start; i < end; i++) {
		printf("%s[%ld]: %p = %d\n", descr, i, &ctl[i], ctl[i]);
        fflush(stdout);
    }
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

void spm_csx32_double_multiply(void *spm, vector_t *in, vector_t *out,
                               spx_scalar_t scale_f,
                               vector_t *local_out /* unused */)
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

    /* uint8_t *ctl_start = ctl; */
    // ctl_print(ctl, 0, csx->ctl_size, "ctl");
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
        /* printf("ctl moved at %zd bytes\n", ctl - ctl_start); */
        /* fflush(stdout); */
	} while (ctl < ctl_end);

	*y_curr += yr;
}
