/*
 * csx_sym_spmv_tmpl.c
 *
 * \brief The CSX-Sym multiplication template
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <sparsex/types.h>
#include <sparsex/internals/CtlUtil.hpp>
#include <sparsex/internals/Vector.hpp>
#include <sparsex/internals/Csx.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

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

#if SPX_DEBUG
static void ctl_print(uint8_t *ctl, spx_index_t start, spx_index_t end,
                      const char *descr)
{
	for (spx_index_t i = start; i < end; i++)
		printf("%s: ctl[%ld] = %d\n", descr, i, ctl[i]);
}

static void deref(void *ptr)
{
	volatile unsigned long val = *((unsigned long *) ptr);
	val++;
}
#endif

${spmv_func_definitions}

void spm_csx_sym_multiply(void *spm, vector_t *in, vector_t *out,
                          spx_value_t scale_f, vector_t *local_out)
{
	csx_sym_matrix_t *csx_sym = (csx_sym_matrix_t *) spm;
	csx_matrix_t *csx = csx_sym->lower_matrix;
	spx_value_t *x = in->elements;
	spx_value_t *y = out->elements;
	spx_value_t *tmp = local_out->elements;
	spx_value_t *v = csx->values;
	spx_value_t *dv = csx_sym->dvalues;
	spx_index_t x_indx = 0;
	spx_index_t y_indx = csx->row_start;
	spx_index_t y_end = csx->row_start + csx->nrows;
    spx_index_t i;
	register spx_value_t yr = 0;
	uint8_t *ctl = csx->ctl;
	uint8_t *ctl_end = ctl + csx->ctl_size;
	uint8_t flags, size, patt_id;
	spx_value_t *cur = tmp;

	do {
		flags = *ctl++;
		size = *ctl++;
		if (test_bit(&flags, CTL_NR_BIT)) {
			y[y_indx] += yr;
			${new_row_hook}
			yr = 0;
			x_indx = 0;
			// Switch Reduction Phase
			cur = tmp;
		}
		
		${next_x}
		// Switch Reduction Phase
		if (cur != y && x_indx >= csx->row_start)
			cur = y;
		patt_id = flags & CTL_PATTERN_MASK;
		${body_hook}
	} while (ctl < ctl_end);
	
	y[y_indx] += yr;
	
	for (spx_index_t i = y_indx; i < y_end; i++) {
		y[i] += x[i] * (*dv) * scale_f;
		dv++;
	}
}
