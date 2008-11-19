
#include <stdio.h>
#include <inttypes.h>

#include "macros.h"
#include "vector.h"
#include "jit/crsvh.h"

#define CI_TYPE UINT_TYPE(CI_BITS)

#define _NAME(name) CON5(spm_crs, CI_BITS, _vhjit_, ELEM_TYPE, name)
#define VHJIT_TYPE CON3(spm_crs, CI_BITS, _vhjit_t)

ELEM_TYPE __huff_decode_hook(unsigned char *hvals, unsigned long *bitcnt);

void _NAME(_mul_template)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	VHJIT_TYPE *vhjit = (VHJIT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	const register CI_TYPE *row_ptr = vhjit->row_ptr;
	const register CI_TYPE *col_ind = vhjit->col_ind;
	const unsigned long rs = vhjit->row_start;
	const unsigned long re = vhjit->row_end;
	register unsigned char *hvals = vhjit->hvals;
	unsigned long bitcnt=0;

	register unsigned i, j=0;
	for (i=rs; i<re; i++){
		register ELEM_TYPE yr = (ELEM_TYPE)0;
		for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
			ELEM_TYPE val = __huff_decode_hook(hvals, &bitcnt);
			yr += (val * x[col_ind[j]]);
		}
		y[i] = yr;
	}
}
