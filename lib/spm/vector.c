/*
 * vector.c
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vector.h"

/*
 * vector: align allocation to 12(=3x4) elements.
 * This (quick-)fixes a bug in bcrs spmv code, when the row block size
 * does not align with the y vector. As a results writes are
 * perfomed beyond the end of the y vector
 */
#define ALIGN_BOUND (12)
#define ALIGN(ul)  (ul + (ul % ALIGN_BOUND))

VECTOR_TYPE *VECTOR_NAME(_create)(unsigned long size)
{
	VECTOR_TYPE *v = malloc(sizeof(VECTOR_TYPE));
	if ( !v ) {
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	v->size = size;
	v->elements = malloc(sizeof(ELEM_TYPE)*(size + 12));
	if ( !v->elements ){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	return v;
}

void VECTOR_NAME(_destroy)(VECTOR_TYPE *v)
{
	free(v->elements);
	free(v);
}

void VECTOR_NAME(_init)(VECTOR_TYPE *v, ELEM_TYPE val)
{
	unsigned long i;
	for (i=0 ; i<v->size; i++){
		v->elements[i] = val;
	}
}

void VECTOR_NAME(_init_rand_range)(VECTOR_TYPE *v, ELEM_TYPE max, ELEM_TYPE min)
{
	unsigned long i;
	ELEM_TYPE val;
	for (i=0; i<v->size; i++){
		val = ((ELEM_TYPE)(rand()+i)/((ELEM_TYPE)RAND_MAX+1));
		v->elements[i] = min + val*(max-min);
	}
}

static inline int elems_neq(ELEM_TYPE a, ELEM_TYPE b)
{
	if ( fabs((double)(a-b)/(double)a)  > 0.0000001 ){
		return 1;
	}
	return 0;
}

int VECTOR_NAME(_compare)(VECTOR_TYPE *v1, VECTOR_TYPE *v2)
{
	unsigned long i;
	if ( v1->size != v2->size ){
		fprintf(stderr, "v1->size=%lu v2->size=%lu differ", v1->size, v2->size);
		return -2;
	}
	for (i=0; i<v1->size; i++){
		if ( elems_neq(v1->elements[i], v2->elements[i])){
			fprintf(stderr, "element %ld differs: %10.20lf != %10.20lf\n",
			        i, (double)v1->elements[i], (double)v2->elements[i]);
			return -1;
		}
	}

	return 0;
}
