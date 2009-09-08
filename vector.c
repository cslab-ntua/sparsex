#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vector.h"

#define ALIGN_BOUND (12)
#define ALIGN(ul)  (ul + (ul % ALIGN_BOUND))


VECTOR_TYPE *VECTOR_NAME(_create)(unsigned long size)
{
	VECTOR_TYPE *v = malloc(sizeof(VECTOR_TYPE));
	if ( !v ) {
		perror("malloc");
		exit(1);
	}

	v->size = size;
	v->elements = malloc(sizeof(ELEM_TYPE)*ALIGN(size));
	if ( !v->elements ){
		perror("malloc");
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
