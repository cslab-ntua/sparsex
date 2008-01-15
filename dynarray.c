#include <stdlib.h>
#include <stdio.h>

#include "dynarray.h"

struct dynarray {
	void *elems;
	unsigned long elems_nr;
	unsigned long elem_size;
	unsigned long elems_idx;
	unsigned long alloc_grain;
};

void *dynarray_create(unsigned long elem_size, unsigned long alloc_grain)
{
	struct dynarray *da = malloc(sizeof(struct dynarray));
	if ( !da ){
		perror("dynarray_init: malloc");
		exit(1);
	}

	da->elems_idx = 0;
	da->elem_size = elem_size;
	da->elems_nr = da->alloc_grain = alloc_grain;

	da->elems = malloc(elem_size*da->elems_nr);
	//printf("new addr: %lu\n", (unsigned long)da->elems);
	if ( !da->elems ){
		perror("dynarray_init: malloc");
		exit(1);
	}

	return da;
}

void *dynarray_get(struct dynarray *da, unsigned long idx)
{
	unsigned long addr = (unsigned long)da->elems;
	if ( idx > da->elems_idx ){
		fprintf(stderr, "dynarray_get: out of bounds idx=%lu elems_idx=%lu\n", idx, da->elems_idx);
		exit(1);
	}
	addr += da->elem_size*idx;
	//printf("dynarray_get: idx:%lu returning %lu\n", idx, addr);
	return (void *)addr;
}

void *dynarray_get_last(struct dynarray *da)
{
	return dynarray_get(da, da->elems_idx - 1);
}

static inline void dynarray_expand(struct dynarray *da)
{
	da->elems_nr += da->alloc_grain;
	//printf("old addr: %lu  ", (unsigned long)da->elems);
	//printf("expand realloc: %lu %lu %lu\n", da->elems_idx, da->elems_nr, (da->elems_idx+1)*da->elem_size);
	da->elems = realloc(da->elems, da->elem_size*da->elems_nr);
	//printf("new addr: %lu\n", (unsigned long)da->elems);
}

void *dynarray_alloc(struct dynarray *da)
{
	void *ret;
	if ( da->elems_idx + 1 == da->elems_nr){
		dynarray_expand(da);
	}
	ret = dynarray_get(da, da->elems_idx);
	da->elems_idx++;

	return ret;
}

void dynarray_dealloc(struct dynarray *da)
{
	da->elems_idx--;
}

void *dynarray_alloc_nr(struct dynarray *da, unsigned long nr)
{
	void *ret;
	while ( da->elems_idx + nr >= da->elems_nr){
		dynarray_expand(da);
	}
	ret = dynarray_get(da, da->elems_idx);
	da->elems_idx += nr;
	
	return ret;
}

void *dynarray_alloc_nr_aligned(struct dynarray *da,
                                unsigned long nr, unsigned long align)
{
	void *ret;
	int nr_padd;

	nr_padd = (align - (da->elems_idx  & (align-1))) & (align -1);
	while (da->elems_idx + nr + nr_padd >= da->elems_nr){
		dynarray_expand(da);
	}

	da->elems_idx += nr_padd;
	ret = dynarray_get(da, da->elems_idx);
	da->elems_idx += nr;

	return ret;
}

void dynarray_dealloc_nr(struct dynarray *da, unsigned long nr)
{
	da->elems_idx -= nr;
}

unsigned long dynarray_size(struct dynarray *da)
{
	return da->elems_idx;
}

void *dynarray_destroy(struct dynarray *da)
{
	void *ret = da->elems;
	//printf("destroy realloc: %lu %lu %lu\n", da->elems_idx, da->elems_nr, (da->elems_idx+1)*da->elem_size);
	ret = realloc(ret, (da->elems_idx+1)*da->elem_size);
	free(da);
	return ret;
}
