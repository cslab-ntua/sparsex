#include <assert.h>
#include <numa.h>
#include <sched.h>
#include <stdlib.h>
#include <stdio.h>

#include "dynarray.h"

struct dynarray {
	void *elems;
	unsigned long elems_nr;
	unsigned long elem_size;
	unsigned long next_idx;
	unsigned long alloc_grain;
};

dynarray_t *dynarray_create(unsigned long elem_size, unsigned long alloc_grain)
{
    struct dynarray *da;
    da = malloc(sizeof(*da));

	if ( !da ) {
		perror("dynarray_create: malloc");
		exit(1);
	}

	da->next_idx = 0;
	da->elem_size = elem_size;
	da->elems_nr = da->alloc_grain = alloc_grain;

    da->elems = malloc(elem_size*da->elems_nr);

	//printf("da->elems: %p\n", da->elems);
	if ( !da->elems ){
		perror("dynarray_create: malloc");
		exit(1);
	}

	return da;
}

/*
 * elems pointer should be subjectible to realloc() calls
 */
dynarray_t *dynarray_init_frombuff(unsigned long elem_size,
								   unsigned long alloc_grain,
								   void *elems, unsigned long elems_nr)
{
	struct dynarray *da = malloc(sizeof(struct dynarray));
	if (!da){
		perror("dynarray_init_frombuff: malloc");
		exit(1);
	}

	da->next_idx = elems_nr;
	da->elem_size = elem_size;
	da->alloc_grain = alloc_grain;
	da->elems = elems;
	da->elems_nr = elems_nr;

	// fix buffer to be aligned with alloc_grain
	unsigned long rem = elems_nr % alloc_grain;
	if (rem != 0) {
		da->elems_nr += alloc_grain - rem;
		da->elems = realloc(da->elems, da->elem_size*da->elems_nr);
		if (!da->elems) {
			perror("dynarray_init_from_buff: realloc failed");
			exit(1);
		}

		// In case da->elems is relocated.
		elems = da->elems;
	}

	return da;
}

void dynarray_seek(struct dynarray *da, unsigned long idx)
{
	if (idx >= da->next_idx) {
		fprintf(stderr, "dynarray_seek: out of bounds idx=%lu next_idx=%lu\n",
				idx, da->next_idx);
		exit(1);
	}

	da->next_idx = idx;
}

void *dynarray_get(struct dynarray *da, unsigned long idx)
{
	unsigned long addr = (unsigned long) da->elems;
	if (idx >= da->next_idx) {
		fprintf(stderr, "dynarray_get: out of bounds idx=%lu next_idx=%lu\n",
                idx, da->next_idx);
		//print_trace();
		exit(1);
	}
	addr += da->elem_size*idx;
	//printf("dynarray_get: addr:%p next_idx:%lu idx:%lu returning 0x%lx\n", da->elems, da->next_idx, idx, addr);
	return (void *)addr;
}

void *dynarray_get_last(struct dynarray *da)
{
	return dynarray_get(da, da->next_idx-1);
}

static inline void dynarray_expand(struct dynarray *da)
{
	da->elems_nr += da->alloc_grain;
	//printf("old addr: %lu	 ", (unsigned long)da->elems);
	//printf("expand realloc: %lu %lu %lu\n", da->next_idx, da->elems_nr, (da->next_idx+1)*da->elem_size);
    da->elems = realloc(da->elems, da->elem_size*da->elems_nr);
        
	if (!da->elems) {
		perror("dynarray_expand: realloc failed");
		exit(1);
	}
	//printf("new addr: %lu\n", (unsigned long)da->elems);
}

void *dynarray_alloc(struct dynarray *da)
{
	void *ret;
	if (da->next_idx >= da->elems_nr) {
		dynarray_expand(da);
	}

	da->next_idx++;
	ret = dynarray_get_last(da);

	return ret;
}

void dynarray_dealloc(struct dynarray *da)
{
	da->next_idx--;
}

void *dynarray_alloc_nr(struct dynarray *da, unsigned long nr)
{
	void *ret;
	while (da->next_idx + nr >= da->elems_nr) {
		dynarray_expand(da);
	}

	unsigned long idx = da->next_idx;
	da->next_idx += nr;
	ret = dynarray_get(da, idx);

	return ret;
}

void dynarray_align(struct dynarray *da, unsigned long align)
{
	int nr_padd = (align - (da->next_idx  & (align-1))) & (align -1);
	while (da->next_idx + nr_padd >= da->elems_nr) {
		dynarray_expand(da);
	}

	da->next_idx += nr_padd;
}

void *dynarray_alloc_nr_aligned(struct dynarray *da,
								unsigned long nr, unsigned long align)
{
	dynarray_align(da, align);
	return dynarray_alloc_nr(da, nr);
}

void dynarray_dealloc_nr(struct dynarray *da, unsigned long nr)
{
	if (da->next_idx < nr){
		fprintf(stderr, "dynarray_dealloc_nr: %lu %lu\n", da->next_idx, nr);
		exit(1);
	}
	da->next_idx -= nr;
}

void dynarray_dealloc_all(struct dynarray *da)
{
	da->next_idx = 0;
}

unsigned long dynarray_size(struct dynarray *da)
{
	return da->next_idx;
}

void *dynarray_destroy(struct dynarray *da)
{
	void *ret = da->elems;
	//printf("destroy realloc: idx:%lu nr:%lu realloc size:%lu\n", da->next_idx, da->elems_nr, (da->next_idx+1)*da->elem_size);
    ret = realloc(ret, da->next_idx*da->elem_size);
    free(da);
	return ret;
}
