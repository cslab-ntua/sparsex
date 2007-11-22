#include <stdlib.h>
#include <stdio.h>

#include "beanbag.h"

struct beanbag {
	void *beans;
	unsigned long beans_nr;
	unsigned long bean_size;
	unsigned long beans_idx;
	unsigned long beans_grain;
};

void *beanbag_create(unsigned long bean_size, unsigned long beans_grain)
{
	struct beanbag *bb = malloc(sizeof(struct beanbag));
	if ( !bb ){
		perror("beanbag_init: malloc");
		exit(1);
	}

	bb->beans_idx = 0;
	bb->bean_size = bean_size;
	bb->beans_nr = bb->beans_grain = beans_grain;

	bb->beans = malloc(bean_size*bb->beans_nr);
	//printf("new addr: %lu\n", (unsigned long)bb->beans);
	if ( !bb->beans ){
		perror("beanbag_init: malloc");
		exit(1);
	}

	return bb;
}

void *beanbag_get(struct beanbag *bb, unsigned long idx)
{
	unsigned long addr = (unsigned long)bb->beans;
	if ( idx > bb->beans_idx ){
		fprintf(stderr, "beanbag_get: out of bounds idx=%lu beans_idx=%lu\n", idx, bb->beans_idx);
		exit(1);
	}
	addr += bb->bean_size*idx;
	//printf("beanbag_get: idx:%lu returning %lu\n", idx, addr);
	return (void *)addr;
}

void *beanbag_get_last(struct beanbag *bb)
{
	return beanbag_get(bb, bb->beans_idx - 1);
}

static inline void beanbag_expand(struct beanbag *bb)
{
	bb->beans_nr += bb->beans_grain;
	//printf("old addr: %lu  ", (unsigned long)bb->beans);
	//printf("expand realloc: %lu %lu %lu\n", bb->beans_idx, bb->beans_nr, (bb->beans_idx+1)*bb->bean_size);
	bb->beans = realloc(bb->beans, bb->bean_size*bb->beans_nr);
	//printf("new addr: %lu\n", (unsigned long)bb->beans);
}

void *beanbag_alloc(struct beanbag *bb)
{
	void *ret;
	if ( bb->beans_idx + 1 == bb->beans_nr){
		beanbag_expand(bb);
	}
	ret = beanbag_get(bb, bb->beans_idx);
	bb->beans_idx++;

	return ret;
}

void *beanbag_alloc_nr(struct beanbag *bb, unsigned long nr)
{
	void *ret;
	while ( bb->beans_idx + nr >= bb->beans_nr){
		beanbag_expand(bb);
	}
	ret = beanbag_get(bb, bb->beans_idx);
	bb->beans_idx += nr;
	
	return ret;
}

unsigned long beanbag_size(struct beanbag *bb)
{
	return bb->beans_idx;
}

void *beanbag_destroy(struct beanbag *bb)
{
	void *ret = bb->beans;
	//printf("destroy realloc: %lu %lu %lu\n", bb->beans_idx, bb->beans_nr, (bb->beans_idx+1)*bb->bean_size);
	ret = realloc(ret, (bb->beans_idx+1)*bb->bean_size);
	free(bb);
	return ret;
}
