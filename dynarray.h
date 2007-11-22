#ifndef __DYNARRAY_H__
#define __DYNARRAY_H__

struct dynarray;
typedef struct dynarray dynarray_t;

void *dynarray_create(unsigned long elem_size, unsigned long alloc_grain);
void *dynarray_alloc(struct dynarray *da);
void dynarray_dealloc(struct dynarray *da);
void *dynarray_alloc_nr(struct dynarray *da, unsigned long nr);
void dynarray_dealloc_nr(struct dynarray *da, unsigned long nr);
void *dynarray_get(struct dynarray *da, unsigned long idx);
void *dynarray_get_last(struct dynarray *da);
unsigned long dynarray_size(struct dynarray *da);
void *dynarray_destroy(struct dynarray *da);

#endif
