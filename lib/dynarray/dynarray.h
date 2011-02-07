#ifndef DYNARRAY_H__
#define DYNARRAY_H__

struct dynarray;
typedef struct dynarray dynarray_t;

#ifdef __cplusplus
extern "C" {
#endif

dynarray_t *dynarray_create(unsigned long elem_size,
                            unsigned long alloc_grain);
dynarray_t *dynarray_init_frombuff(unsigned long elem_size,
                                   unsigned long alloc_grain,
                                   void *elems, unsigned long elems_nr);

void *dynarray_alloc(struct dynarray *da);
void dynarray_dealloc(struct dynarray *da);
void *dynarray_alloc_nr(struct dynarray *da, unsigned long nr);
void dynarray_align(struct dynarray *da, unsigned long align);
void *dynarray_alloc_nr_aligned(struct dynarray *da, unsigned long nr,
                                unsigned long align);
void dynarray_dealloc(struct dynarray *da);
void dynarray_dealloc_nr(struct dynarray *da, unsigned long nr);
void dynarray_dealloc_all(struct dynarray *da);
void *dynarray_get(struct dynarray *da, unsigned long idx);
void *dynarray_get_last(struct dynarray *da);
void dynarray_seek(struct dynarray *da, unsigned long idx);
unsigned long dynarray_size(struct dynarray *da);
void *dynarray_destroy(struct dynarray *da);

#ifdef __cplusplus
}
#endif

#endif	/* DYNARRAY_H__ */
