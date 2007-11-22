#ifndef __VECTOR_H__
#define __VECTOR_H__

/* helpers for cpp abuse */
#ifndef ELEM_TYPE
#define ELEM_TYPE double
#endif

#define VECTOR_CON3_(a,b,c) a ## b ## c
#define VECTOR_CON3(a,b,c) VECTOR_CON3_(a,b,c)
#define VECTOR_NAME(name) VECTOR_CON3(vector_, ELEM_TYPE, name)
#define VECTOR_TYPE VECTOR_NAME(_t)

typedef struct {
	float *elements;
	unsigned long size;
} vector_float_t;

vector_float_t *vector_float_create(unsigned long size);
void vector_float_destroy(vector_float_t *v);
void vector_float_init(vector_float_t *v, float val);
void vector_float_init_rand_range(vector_float_t *v, float max, float min);
int vector_float_compare(vector_float_t *v1, vector_float_t *v2);

typedef struct {
	double *elements;
	unsigned long size;
} vector_double_t;

vector_double_t *vector_double_create(unsigned long size);
void vector_double_destroy(vector_double_t *v);
void vector_double_init(vector_double_t *v, double val);
void vector_double_init_rand_range(vector_double_t *v, double max, double min);
int vector_double_compare(vector_double_t *v1, vector_double_t *v2);

#endif /* __VECTOR_H__ */
