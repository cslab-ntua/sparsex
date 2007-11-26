#ifndef __SPMV_METHOD_H__
#define __SPMV_METHOD_H__

#include "method.h"
#include "vector.h"

#define SPMV_NAME(name) VECTOR_CON3(spmv_, ELEM_TYPE, name)

typedef struct {
	void *mmf_init;
	void *size;
	int elem_size;
} spmv_method_t;

typedef void spmv_double_fn_t(void *matrix, vector_double_t *in, vector_double_t *out);
typedef void spmv_float_fn_t(void *matrix, vector_float_t *in, vector_float_t *out);
typedef void *spmv_load_fn_t(char *mmf_file, 
                             unsigned long *nrows, unsigned long *ncols, 
                             unsigned long *nnz);

typedef unsigned long spmv_size_fn_t(void *matrix);

#define SPMV_M_INIT(_mmf_init, _size, _elem_size)      \
{                                               \
	.mmf_init = _mmf_init,                  \
	.size = _size,                          \
	.elem_size = _elem_size,                \
}


#define SPMV_METH_INIT(fn, _mmf_init, _size, _elem_size)               \
	static spmv_method_t _spmv_meth_ ## fn = SPMV_M_INIT(_mmf_init, _size, _elem_size); \
	METHOD_INIT(fn, &_spmv_meth_ ## fn)                            \

#endif /* __SPMV_METHOD_H__ */
