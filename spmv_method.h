#ifndef __SPMV_METHOD_H__
#define __SPMV_METHOD_H__

#include <inttypes.h>

#include "method.h"
#include "vector.h"

#define SPMV_NAME(name) CON3(spmv_, ELEM_TYPE, name)

typedef void spmv_double_fn_t(void *matrix, vector_double_t *in, vector_double_t *out);
typedef void spmv_float_fn_t(void *matrix, vector_float_t *in, vector_float_t *out);
typedef void *spm_load_fn_t(char *mmf, uint64_t *nrows, uint64_t *ncols, uint64_t *nnz);

typedef uint64_t spm_size_fn_t(void *matrix);
typedef void spm_destroy_fn_t(void *matrix);

typedef struct {
	spm_load_fn_t    *mmf_init_fn;
	spm_size_fn_t    *size_fn;
	spm_destroy_fn_t *destroy_fn;
	int               elem_size;
	unsigned          mt_flag:1; // multithreaded flag
} spmv_method_t;



#define __SPMV_METH_INIT(_mmf_init, _size, _destroy, _elem_size, _mt) \
{                                      \
	.mmf_init_fn  = _mmf_init,     \
	.size_fn      = _size,         \
	.destroy_fn   = _destroy,      \
	.elem_size    = _elem_size,    \
	.mt_flag      = _mt            \
}

#define _SPMV_METH_INIT(a, b, c, d)    __SPMV_METH_INIT(a, b, c, d, 0)
#define _SPMV_MT_METH_INIT(a, b, c, d) __SPMV_METH_INIT(a, b, c, d, 1)

#define SPMV_METH_INIT(fn, _mmf_init, _size, _destroy, _elem_size) \
	static spmv_method_t _spmv_meth_ ## fn = _SPMV_METH_INIT(_mmf_init, _size, _destroy, _elem_size); \
	METHOD_INIT(fn, &_spmv_meth_ ## fn) \

#define XSPMV_METH_INIT(fn, mmf_init, size, destroy, elem_size) \
         SPMV_METH_INIT(fn, mmf_init, size, destroy, elem_size)


#endif /* __SPMV_METHOD_H__ */
