#ifndef __SPMV_METHOD_H__
#define __SPMV_METHOD_H__

#include "method.h"
#include "vector.h"

#define SPMV_NAME(name) VECTOR_CON3(spmv_, ELEM_TYPE, name)

typedef struct {
	void *mmf_init;
	void *spmv_loop;
} spmv_method_t;

typedef void spmv_double_fn_t(void *matrix, vector_double_t *in, vector_double_t *out);
typedef void spmv_float_fn_t(void *matrix, vector_float_t *in, vector_float_t *out);
typedef void *spmv_load_fn_t(char *mmf_file, 
                             unsigned long *nrows, unsigned long *ncols, 
                             unsigned long *nnz);

#define SPMV_METH_INIT(function, _mmf_init, _spmv_loop)  \
	METHOD_INIT(function,                            \
		    ( (spmv_method_t *) {                \
			(void *)_mmf_init,             \
			(void *)_spmv_loop            \
		    })                                   \
	)

#endif /* __SPMV_METHOD_H__ */
