#ifndef __SPMV_LOOPS_MT__
#define __SPMV_LOOPS_MT__

#include "spmv_method.h"

float spmv_double_bench_mt_loop(spmv_double_fn_t *fn, void *matrix,
                                unsigned long loops, unsigned long cols_nr);

float spmv_float_bench_mt_loop(spmv_float_fn_t *fn, void *matrix,
                               unsigned long loops, unsigned long cols_nr);

#endif

