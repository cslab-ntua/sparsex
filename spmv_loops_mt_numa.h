#ifndef SPMV_LOOPS_MT_NUMA_H__
#define SPMV_LOOPS_MT_NUMA_H__

#include "spmv_method.h"
#include "spm_mt.h"

float spmv_double_bench_mt_loop_numa(spm_mt_t *spm_mt,
                                     unsigned long loops,
                                     unsigned long rows_nr,
                                     unsigned long cols_nr,
                                     SPMV_NAME(_fn_t) *fn);

void spmv_double_check_mt_loop_numa(void *spm,
                                    spm_mt_t *spm_mt,
                                    SPMV_NAME(_fn_t) *fn,
                                    unsigned long loops,
                                    unsigned long rows_nr,
                                    unsigned long cols_nr,
                                    SPMV_NAME(_fn_t) *mt_fn);

#endif  /* SPMV_LOOPS_MT_NUMA_H__ */
