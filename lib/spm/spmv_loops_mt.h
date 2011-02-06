#ifndef __SPMV_LOOPS_MT__
#define __SPMV_LOOPS_MT__

#include "spmv_method.h"
#include "spm_mt.h"

float spmv_double_bench_mt_loop(spm_mt_t *spm_mt, unsigned long loops,
                                unsigned long rows_nr, unsigned long cols_nr,
                                spmv_double_fn_t *fn);

float spmv_float_bench_mt_loop(spm_mt_t *spm_mt, unsigned long loops,
                               unsigned long rows_nr, unsigned long cols_nr,
                               spmv_float_fn_t *fn);

void spmv_double_check_mt_loop(void *spm, spm_mt_t *spm_mt,
                               spmv_double_fn_t *fn, unsigned long loops,
                               unsigned long rows_nr, unsigned long cols_nr,
                               spmv_double_fn_t *mt_fn);

void spmv_float_check_mt_loop(void *spm, spm_mt_t *spm_mt,
                              spmv_float_fn_t *fn, unsigned long loops,
                              unsigned long rows_nr, unsigned long cols_nr,
                              spmv_float_fn_t *mt_fn);

void spmv_double_check_mt_loop_serial(void *spm, spm_mt_t *spm_mt,
                                      spmv_double_fn_t *fn, unsigned long loops,
                                      unsigned long rows_nr, unsigned long cols_nr,
                                      spmv_double_fn_t *mt_fn);

void spmv_float_check_mt_loop_serial(void *spm, spm_mt_t *spm_mt,
                                     spmv_float_fn_t *fn, unsigned long loops,
                                     unsigned long rows_nr, unsigned long cols_nr,
                                     spmv_float_fn_t *mt_fn);
#endif

