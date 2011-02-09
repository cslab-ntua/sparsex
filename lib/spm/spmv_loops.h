/*
 * spmv_loops.h
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __SPMV_LOOPS_H__
#define __SPMV_LOOPS_H__

#include "spmv_method.h"

float spmv_double_bench_loop(spmv_double_fn_t *fn, void *matrix,
                             unsigned long loops,
                             unsigned long rows_nr, unsigned long cols_nr);

float spmv_float_bench_loop(spmv_float_fn_t *fn, void *matrix,
                            unsigned long loops,
                            unsigned long rows_nr, unsigned long cols_nr);

void spmv_double_check_loop(void *m1, void *m2,
                            spmv_double_fn_t *fn1, spmv_double_fn_t *fn2,
                            unsigned long loops,
                            unsigned long rows_nr, unsigned long cols_nr,
                            unsigned long nz_nr);

void spmv_float_check_loop(void *m1, void *m2,
                           spmv_float_fn_t *fn1, spmv_float_fn_t *fn2,
                           unsigned long loops,
                           unsigned long rows_nr, unsigned long cols_nr,
                           unsigned long nz_nr);
#endif
