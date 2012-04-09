/*
 * spmv_loops_sym_mt.h -- Multithreaded spmv loops for symmetric matrices.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
 
#ifndef __SPMV_LOOPS_SYM_MT_H__
#define __SPMV_LOOPS_SYM_MT_H__

#include <stdlib.h>
#include <pthread.h>

#include "map.h"
#include "vector.h"
#include "mt_lib.h"
#include "spm_mt.h"
#include "spmv_method.h"
#include "spm_crs_sym_mt.h"
#include "tsc.h"
#ifdef SPMV_PRFCNT
#include "prfcnt.h"
#endif /* SPMV_PRFCNT */

float spmv_double_bench_sym_mt_loop(spm_mt_t *spm_mt, unsigned long loops,
                                    unsigned long nrows, unsigned long ncols,
                                    spmv_double_fn_t *fn);

float spmv_float_bench_sym_mt_loop(spm_mt_t *spm_mt, unsigned long loops,
                                   unsigned long nrows, unsigned long ncols,
                                   spmv_float_fn_t *fn);

void spmv_double_check_sym_mt_loop(void *spm, spm_mt_t *spm_mt,
                                   spmv_double_fn_t *fn, unsigned long loops,
                                   unsigned long nrows, unsigned long ncols,
                                   spmv_double_fn_t *mt_fn);

void spmv_float_check_sym_mt_loop(void *spm, spm_mt_t *spm_mt,
                                  spmv_float_fn_t *fn, unsigned long loops,
                                  unsigned long nrows, unsigned long ncols,
                                  spmv_float_fn_t *mt_fn);

void spmv_double_check_sym_mt_loop_serial(void *spm, spm_mt_t *spm_mt,
                                          spmv_double_fn_t *fn,
                                          unsigned long loops,
                                          unsigned long nrows,
                                          unsigned long ncols,
                                          spmv_double_fn_t *mt_fn);

void spmv_float_check_sym_mt_loop_serial(void *spm, spm_mt_t *spm_mt,
                                         spmv_float_fn_t *fn,
                                         unsigned long loops,
                                         unsigned long nrows,
                                         unsigned long ncols,
                                         spmv_float_fn_t *mt_fn);
                                         		     
#endif /* __SPMV_LOOPS_SYM_MT_H__ */
