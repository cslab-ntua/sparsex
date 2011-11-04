/* 
 * cg.h -- The CG Manager Interface.
 *
 * Copyright (C) 2011,      Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef CG_H_
#define CG_H_

#include <pthread.h>
#include <sys/time.h>
#include <iostream>

extern "C" {
#include "mt_lib.h"
#include "spm_crs_mt.h"
#include "spm_crs_sym_mt.h"
#include "timer.h"
}

#include "cg_vector.h"

xtimer_t cg_timer;
xtimer_t spmv_timer;
xtimer_t reduction_timer;

typedef struct cg_params
{
    uint64_t            nloops;
    uint64_t            ncpus;
    spm_mt_thread_t     *spm_thread;
    vector_double_t     *p;
    vector_double_t     *temp;
    vector_double_t     *r;
    vector_double_t     *x;
    double              *rr;
    double              *tp;
    double              *rr_new;
    double              *ai;
    double              *bi;
    pthread_barrier_t   *barrier;
    uint64_t            start;
    uint64_t            end;    
    vector_double_t     **sub_vectors;
} cg_params;

void FindSolution(spm_mt_t *spm_mt, vector_double_t *sol, vector_double_t *b,
                  vector_double_t *temp);
void FindSymSolution(spm_mt_t *spm_mt, vector_double_t *sol,
                     vector_double_t *b, vector_double_t *temp);
                     
void InitializeCg(spm_mt_t *spm_mt, vector_double_t *x, vector_double_t *r,
                  vector_double_t *p, vector_double_t *b,
                  vector_double_t *temp);
void InitializeSymCg(spm_mt_t *spm_mt, vector_double_t *x, vector_double_t *r,
                     vector_double_t *p, vector_double_t *b,
                     vector_double_t *temp);                  

void *NormalCgSideThread(void *arg);
void *SymCgSideThread(void *arg);

void NormalCgMainThread(cg_params *params, double * cg_time, double *spmv_time,
                        double * red_time);
void SymCgMainThread(cg_params *params, double * cg_time, double *spmv_time,
                     double * red_time);

#endif /* CG_H_ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
