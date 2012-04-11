/* 
 * cg.h -- The CG Manager Interface.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
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

///> Timers for CG breakdown to three different parts.
xtimer_t cg_timer;
xtimer_t spmv_timer;
xtimer_t reduction_timer;

///> Parameters needed for multithreaded CG execution.
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

/**
 *  Find b vector for pre-determined solution in the equation (Ax = b).
 *
 *  @param spm_mt symmetric sparse matrix in the appropriate format.
 *  @param sol    pre-determined solution x (used to test performance of CG).
 *  @param b      output vector b.
 *  @param temp   helpful vector for this operation in order to avoid redudant
 *                allocations.
 */
void FindSolution(spm_mt_t *spm_mt, vector_double_t *sol, vector_double_t *b,
                  vector_double_t *temp);
void FindSymSolution(spm_mt_t *spm_mt, vector_double_t *sol,
                     vector_double_t *b, vector_double_t *temp);

/**
 *  Initialize parameters of CG.
 *
 *  @param spm_mt symmetric sparse in the appropriate format.
 *  @param x      initial estimated solution x.
 *  @param r      initial residual.
 *  @param p      vector p used in CG algorithm. Initially is equal to residual.
 *  @param b      the product in the equation Αχ = b.
 *  @param temp   helpful vector for the initialization in order to avoid
 *                redudant allocations.
 */                     
void InitializeCg(spm_mt_t *spm_mt, vector_double_t *x, vector_double_t *r,
                  vector_double_t *p, vector_double_t *b,
                  vector_double_t *temp);
void InitializeSymCg(spm_mt_t *spm_mt, vector_double_t *x, vector_double_t *r,
                     vector_double_t *p, vector_double_t *b,
                     vector_double_t *temp);                  

/**
 *  CG algorithm for regular and symmetric sparse matrix formats.
 *
 *  @param arg       parameters needed for multithreaded execution of CG.
 *  @param params    parameters needed for multithreaded execution of CG.
 *  @param cg_time   total execution time for CG.
 *  @param spmv_time SpMxV's multiplication operation time for CG.
 *  @param red_time  SpMxV's reduction operation time for CG.
 */
void *NormalCgSideThread(void *arg);
void *SymCgSideThread(void *arg);
void NormalCgMainThread(cg_params *params, double * cg_time, double *spmv_time,
                        double * red_time);
void SymCgMainThread(cg_params *params, double * cg_time, double *spmv_time,
                     double * red_time);

#endif /* CG_H_ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
