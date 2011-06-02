/*
 * cg.cc -- The CG Manager Implementation.
 *
 * Copyright (C) 2011,      Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
 
#include "cg.h"

void *cg_side_thread(void *arg)
{
    uint64_t i;
    cg_params *params = (cg_params *) arg;
    uint64_t nr_loops = params->nr_loops;
    spm_mt_thread_t *spm_thread = params->spm_thread;
    spmv_double_fn_t *fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
    vector_double_t *in = params->in;
    vector_double_t *out = params->out;
    pthread_barrier_t *barrier = params->barrier;
    
    
    ///> Set thread to the appropriate cpu.
    setaffinity_oncpu(spm_thread->cpu);
    
    ///> Do nr_loops SpMxVs.
    for (i = 0; i < nr_loops; i++) {
        pthread_barrier_wait(barrier);
        fn(spm_thread->spm, in, out);
        pthread_barrier_wait(barrier);
    }
    return NULL;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
