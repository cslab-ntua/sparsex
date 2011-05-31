/* -*- C -*-
 *
 * cg.h -- The CG Manager Implementation.
 *
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
    spm_crs32_double_mt_t *csr_mt = (spm_crs32_double_mt_t *) spm_thread->spm;
    vector_double_t *in = params->in;
    vector_double_t *out = params->out;
    
    ///> Set thread to the appropriate cpu.
    setaffinity_oncpu(spm_thread->cpu);
    
    ///> Do nr_loops SpMxVs.
    for (i = 0; i < nr_loops; i++) {
        pthread_barrier_wait(&barrier);
        spm_crs32_double_mt_multiply((void *) &csr_mt->crs, in, out);
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
