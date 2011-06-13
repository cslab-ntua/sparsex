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

extern "C" {
#include "mt_lib.h"
#include "spm_crs_mt.h"
}
#include "cg_vector.h"

typedef struct cg_params
{
    uint64_t            nloops;
    spm_mt_thread_t     *spm_thread;
    vector_double_t     *in;
    vector_double_t     *out;
    pthread_barrier_t   *barrier;
} cg_params;

void *cg_side_thread(void *arg);

#endif /* CG_H_ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
