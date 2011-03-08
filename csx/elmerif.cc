/*
 * elmerif.cc -- Interface for Elmer integration
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include "spm.h"
#include "spmv.h"

extern "C" {
#include "elmerif.h"
#include "mt_lib.h"
#include "spmv_matvec_mt.h"
}

void elmer_matvec(void **tuned, void *n, void *rowptr, void *colind,
                  void *values, void *x, void *y, void *reinit)
{
    unsigned int nr_threads;
    unsigned int *cpus __attribute__ ((unused));
    csx::SPM *spms;
    spm_mt_t *spm_mt;
    unsigned long n_ = *((unsigned long *) n);
    uint64_t *rowptr_ = (uint64_t *) rowptr;
    uint64_t *colind_ = (uint64_t *) colind;
    double *values_ = (double *) values;
    bool reinit_ = *((bool *) reinit);
    double *x_ = (double *) x;
    double *y_ = (double *) y;

    if (!*tuned || reinit_) {
        mt_get_options(&nr_threads, &cpus);
        spms = csx::SPM::LoadFromCSR_mt(rowptr_, colind_, values_, n_, n_,
                                        nr_threads);
        *tuned = spm_mt = GetSpmMt(NULL, spms);
    } else {
        spm_mt = (spm_mt_t *) *tuned;
    }

    vector_double_t *vec_x = vector_double_create_from_buff(x_, n_);
    vector_double_t *vec_y = vector_double_create_from_buff(y_, n_);
    spmv_double_matvec_mt(spm_mt, vec_x, vec_y);
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
