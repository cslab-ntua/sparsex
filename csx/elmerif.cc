/*
 * elmerif.cc -- Interface for Elmer integration
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <iostream>
#include "spm.h"
#include "spmv.h"
#include <sys/time.h>

extern "C" {
#include "elmerif.h"
#include "mt_lib.h"
#include "spmv_matvec_mt.h"
#include "timer.h"
}

xtimer_t load_timer, preproc_timer, spmv_timer, total_timer;

static void matvec_csr(const elmer_index_t *rowptr,
                       const elmer_index_t *colind,
                       const elmer_value_t *values,
                       elmer_index_t nr_rows,
                       const elmer_value_t *x,
                       elmer_value_t *y)
{
    elmer_index_t i, j;
    for (i = 0; i < nr_rows; ++i) {
        register elmer_value_t yi_ = 0;
        for (j = rowptr[i] - 1; j < rowptr[i+1] - 1; ++j) {
            yi_ += values[j]*x[colind[j] - 1];
        }

        y[i] = yi_;
    }
}

void elmer_matvec_(void **tuned, void *n, void *rowptr, void *colind,
                   void *values, void *u, void *v, void *reinit)
{
    unsigned int nr_threads;
    unsigned int *cpus __attribute__ ((unused));
    csx::SPM *spms;
    spm_mt_t *spm_mt;
    elmer_index_t n_ = *((elmer_index_t *) n);
    elmer_index_t *rowptr_ = (elmer_index_t *) rowptr;
    elmer_index_t *colind_ = (elmer_index_t *) colind;
    elmer_value_t *values_ = (elmer_value_t *) values;
    bool reinit_ = *((bool *) reinit);
    elmer_value_t *x_ = (elmer_value_t *) u;
    elmer_value_t *y_ = (elmer_value_t *) v;

    timer_start(&total_timer);
#ifndef CSR
    if (!*tuned || reinit_) {
        uint64_t nr_nzeros = rowptr_[n_] - 1;
        uint64_t ws = (n_ + nr_nzeros)*sizeof(elmer_index_t) +
            nr_nzeros*sizeof(elmer_value_t);
        std::cout << "Rows: " << n_ << " "
                  << "Non-zeros: " << nr_nzeros << " "
                  << "Working set size (MB): " << ((double) ws) / (1024*1024)
                  << std::endl;
        timer_start(&load_timer);
        mt_get_options(&nr_threads, &cpus);
        spms = csx::SPM::LoadCSR_mt<elmer_index_t, elmer_value_t>
            (rowptr_, colind_, values_, n_, n_, false, nr_threads);
        timer_pause(&load_timer);

        timer_start(&preproc_timer);
        *tuned = spm_mt = GetSpmMt(NULL, spms);
        timer_pause(&preproc_timer);
    } else {
        spm_mt = (spm_mt_t *) *tuned;
    }
#endif
     timer_start(&spmv_timer);
     
#ifndef CSR
    // FIXME: Although Elmer uses doubles by default, this is not portable
    vector_double_t *vec_x = vector_double_create_from_buff(x_, n_);
    vector_double_t *vec_y = vector_double_create_from_buff(y_, n_);
    vector_double_init(vec_y, 0);
    spmv_double_matvec_mt(spm_mt, vec_x, vec_y);
#else
    if (!*tuned || reinit_) {
        *tuned = rowptr_;
    }

    matvec_csr(rowptr_, colind_, values_, n_, x_, y_);
#endif
    timer_pause(&spmv_timer);
    timer_pause(&total_timer);
#ifndef CSR
    std::cout << "Load CSR: " << timer_secs(&load_timer) << " "
              << "CSX Preproc.: " << timer_secs(&preproc_timer) << " "
              << "SpMV: " << timer_secs(&spmv_timer) << " "
              << "Total: " << timer_secs(&total_timer) << std::endl;
#else
    std::cout << "SpMV: " << timer_secs(&spmv_timer) << " "
              << "Total: " << timer_secs(&total_timer) << std::endl;
#endif
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
