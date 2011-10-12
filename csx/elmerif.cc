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
#include <math.h>
#include <cfloat>
#include "spm.h"
#include "spmv.h"

extern "C" {
#include "elmerif.h"
#include "mt_lib.h"
#include "spmv_matvec_mt.h"
#include "spm_crs_mt.h"
#include "timer.h"
}

enum {
    TIMER_CONSTRUCT_INTERN = 0,
    TIMER_CONSTRUCT_CSX,
    TIMER_SPMV,
    TIMER_TOTAL,
    TIMER_END,
};

static uint64_t nr_calls = 0;
xtimer_t timers[TIMER_END];

const char *timer_desc[] = {
    "Convert to internal repr.",
    "Convert to CSX",
    "Multiplication",
    "Total time in library",
};

void __attribute__ ((constructor)) init() {
    for (int i = 0; i < TIMER_END; ++i) {
        timer_init(&timers[i]);
    }
}

void __attribute__ ((destructor)) finalize() {
    std::cout << "libcsx called " << nr_calls
              << " time(s)" << std::endl;

    // Print timer statistics
    for (int i = 0; i < TIMER_END; ++i) {
        std::cout << timer_desc[i] << ": "
                  << timer_secs(&timers[i]) << std::endl;
    }
}

#ifdef _DEBUG_
elmer_value_t *y_copy = 0;
elmer_value_t *values_last = 0;

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

static bool equal(elmer_value_t *v1, elmer_value_t *v2,
                    elmer_index_t n)
{
    for (elmer_index_t i = 0; i < n; ++i) {
        if (fabs(v1[i] - v2[i]) > DBL_EPSILON) {
            std::cout << "Element " << i << " differs: "
                      << v1[i] << " != " << v2[i] << std::endl;
            return false;
        }
    }

    return true;
}

#endif

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

    ++nr_calls;
    timer_start(&timers[TIMER_TOTAL]);
#ifdef _DEBUG_
    if (!y_copy)
        y_copy = (elmer_value_t *) malloc(n_*sizeof(*y_copy));
    memcpy(y_copy, y_, n_*sizeof(*y_copy));
    
#if 0
    size_t nnz = rowptr_[n_] - 1;
    if (!values_last) {
        values_last = (elmer_value_t *) malloc(nnz*sizeof(*values_last));
    } else {
        std::cout << "Checking values ...";
        if (!equal(values_, values_last, nnz)) {
            std::cout << "FAILED" << std::endl;
        } else {
            std::cout << "DONE" << std::endl;
        }
    }

    memcpy(values_last, values_, nnz*sizeof(*values_last));
#endif
    matvec_csr(rowptr_, colind_, values_, n_, x_, y_copy);
#endif

    if (!*tuned || reinit_) {
        uint64_t nr_nzeros = rowptr_[n_] - 1;
        uint64_t ws = (n_ + nr_nzeros)*sizeof(elmer_index_t) +
            nr_nzeros*sizeof(elmer_value_t);
        std::cout << "Rows: " << n_ << " "
                  << "Non-zeros: " << nr_nzeros << " "
                  << "Working set size (MB): " << ((double) ws) / (1024*1024)
                  << std::endl;
        mt_get_options(&nr_threads, &cpus);
#ifdef _CSR_
        spm_mt = spm_crs32_double_mt_get_spm((uint32_t *) rowptr_,
                                             (uint32_t *) colind_,
                                             values_, n_, nr_threads, cpus);
#else
        timer_start(&timers[TIMER_CONSTRUCT_INTERN]);
        spms = csx::SPM::LoadCSR_mt<elmer_index_t, elmer_value_t>
            (rowptr_, colind_, values_, n_, n_, false, nr_threads);
        timer_pause(&timers[TIMER_CONSTRUCT_INTERN]);

        timer_start(&timers[TIMER_CONSTRUCT_CSX]);
        spm_mt = GetSpmMt(NULL, spms);
        timer_pause(&timers[TIMER_CONSTRUCT_CSX]);
#endif        
        *tuned = spm_mt;
    } else {
        spm_mt = (spm_mt_t *) *tuned;
    }
     
    // FIXME: Although Elmer uses doubles by default, this is not portable
    vector_double_t *vec_x = vector_double_create_from_buff(x_, n_);
    vector_double_t *vec_y = vector_double_create_from_buff(y_, n_);
    timer_start(&timers[TIMER_SPMV]);
    spmv_double_matvec_mt(spm_mt, vec_x, vec_y);
    timer_pause(&timers[TIMER_SPMV]);

#ifdef _DEBUG_
    std::cout << "Checking result... ";
    if (equal(vec_y->elements, y_copy, n_))
        std::cout << "DONE" << std::endl;
    else
        std::cout << "FAILED" << std::endl;
#endif
    timer_pause(&timers[TIMER_TOTAL]);
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
