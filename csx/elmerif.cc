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
#include <cassert>
#include <cfloat>
#include "SparseInternal.h"
#include "spmv.h"
#include "jit.h"
#include "TimingFramework.h"

extern "C" {
#include "elmerif.h"
#include "mt_lib.h"
#include "spmv_matvec_mt.h"
#include "spm_crs_mt.h"
//#include "timer.h"
}

//
// Assume zero-indexing by default
// 
#ifndef CSX_IDX_OFFSET
#   define CSX_IDX_OFFSET 0
#endif

enum {
    TIMER_CONSTRUCT_INTERN = 0,
    TIMER_CONSTRUCT_CSX,
    TIMER_SPMV,
    TIMER_TOTAL,
    TIMER_ENDD,
};

static uint64_t nr_calls = 0;
//xtimer_t timers[TIMER_END];
csx::Timer timers[TIMER_ENDD];

const char *timer_desc[] = {
    "Convert to internal repr.",
    "Convert to CSX",
    "Multiplication",
    "Total time in library",
};

void __attribute__ ((constructor)) init() {
    for (int i = 0; i < TIMER_ENDD; ++i) {
        //timer_init(&timers[i]);
        timers[i].SetDescription(timer_desc[i]);
    }
}

void __attribute__ ((destructor)) finalize() {
    std::cout << "libcsx called " << nr_calls
              << " time(s)" << std::endl;

    // Print timer statistics
    for (int i = 0; i < TIMER_ENDD; ++i) {
        //std::cout << timer_desc[i] << ": "
        //          << timer_secs(&timers[i]) << std::endl;
        std::cout << timers[i].GetDescription() << ": "
                  << timers[i].ElapsedTime() << std::endl;
    }
}

#ifdef _DEBUG_
static __thread elmer_value_t *y_copy = 0;
static __thread elmer_value_t *values_last = 0;
static __thread struct csr_s {
    elmer_index_t *rowptr_;
    elmer_index_t *colind_;
    elmer_value_t *values_;
    elmer_index_t nr_rows_;
} shadow_csr_;

static void matvec_csr(struct csr_s *csr,
                       const elmer_value_t *x,
                       elmer_value_t *y)
{
    elmer_index_t i, j;
    for (i = 0; i < csr->nr_rows_; ++i) {
        register elmer_value_t yi_ = 0;
        for (j = csr->rowptr_[i] - CSX_IDX_OFFSET;
             j < csr->rowptr_[i+1] - CSX_IDX_OFFSET; ++j) {
            yi_ += csr->values_[j]*x[csr->colind_[j] - CSX_IDX_OFFSET];
        }

        y[i] = yi_;
    }
}

static bool equal(elmer_value_t *v1, elmer_value_t *v2, elmer_index_t n)
{
    for (elmer_index_t i = 0; i < n; ++i) {
        if (fabs(v1[i] - v2[i]) > 1e-6) {
            printf("Element %d differs: %.6lf != %.6lf\n", i, v1[i], v2[i]);
            return false;
        }
    }

    return true;
}

#endif  // _DEBUG_

void *csx_mattune(elmer_index_t *rowptr,
                  elmer_index_t *colind,
                  elmer_value_t *values,
                  elmer_index_t nr_rows,
                  elmer_index_t nr_cols)
{
    unsigned int nr_threads;
    unsigned int *cpus __attribute__ ((unused));
    csx::SparseInternal<elmer_index_t, elmer_value_t> *spms;
    spm_mt_t *spm_mt;
    uint64_t nr_nzeros = rowptr[nr_rows] - CSX_IDX_OFFSET;
    uint64_t ws = (nr_rows + nr_nzeros + 1)*sizeof(elmer_index_t) +
        nr_nzeros*sizeof(elmer_value_t);
    std::cout << "Rows: " << nr_rows << " "
              << "Cols: " << nr_cols << " "
              << "Non-zeros: " << nr_nzeros << " "
              << "Working set size (MB): " << ((double) ws) / (1024*1024)
              << std::endl;
    mt_get_options(&nr_threads, &cpus);
#ifdef _CSR_
    spm_mt = spm_crs32_double_mt_get_spm((uint32_t *) rowptr,
                                         (uint32_t *) colind,
                                         values, nr_rows, nr_threads, cpus);
#else
    // Initialization of runtime configuration
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    CsxContext csx_context;
    
    Configuration config;
    config = ConfigFromEnv(config, false, true);

    rt_context.SetRuntimeContext(config);
    csx_context.SetCsxContext(config);

    //timer_start(&timers[TIMER_CONSTRUCT_INTERN]);
    timers[TIMER_CONSTRUCT_INTERN].Start();
    spms = csx::LoadCSR_mt<elmer_index_t, elmer_value_t>
        (rowptr, colind, values, nr_rows, nr_cols, !CSX_IDX_OFFSET, nr_threads);
    //timer_pause(&timers[TIMER_CONSTRUCT_INTERN]);
    timers[TIMER_CONSTRUCT_INTERN].Pause();

    //timer_start(&timers[TIMER_CONSTRUCT_CSX]);
    timers[TIMER_CONSTRUCT_CSX].Start();
    //CsxExecutionEngine &engine = CsxJitInit();
    //spm_mt = GetSpmMt(NULL, engine, true, false, spms);
    double pre_time;
    spm_mt = BuildCsx<elmer_index_t, elmer_value_t>(spms, rt_context, csx_context,
                                                    pre_time);
    //timer_pause(&timers[TIMER_CONSTRUCT_CSX]);
    timers[TIMER_CONSTRUCT_CSX].Pause();
    delete[] spms;
#endif        

#ifdef _DEBUG_
    // Initialize the shadow CSR for debugging
    shadow_csr_.rowptr_ = rowptr;
    shadow_csr_.colind_ = colind;
    shadow_csr_.values_ = values;
    shadow_csr_.nr_rows_ = nr_rows;
#endif  // _DEBUG_
    return spm_mt;
}

void csx_matvec(void *spm, elmer_value_t *x, elmer_index_t nr_x,
                elmer_value_t *y, elmer_index_t nr_y)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    // FIXME: Although Elmer uses doubles by default, this is not portable
    vector_double_t *vec_x = vector_double_create_from_buff(x, nr_x);
    vector_double_t *vec_y = vector_double_create_from_buff(y, nr_y);

#ifdef _DEBUG_
    if (!y_copy)
        y_copy = (elmer_value_t *) malloc(nr_y*sizeof(*y_copy));
    memcpy(y_copy, y, nr_y*sizeof(*y_copy));
    matvec_csr(&shadow_csr_, x, y_copy);
#endif  // _DEBUG_

    ++nr_calls;
    //timer_start(&timers[TIMER_SPMV]);
    timers[TIMER_SPMV].Start();
    spmv_double_matvec_mt(spm_mt, vec_x, vec_y);
    //timer_pause(&timers[TIMER_SPMV]);
    timers[TIMER_SPMV].Pause();

#ifdef _DEBUG_
    std::cout << "Checking result... ";
    if (equal(vec_y->elements, y_copy, nr_y))
        std::cout << "DONE" << std::endl;
    else
        std::cout << "FAILED" << std::endl;
#endif  // _DEBUG_
}

void elmer_matvec_(void **tuned, void *n, void *rowptr, void *colind,
                   void *values, void *u, void *v, void *reinit)
{
    elmer_index_t n_ = *((elmer_index_t *) n);
    elmer_index_t *rowptr_ = (elmer_index_t *) rowptr;
    elmer_index_t *colind_ = (elmer_index_t *) colind;
    elmer_value_t *values_ = (elmer_value_t *) values;
    bool reinit_ = *((bool *) reinit);
    elmer_value_t *x_ = (elmer_value_t *) u;
    elmer_value_t *y_ = (elmer_value_t *) v;

    //timer_start(&timers[TIMER_TOTAL]);
    timers[TIMER_TOTAL].Start();
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

    if (!*tuned || reinit_) {
        *tuned = csx_mattune(rowptr_, colind_, values_, n_, n_);
    }

    csx_matvec(*tuned, x_, n_, y_, n_);
    //timer_pause(&timers[TIMER_TOTAL]);
    timers[TIMER_SPMV].Pause();
}



// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
