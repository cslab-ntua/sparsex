/*
 * spmv.cc -- Front-end utilities for invoking CSX.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "spmv.h"
#include "CsxUtil.h"

// malloc wrapper
#define xmalloc(x)                         \
({                                         \
    void *ret_;                            \
    ret_ = malloc(x);                      \
    if (ret_ == NULL){                     \
        std::cerr << __FUNCTION__          \
                  << " " << __FILE__       \
                  << ":" << __LINE__       \
                  << ": malloc failed\n";  \
        exit(1);                           \
    }                                      \
    ret_;                                  \
})

static int GetOptionOuterLoops()
{
    const char *loops_env = getenv("OUTER_LOOPS");
    int ret = 1;
    
    if (loops_env) {
        ret = atoi(loops_env);
        if (ret < 0)
            ret = 0;
    }
    
    return ret;
}

uint64_t MapSize(void *spm)
{
    unsigned int i;
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    spm_mt_thread_t *spm_thread;
    uint64_t size = 0;
    
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        size += spm_thread->map->length * sizeof(uint32_t);
        size += spm_thread->map->length * sizeof(uint32_t);
        size += spm_thread->map->length * sizeof(double);
    }
    
    return size;
}

spm_mt_t *PrepareSpmMt(const RuntimeContext &rt_config,
                       const CsxContext &csx_config)
{
    spm_mt_t *spm_mt;
    spm_mt = (spm_mt_t *) xmalloc(sizeof(spm_mt_t));

    spm_mt->nr_threads = rt_config.GetNrThreads();
    spm_mt->symmetric = csx_config.IsSymmetric();
    spm_mt->spm_threads = (spm_mt_thread_t *) xmalloc
        (sizeof(spm_mt_thread_t) * rt_config.GetNrThreads());

    for (size_t i = 0; i < rt_config.GetNrThreads(); i++) {
        spm_mt->spm_threads[i].cpu = rt_config.GetAffinity(i);
        spm_mt->spm_threads[i].node =
            numa_node_of_cpu(rt_config.GetAffinity(i));
        spm_mt->spm_threads[i].id = i;
    }

    return spm_mt;
}

void CheckLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    void *crs;
    uint64_t nrows, ncols, nnz;

    crs = spm_crs32_double_init_mmf(mmf_name, &nrows, &ncols, &nnz, NULL);
    std::cout << "Checking ... " << std::flush;
    if (!spm_mt->symmetric)
        SPMV_CHECK_FN(crs, spm_mt, spm_crs32_double_multiply, 1, nrows, ncols,
                      NULL);
    else
        SPMV_CHECK_SYM_FN(crs, spm_mt, spm_crs32_double_multiply, 1, nrows,
                          ncols, NULL);

    spm_crs32_double_destroy(crs);
    std::cout << "Check Passed" << std::endl;
}

void BenchLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    uint64_t nrows, ncols, nnz;
    double secs, flops;
    long loops_nr = 128;

    ReadMmfSizeLine(mmf_name, nrows, ncols, nnz);
    int nr_outer_loops = GetOptionOuterLoops();
    
    for (int i = 0; i < nr_outer_loops; ++i) {
        if (!spm_mt->symmetric) {
            secs = SPMV_BENCH_FN(spm_mt, loops_nr, nrows, ncols, NULL);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            printf("m:%s f:%s s:%lu pt:%lf t:%lf r:%lf\n", "csx",
                   basename(mmf_name), CsxSize<int, double>(spm_mt),
                   pre_time, secs, flops);
        } else {
            secs = SPMV_BENCH_SYM_FN(spm_mt, loops_nr, nrows, ncols, NULL);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            // Switch Reduction Phase
            printf("m:%s f:%s ms:%lu s:%lu pt:%lf t:%lf r:%lf\n", "csx-sym",
                   basename(mmf_name), MapSize(spm_mt),
                   CsxSymSize<int, double>(spm_mt), pre_time, secs, flops);
        }
    }
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
