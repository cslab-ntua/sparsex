/*
 * csx_build.cc -- Front-end utilities for building CSX.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "csx_build.h"

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

spm_mt_t *PrepareSpmMt(const CsxContext &csx_config)
{
    spm_mt_t *spm_mt;
    RuntimeContext &rt_config = RuntimeContext::GetInstance();

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

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
