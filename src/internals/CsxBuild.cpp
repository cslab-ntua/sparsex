/*
 * \file CsxBuild.cpp
 *
 * \brief Front-end utilities for building CSX
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/internals/CsxBuild.hpp>

namespace sparsex {
namespace csx {

spm_mt_t *PrepareSpmMt()
{
    spm_mt_t *spm_mt;
    RuntimeConfiguration &rt_config = RuntimeConfiguration::GetInstance();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();

    spm_mt = new spm_mt_t;
    spm_mt->nr_threads = rt_context.GetNrThreads();
    spm_mt->local_buffers = NULL;
#if SPX_USE_NUMA
    spm_mt->interleaved = false;
#endif
    spm_mt->symmetric =
        rt_config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric);
    spm_mt->spm_threads = new spm_mt_thread_t[rt_context.GetNrThreads()];
    for (size_t i = 0; i < rt_context.GetNrThreads(); i++) {
        spm_mt->spm_threads[i].cpu = rt_context.GetAffinity(i);
        spm_mt->spm_threads[i].node =
            numa_node_of_cpu(rt_context.GetAffinity(i));
        spm_mt->spm_threads[i].id = i;
        spm_mt->spm_threads[i].sense = NULL;
    }

    return spm_mt;
}

} // end of namespace csx
} // end of namespace sparsex

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
