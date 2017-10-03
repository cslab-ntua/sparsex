/*
 * Copyright (C) 2011-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file CsxBuild.cpp
 * \brief Front-end utilities for building CSX
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/CsxBuild.hpp>

namespace sparsex {
  namespace csx {

    spm_mt_t *PrepareSpmMt()
    {
      spm_mt_t *spm_mt;
      RtConfig &rt_config = RtConfig::GetInstance();
      RtCtx &rt_context = RtCtx::GetInstance();

      spm_mt = new spm_mt_t;
      spm_mt->nr_threads = rt_context.GetNrThreads();
      spm_mt->local_buffers = NULL;
#if SPX_USE_NUMA
      spm_mt->interleaved = false;
#endif
      spm_mt->symmetric =
        rt_config.GetProperty<bool>(RtConfig::MatrixSymmetric);
      spm_mt->spm_threads = new spm_mt_thread_t[rt_context.GetNrThreads()];
      for (size_t i = 0; i < rt_context.GetNrThreads(); i++) {
        spm_mt->spm_threads[i].cpu = rt_context.GetAffinity(i);
        spm_mt->spm_threads[i].node =
	  numa_node_of_cpu(rt_context.GetAffinity(i));
        spm_mt->spm_threads[i].id = i;
        spm_mt->spm_threads[i].sense = 0;
      }

      return spm_mt;
    }

  } // end of namespace csx
} // end of namespace sparsex
