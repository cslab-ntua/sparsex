/*
 * Copyright (C) 2007-2014, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Affinity.cpp
 * \brief Routine that sets the cpu affinity
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/Affinity.hpp>
#include <sparsex/internals/logger/Logger.hpp>
#include <stdlib.h>
#include <sched.h>

namespace sparsex {
  namespace runtime {

    void setaffinity_oncpu(unsigned int cpu)
    {
      cpu_set_t cpu_mask;

      CPU_ZERO(&cpu_mask);
      CPU_SET(cpu, &cpu_mask);

      int err = sched_setaffinity(0, sizeof(cpu_set_t), &cpu_mask);
      if (err) {
	LOG_ERROR << "sched_setaffinity failed\n";
	exit(1);
      }
    }

    int get_cpu_count()
    {
      cpu_set_t cpu_mask;

      CPU_ZERO(&cpu_mask);
      int err = sched_getaffinity(0, sizeof(cpu_set_t), &cpu_mask);
      if (err) {
	LOG_ERROR << "sched_getaffinity failed\n";
	exit(1);
      }

      int count = CPU_COUNT(&cpu_mask);
      printf("%d\n", count);
    
      return count;
    }

  } // end of namespace runtime
} // end of namespace sparsex
