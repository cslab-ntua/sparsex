/* -*- C++ -*-
 *
 * affinity.cc -- Routine that sets the cpu affinity.
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "affinity.h"

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
