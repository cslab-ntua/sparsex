/*
 * Affinity.cpp -- CPU related utilities.
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "Affinity.hpp"

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

// int get_cpu_count2()
// {
// #ifdef _SC
//     return ;
// #else
//     int count = sysconf(_SC_NR_PROCESSORS_ONLN);
//     return count;
// #endif
// }
