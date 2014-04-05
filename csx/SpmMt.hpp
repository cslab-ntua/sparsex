/*
 * SpmMt.hpp
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPM_MT_HPP
#define SPM_MT_HPP

#include "Map.hpp"

#include <stdbool.h>
#include <stdint.h>

struct vec;

struct spm_mt_thread {
	void *spm;
	void *spmv_fn;
	unsigned int cpu;
	unsigned int id;
	int node;
	uint64_t row_start;
	uint64_t end;
	uint64_t nr_rows;
	// uint64_t *col_map;
	map_t *map;
	struct vec *x, *y;
	uint64_t size_assigned;
	double secs;
	void *data;
    bool local_sense;
};
typedef struct spm_mt_thread spm_mt_thread_t;

struct spm_mt {
	spm_mt_thread_t *spm_threads;
	unsigned int nr_threads;
	bool symmetric;
#ifdef SPM_NUMA
    bool interleaved;
#endif
    struct vec **local_buffers;
};
typedef struct spm_mt spm_mt_t;

#endif /* SPM_MT_HPP */
