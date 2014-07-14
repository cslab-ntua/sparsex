/*
 * \file SpmMt.hpp
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

#ifndef SPARSEX_INTERNALS_SPM_MT_HPP
#define SPARSEX_INTERNALS_SPM_MT_HPP

#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/Map.hpp>
#include <sparsex/internals/Vector.hpp>
#include <sparsex/internals/SpmvMethod.hpp>
#include <stdbool.h>
#include <stdint.h>

struct spm_mt_thread {
	void *spm;
	spmv_fn_t spmv_fn;
	unsigned int cpu;
	unsigned int id;
	int node;
	size_t row_start;
	size_t end;
	size_t nr_rows;
	// uint64_t *col_map;
	map_t *map;
	vector_t *x, *y;
	size_t size_assigned;
	double secs;
	void *data;
    int *sense;
};
typedef struct spm_mt_thread spm_mt_thread_t;

struct spm_mt {
	spm_mt_thread_t *spm_threads;
	unsigned int nr_threads;
	bool symmetric;
#if SPX_USE_NUMA
    bool interleaved;
#endif
    vector_t **local_buffers;
};

typedef struct spm_mt spm_mt_t;

#endif // SPARSEX_INTERNALS_SPM_MT_HPP
