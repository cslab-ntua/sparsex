/*
 * Copyright (C) 2012-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file CsxSpmv.cpp
 * \brief Thread-level SpMV routines
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/CsxSpmv.hpp>

using namespace sparsex::runtime;

extern vector_t **temp;
extern size_t nr_threads;

namespace sparsex {

void do_mv_thread(void *args)
{
	spm_mt_thread_t *spm_thread = (spm_mt_thread_t *) args;
	setaffinity_oncpu(spm_thread->cpu);
	spmv_fn_t fn = spm_thread->spmv_fn;
    double ALPHA = spm_thread->alpha;

	fn(spm_thread->csx, spm_thread->x, spm_thread->y, ALPHA, NULL);
}

void do_mv_sym_thread(void *args)
{
	spm_mt_thread_t *spm_thread = (spm_mt_thread_t *) args;
	setaffinity_oncpu(spm_thread->cpu);
	spmv_fn_t fn = spm_thread->spmv_fn;
    double ALPHA = spm_thread->alpha;
	int id = spm_thread->id;

    VecInitFromMap(temp, 0, spm_thread->map);
    centralized_barrier(&(spm_thread->sense), nr_threads);
    fn(spm_thread->csx, spm_thread->x, spm_thread->y, ALPHA, temp[id]);
    centralized_barrier(&(spm_thread->sense), nr_threads);
    /* Switch Reduction Phase */
    VecAddFromMap(spm_thread->y, temp, spm_thread->y, spm_thread->map);
}

void do_kernel_thread(void *args)
{
	spm_mt_thread_t *spm_thread = (spm_mt_thread_t *) args;
	setaffinity_oncpu(spm_thread->cpu);
	spmv_fn_t fn = spm_thread->spmv_fn;
    double ALPHA = spm_thread->alpha;
    double BETA = spm_thread->beta;
    int start = spm_thread->row_start;
    int end = start + spm_thread->nr_rows;

	if (BETA != 1)
        VecScalePart(spm_thread->y, spm_thread->y, BETA, start, end);

	fn(spm_thread->csx, spm_thread->x, spm_thread->y, ALPHA, NULL);
}

void do_kernel_sym_thread(void *args)
{
    spm_mt_thread_t *spm_thread = (spm_mt_thread_t *) args;
	setaffinity_oncpu(spm_thread->cpu);
	spmv_fn_t fn = spm_thread->spmv_fn;
    double ALPHA = spm_thread->alpha;
    double BETA = spm_thread->beta;
    int start = spm_thread->row_start;
    int end = start + spm_thread->nr_rows;
    int id = spm_thread->id;

    VecInitFromMap(temp, 0, spm_thread->map);
    centralized_barrier(&(spm_thread->sense), nr_threads);
	if (BETA != 1)
        VecScalePart(spm_thread->y, spm_thread->y, BETA, start, end);

    fn(spm_thread->csx, spm_thread->x, spm_thread->y, ALPHA, temp[id]);
    centralized_barrier(&(spm_thread->sense), nr_threads);
    /* Switch Reduction Phase */
    VecAddFromMap(spm_thread->y, temp, spm_thread->y, spm_thread->map);
}

} // end of namespace sparsex
