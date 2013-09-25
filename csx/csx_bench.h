/* -*- C++ -*-
 * 
 * csx_bench.h -- Benchmark utilities.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_BENCH_H__
#define CSX_BENCH_H__

#include "mmf.h"
#include "spm_mt.h"
#include "csx_util.h"
#include "csx_spmv_mt.h"
#include "csxsym_spmv_mt.h"

#include <libgen.h>

#define SPMV_CHECK_FN spmv_check_mt
#define SPMV_BENCH_FN spmv_bench_mt
#define SPMV_CHECK_SYM_FN spmv_check_sym_mt
#define SPMV_BENCH_SYM_FN spmv_bench_sym_mt

/**
 *  Check the CSX SpMV result against the baseline single-thread CSR
 *  implementation.
 *
 *  @param the (mulithreaded) CSX matrix.
 *  @param the MMF input file.
 */
void CheckLoop(spm_mt_t *spm_mt, char *mmf_name);

/**
 *  Run CSX SpMV and record the performance information.
 */
void BenchLoop(spm_mt_t *spm_mt, char *mmf_name);

#endif // CSX_BENCH_H__
