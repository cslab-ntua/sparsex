/* -*- C++ -*-
 * 
 * spmv.h -- Front-end utilities for invoking CSX.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPMV_H__
#define SPMV_H__

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <numa.h>
#include <sched.h>

#include "spm.h"
#include "mmf.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"
#include "llvm_jit_help.h"
#include "spmv.h"

#include <numa.h>
#include <numaif.h>

extern "C" {
#include <libgen.h>
#include "mt_lib.h"
#include "spmv_method.h"
#include "spm_crs.h"
#include "spm_mt.h"
#ifdef SPM_NUMA
    #include "spmv_loops_mt_numa.h"
    #include "spmv_loops_sym_mt_numa.h"
    #define SPMV_CHECK_FN spmv_double_check_mt_loop_numa
    #define SPMV_BENCH_FN spmv_double_bench_mt_loop_numa
    #define SPMV_CHECK_SYM_FN spmv_double_check_sym_mt_loop_numa
    #define SPMV_BENCH_SYM_FN spmv_double_bench_sym_mt_loop_numa
#else
    #include "spmv_loops_mt.h"
    #include "spmv_loops_sym_mt.h"
    #define SPMV_CHECK_FN spmv_double_check_mt_loop
    #define SPMV_BENCH_FN spmv_double_bench_mt_loop
    #define SPMV_CHECK_SYM_FN spmv_double_check_sym_mt_loop
    #define SPMV_BENCH_SYM_FN spmv_double_bench_sym_mt_loop
#endif // SPM_NUMA
#include "timer.h"
#include "ctl_ll.h"
#include <libgen.h>
}

using namespace csx;

double pre_time;

///> Thread data essential for parallel preprocessing
typedef struct thread_info {
    unsigned int thread_no;
    unsigned int cpu;
    csx::SPM *spm;
    csx::SPMSym *spm_sym;
    spm_mt_thread_t *spm_encoded;
    csx::CsxManager *csxmg;
    uint64_t wsize;
    std::ostringstream buffer;
    int *xform_buf;
    double sampling_prob;
    uint64_t samples_max;
    double sampling_portion;
    bool split_blocks;
    bool symmetric;
    int **deltas;
} thread_info_t;

/**
 *  Parallel Preprocessing.
 *
 *  @param thread_info parameters of matrix needed by each thread.
 */
void *PreprocessThread(void *thread_info);

/**
 *  Routine responsible for making a map for the symmetric representation of
 *  sparse matrices.
 *  
 *  @param spm_mt  parameters of the multithreaded execution.
 *  @param spm_sym the multithreaded CSX-Sym matrix.
 */
void MakeMap(spm_mt_t *spm_mt, SPMSym *spm_sym);

/**
 *  Find the size of the map including the values accessed by it.
 *
 *  @param spm_mt  the complete multithreaded CSX-Sym matrix.
 *  @return        the size of the map.
 */
uint64_t MapSize(void *spm);

/**
 *  Routine responsible for retrieving the CSX or CSX-Sym sparse matrix format
 *  according to the command line parameters.
 *  
 *  @param mmf_fname    name of the sparse matrix file.
 *  @param engine       """bkk"""
 *  @param split_blocks whether or not split-blocks update is on/off.
 *  @param symmetric    true if CSX-Sym format is used, false otherwise. 
 *  @param spms         an initial sparse matrix format if it exists.
 *  @return             the (multithreaded) CSX or CSX-Sym sparse matrix.
 */
spm_mt_t *GetSpmMt(char *mmf_fname, CsxExecutionEngine &engine,
                   bool split_blocks, bool symmetric, SPM *spms = NULL);
                   
/**
 *  Deallocation of CSX or CSX-Sym sparse matrix.
 *  
 *  @param spm_mt the (multithreaded) CSX or CSX-Sym sparse matrix.
 */
void PutSpmMt(spm_mt_t *spm_mt);

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

/**
 *  Compute the size (in bytes) of the compressed matrix in CSX form.
 *
 *  @param spm_mt  the sparse matrix in CSX format.
 */
uint64_t CsxSize(void *spm_mt);

/**
 *  Compute the size (in bytes) of the compressed matrix in CSX-Sym form.
 *
 *  @param spm_mt  the sparse matrix in CSX-Sym format.
 */
uint64_t CsxSymSize(void *spm_mt);

#endif // SPMV_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
