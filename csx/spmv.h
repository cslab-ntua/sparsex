/* -*- C++ -*-
 *
 * spmv.h --  Utility routines for setting up the SpMV execution.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * Copyright (C) 2011, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPMV_H__
#define SPMV_H__

#include <stdint.h>
#include <pthread.h>
#include <sstream>

#include "spm.h"
#include "csx.h"

extern "C" {
#include "spm_mt.h"
}

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

int SplitString(char *str, char **str_buf, const char *start_sep,
                const char *end_sep);

/*
 *  Utility routines to retrieve the values of the controlling environment
 *  variables.
 */
void GetOptionXform(int **xform_buf);
void GetOptionEncodeDeltas(int ***deltas);
uint64_t GetOptionWindowSize();
uint64_t GetOptionSamples();
double GetOptionProbability();
bool GetOptionSplitBlocks();


/**
 *  Parallel Preprocessing.
 *
 *  @param info parameters of matrix needed by each thread.
 */
void *PreprocessThread(void *thread_info);

/**
 *  Routine responsible for making a map for the symmetric representation of
 *  sparse matrices.
 */
void MakeMap(spm_mt_t *spm_mt, csx::SPMSym *spm_sym);
uint64_t MapSize(void *spm);

#ifdef __cplusplus                  /* included in a C++ program */
spm_mt_t *GetSpmMt(char *mmf_fname, csx::SPM *Spms = 0);
#else                               /* included in a C program */
spm_mt_t *GetSpmMt(char *mmf_fname, csx::SPM *Spms = 0);
#endif
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
unsigned long CsxSize(spm_mt_t *spm_mt);
uint64_t CsxSymSize(void *spm_mt);
unsigned long CsxSymSize(spm_mt_t *spm_mt);

#endif // SPMV_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
