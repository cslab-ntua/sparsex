/*
 * spmv.h -- Front-end utilities for invoking CSX.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPMV_H
#define SPMV_H

#include "jit.h"

extern "C" {
#include "spm_mt.h"
}

using namespace csx;

spm_mt_t *GetSpmMt(char *mmf_fname, CsxExecutionEngine &engine);

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

#endif // SPMV_H

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
