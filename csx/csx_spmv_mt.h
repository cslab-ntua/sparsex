/*
 * csx_spmv_mt.h
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_SPMV_MT_H__
#define CSX_SPMV_MT_H__

#include "spmv_method.h"
#include "vector.h"
#include "spm_mt.h"
#include "affinity.h"
#include "tsc.h"
#include "csr.h"

#ifdef SPMV_PRFCNT
#include "prfcnt.h"
#endif /* SPMV_PRFCNT */

#include <pthread.h>

float spmv_bench_mt(spm_mt_t *spm_mt, unsigned long loops,
                    unsigned long rows_nr, unsigned long cols_nr);
void spmv_check_mt(csx::CSR<uint64_t, double> *spm, spm_mt_t *spm_mt,
                   unsigned long loops, unsigned long rows_nr,
                   unsigned long cols_nr);

#endif  // CSX_SPMV_MT_H__

