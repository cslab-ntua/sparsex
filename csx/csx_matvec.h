/* -*- C++ -*-
 *
 * csx_matvec.h -- Multithreaded kernel y <-- alpha*A*x + beta*y
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_MATVEC_H__
#define CSX_MATVEC_H__

#include "spmv_method.h"
#include "spm_mt.h"
#include "vector.h"

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct matvec_params
{
    spm_mt_thread_t     *spm_thread;
    vector_t            *tmp;
    int                 start;
    int                 end;
} matvec_params;

void matvec_mt(spm_mt_t *spm_mt, vector_t *x, double alpha, vector_t *y,
               double beta);
void matvec_sym_mt(spm_mt_t *spm_mt, vector_t *x, double alpha, vector_t *y,
                   double beta);

#ifdef __cplusplus
}
#endif

#endif // CSX_MATVEC_H__
