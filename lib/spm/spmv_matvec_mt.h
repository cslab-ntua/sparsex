/*
 * spmv_matvec_mt.h -- Multithread matvec product routines.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPMV_MATVEC_MT_H__
#define SPMV_MATVEC_MT_H__

#include "spm_mt.h"
#include "vector.h"

void spmv_double_matvec_mt(spm_mt_t *spm_mt, VECTOR_TYPE *x, VECTOR_TYPE *y);
void spmv_float_matvec_mt(spm_mt_t *spm_mt, VECTOR_TYPE *x, VECTOR_TYPE *y);

#endif  /* SPMV_MATVEC_MT_H__ */
