/*
 * spm_bcsr.h -- BCSR format
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __SPM_BCSR_H__
#define __SPM_BCSR_H__

#include <inttypes.h>
#include "spmv_method.h"
#include "spm_crs.h"
#include "export.h"

#define BCSR_PAD_MAX(r, c)  (r*c - 1)

#define BCSR_NAME(val_type, ci_bits, name) spm_bcsr ## ci_bits ## _ ## val_type ## name
#define BCSR_TYPE(val_type, ci_bits) BCSR_NAME(val_type, ci_bits, _t)

#define SPM_BCSR_DECLARE(val_type, ci_bits) \
    typedef struct {                                \
        val_type            *bvalues;               \
        UINT_TYPE(ci_bits)  *bcol_ind, *brow_ptr;   \
        uint64_t            br, bc;                 \
        uint64_t            nz, nrows, ncols;       \
        uint64_t            nblocks, nbrows;        \
        int                 storage;                \
    } BCSR_TYPE(val_type, ci_bits);                 \
\
void *\
BCSR_NAME(val_type, ci_bits, _init_mmf)(char *mmf_file,                 \
                                        uint64_t *rows_nr, uint64_t *cols_nr, \
                                        uint64_t *nz_nr, void *metadata); \
void BCSR_NAME(val_type, ci_bits, _destroy)(void *crs); \
uint64_t BCSR_NAME(val_type, ci_bits, _size)(void *spm); \
spmv_  ## val_type ## _fn_t BCSR_NAME(val_type, ci_bits, _multiply);

SPM_BCSR_DECLARE(double, 32)
SPM_BCSR_DECLARE(double, 64)
SPM_BCSR_DECLARE(float, 32)
SPM_BCSR_DECLARE(float, 64)

#include "macros.h"
#define SPM_BCSR_NAME(name) CON5(spm_bcsr, SPM_CRS_BITS, _, ELEM_TYPE, name)
#define SPM_BCSR_TYPE SPM_BCSR_NAME(_t)

SPM_BCSR_TYPE *SPM_BCSR_NAME(_init_crs)(const SPM_CRS_TYPE *crs, uint64_t r, uint64_t c, uint64_t pad);

#endif
