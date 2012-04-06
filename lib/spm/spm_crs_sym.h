/*
 * spm_crs_sym.h -- Interface of expansion of CSR for symmetric sparse matrices.
 *
 * Copyright (C) 2011,      Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPM_CRS_SYM_H_
#define SPM_CRS_SYM_H_

#include <inttypes.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

#include "macros.h"
#include "vector.h"
#include "mmf.h"
#include "spmv_method.h"
#include "spm_crs.h"

#define CRS_SYM_NAME(val_type, ci_bits, name) spm_crs ## ci_bits ## _ ## val_type ## _sym ## name
#define CRS_SYM_TYPE(val_type, ci_bits) CRS_SYM_NAME(val_type, ci_bits, _t)

#define SPM_CRS_SYM_DECLARE(val_type, ci_bits) \
typedef struct { \
    val_type            *values; \
    UINT_TYPE(ci_bits)  *col_ind, *row_ptr; \
    val_type            *dvalues; \
    uint64_t            nnz, n; \
} CRS_SYM_TYPE(val_type, ci_bits); \
\
void *\
CRS_SYM_NAME(val_type, ci_bits, _init_mmf)(char *mmf_file, \
                                           uint64_t *nrows, uint64_t *ncols, \
                                           uint64_t *nnz, void *metadata); \
void CRS_SYM_NAME(val_type, ci_bits, _destroy)(void *crs); \
uint64_t CRS_SYM_NAME(val_type, ci_bits, _size)(void *spm); \
spmv_  ## val_type ## _fn_t CRS_SYM_NAME(val_type, ci_bits, _multiply);

SPM_CRS_SYM_DECLARE(double, 32)
SPM_CRS_SYM_DECLARE(double, 64)
SPM_CRS_SYM_DECLARE(float, 32)
SPM_CRS_SYM_DECLARE(float, 64)

#define SPM_CRS_SYM_NAME(name) CON6(spm_crs, SPM_CRS_BITS, _, ELEM_TYPE, _sym,  name)
#define SPM_CRS_SYM_TYPE SPM_CRS_SYM_NAME(_t)
#define SPM_CRS_SYM_IDX_TYPE UINT_TYPE(SPM_CRS_BITS)

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
