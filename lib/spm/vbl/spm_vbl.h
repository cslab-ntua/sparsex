/*
 * spm_vbl.h -- VBL format
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __SPM_VBL_H__
#define __SPM_VBL_H__

#include <inttypes.h>
#include "spmv_method.h"

#define VBL_NAME(val_type, ci_bits, name) spm_vbl ## ci_bits ## _ ## val_type ## name
#define VBL_TYPE(val_type, ci_bits) VBL_NAME(val_type, ci_bits, _t)

#define SPM_VBL_DECLARE(val_type, ci_bits) \
    typedef struct {                                \
        val_type            *values;                \
        UINT_TYPE(ci_bits)  *bcol_ind, *row_ptr;    \
        uint8_t             *bsize;                 \
        uint64_t            nz, nrows, ncols;       \
        uint64_t            nblocks;                \
    } VBL_TYPE(val_type, ci_bits);                  \
\
void *\
VBL_NAME(val_type, ci_bits, _init_mmf)(char *mmf_file, \
                                       uint64_t *rows_nr, uint64_t *cols_nr, \
                                       uint64_t *nz_nr, void *metadata);               \
void VBL_NAME(val_type, ci_bits, _destroy)(void *crs); \
uint64_t VBL_NAME(val_type, ci_bits, _size)(void *spm); \
spmv_  ## val_type ## _fn_t VBL_NAME(val_type, ci_bits, _multiply);

SPM_VBL_DECLARE(double, 32)
SPM_VBL_DECLARE(double, 64)
SPM_VBL_DECLARE(float, 32)
SPM_VBL_DECLARE(float, 64)

#include "macros.h"
#define SPM_VBL_NAME(name) CON5(spm_vbl, SPM_CRS_BITS, _, ELEM_TYPE, name)
#define SPM_VBL_TYPE SPM_VBL_NAME(_t)
#define SPM_VBL_IDX_TYPE UINT_TYPE(SPM_CRS_BITS)

SPM_VBL_TYPE *SPM_VBL_NAME(_init_crs)(SPM_CRS_TYPE *crs);

#endif
