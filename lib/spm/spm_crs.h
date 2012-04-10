/*
 * spm_crs.h
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __SPM_CRS_H__
#define __SPM_CRS_H__

#include <inttypes.h>
#include "spmv_method.h"

#ifndef SPM_CRS_BITS
#define SPM_CRS_BITS 64
#endif

#define CRS_NAME(val_type, ci_bits, name)       \
    spm_crs ## ci_bits ## _ ## val_type ## name
#define CRS_TYPE(val_type, ci_bits) CRS_NAME(val_type, ci_bits, _t)

#define SPM_CRS_DECLARE(val_type, ci_bits) \
typedef struct { \
	val_type            *values; \
	UINT_TYPE(ci_bits)  *col_ind, *row_ptr; \
	uint64_t            nz, nrows, ncols; \
} CRS_TYPE(val_type, ci_bits); \
\
void *\
CRS_NAME(val_type, ci_bits, _init_mmf)(char *mmf_file,                  \
                                       uint64_t *rows_nr, uint64_t *cols_nr, \
                                       uint64_t *nz_nr, void *metadata); \
void CRS_NAME(val_type, ci_bits, _destroy)(void *crs); \
uint64_t CRS_NAME(val_type, ci_bits, _size)(void *spm); \
spmv_  ## val_type ## _fn_t CRS_NAME(val_type, ci_bits, _multiply);

SPM_CRS_DECLARE(double, 32)
SPM_CRS_DECLARE(double, 64)
SPM_CRS_DECLARE(float, 32)
SPM_CRS_DECLARE(float, 64)

#include "macros.h"
#define SPM_CRS_NAME(name) CON5(spm_crs, SPM_CRS_BITS, _, ELEM_TYPE, name)
#define SPM_CRS_TYPE SPM_CRS_NAME(_t)
#define SPM_CRS_IDX_TYPE UINT_TYPE(SPM_CRS_BITS)

#endif
