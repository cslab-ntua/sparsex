/*
 * spm_vbl_mt.h -- multithreaded VBL
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __SPM_VBL_MT_H__
#define __SPM_VBL_MT_H__

#include <inttypes.h>

#include "spm_vbl.h"
#include "../spm_mt.h"
#include "../spmv_method.h"

#define SPM_VBL_MT_DECLARE(__idx_bits, __elem_type) \
struct spm_vbl ## __idx_bits ## _ ## __elem_type ## _mt { \
	spm_vbl ## __idx_bits ## _ ## __elem_type ## _t    *vbl; \
	uint64_t row_start, row_end; \
    uint64_t bstart; \
    uint64_t bend;   \
}; \
typedef struct spm_vbl ## __idx_bits ## _ ## __elem_type ## _mt spm_vbl ## __idx_bits ## _ ## __elem_type ## _mt ## _t; \
\
void * \
spm_vbl ## __idx_bits ## _ ## __elem_type ## _mt_init_mmf( \
	char *mmf_file, \
	uint64_t *rows_nr, uint64_t *cols_nr, \
	uint64_t *nz_nr, void *metadata);     \
\
/* XXX: Destroy */ \
\
spmv_ ## __elem_type ## _fn_t spm_vbl ## __idx_bits ## _ ## __elem_type ## _mt_multiply;

SPM_VBL_MT_DECLARE(32, double)
SPM_VBL_MT_DECLARE(64, double)
SPM_VBL_MT_DECLARE(32, float)
SPM_VBL_MT_DECLARE(64, float)

#include "../macros.h"
#define SPM_VBL_MT_NAME(name) CON6(spm_vbl, SPM_CRS_BITS, _,ELEM_TYPE,_mt,name)
#define SPM_VBL_MT_TYPE SPM_VBL_MT_NAME(_t)

#endif
