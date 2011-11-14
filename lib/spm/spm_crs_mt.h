/*
 * spm_crs_mt.h -- multithreaded CRS
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __SPM_CRS_MT_H__
#define __SPM_CRS_MT_H__

#include <inttypes.h>

#include "spm_crs.h"
#include "spm_mt.h"
#include "spmv_method.h"
   
#define SPM_CRS_MT_DECLARE(__idx_bits, __elem_type) \
struct spm_crs ## __idx_bits ## _ ## __elem_type ## _mt { \
	spm_crs ## __idx_bits ## _ ## __elem_type ## _t    *crs; \
	uint64_t row_start, row_end; \
}; \
typedef struct spm_crs ## __idx_bits ## _ ## __elem_type ## _mt spm_crs ## __idx_bits ## _ ## __elem_type ## _mt ## _t; \
\
void * \
spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_init_mmf( \
	char *mmf_file, \
	uint64_t *rows_nr, uint64_t *cols_nr, \
	uint64_t *nz_nr, void *metadata);     \
\
uint64_t \
spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_size(void *spm); \
\
/* XXX: Destroy */ \
\
spmv_ ## __elem_type ## _fn_t spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_multiply;

SPM_CRS_MT_DECLARE(32, double)
SPM_CRS_MT_DECLARE(64, double)
SPM_CRS_MT_DECLARE(32, float)
SPM_CRS_MT_DECLARE(64, float)

#include "macros.h"
#define SPM_CRS_MT_NAME(name) CON6(spm_crs, SPM_CRS_BITS, _,ELEM_TYPE,_mt,name)
#define SPM_CRS_MT_TYPE SPM_CRS_MT_NAME(_t)

#endif
