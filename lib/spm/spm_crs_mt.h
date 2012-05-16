/*
 * spm_crs_mt.h -- Multithreaded CSR
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
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
	spm_crs ## __idx_bits ## _ ## __elem_type ## _t *crs; \
	uint64_t row_start, row_end; \
}; \
typedef struct spm_crs ## __idx_bits ## _ ## __elem_type ## _mt \
               spm_crs ## __idx_bits ## _ ## __elem_type ## _mt ## _t; \
\
void * spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_init_mmf( \
    char *mmf_file, uint64_t *rows_nr, uint64_t *cols_nr, uint64_t *nz_nr, \
    void *metadata); \
\
uint64_t spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_size(void *spm); \
\
void spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_destroy(void *spm); \
\
spmv_ ## __elem_type ## _fn_t \
    spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_multiply;

SPM_CRS_MT_DECLARE(32, double)
SPM_CRS_MT_DECLARE(64, double)
SPM_CRS_MT_DECLARE(32, float)
SPM_CRS_MT_DECLARE(64, float)

#ifdef SPM_NUMA
#define SPM_CRS_MT_NUMA_DECLARE(__idx_bits, __elem_type) \
struct spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_numa { \
	spm_crs ## __idx_bits ## _ ## __elem_type ## _t *crs; \
	uint64_t row_start, row_end; \
}; \
typedef struct spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_numa \
               spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_numa ## _t; \
\
void * spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_numa_init_mmf( \
           char *mmf_file, uint64_t *rows_nr, uint64_t *cols_nr, \
           uint64_t *nz_nr, void *metadata); \
\
uint64_t \
spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_numa_size(void *spm); \
\
void spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_numa_destroy(void *spm); \
\
spmv_ ## __elem_type ## _fn_t \
    spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_numa_multiply;
// void
// spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_numa_make_col_map(void *spm);

SPM_CRS_MT_NUMA_DECLARE(32, double)
SPM_CRS_MT_NUMA_DECLARE(64, double)
SPM_CRS_MT_NUMA_DECLARE(32, float)
SPM_CRS_MT_NUMA_DECLARE(64, float)
#endif /* SPM_NUMA */

#include "macros.h"
#define SPM_CRS_MT_NAME(name) CON6(spm_crs, SPM_CRS_BITS, _,ELEM_TYPE,_mt,name)
#define SPM_CRS_MT_TYPE SPM_CRS_MT_NAME(_t)

spm_mt_t *SPM_CRS_MT_NAME(_get_spm)(SPM_CRS_IDX_TYPE *rowptr,
                                    SPM_CRS_IDX_TYPE *colind,
                                    ELEM_TYPE *values,
                                    SPM_CRS_IDX_TYPE nr_rows,
                                    unsigned int nr_threads,
                                    unsigned int *cpus);

void SPM_CRS_MT_NAME(_multiply_base_one)(void *spm, VECTOR_TYPE *in,
                                         VECTOR_TYPE *out);

#endif /* __SPM_CRS_MT_H__ */
