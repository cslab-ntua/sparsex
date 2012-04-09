/*
 * spm_crs_sym_mt.h -- Interface of expansion of CSR for symmetric sparse
 *                     matrices.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef __SPM_CRS_SYM_MT_H__
#define __SPM_CRS_SYM_MT_H__

#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <pthread.h>

#include "macros.h"
#include "mt_lib.h"
#include "spm_mt.h"
#include "spm_crs_sym.h"
#include "spmv_method.h"

uint64_t map_size;

#define SPM_CRS_SYM_MT_DECLARE(__idx_bits, __elem_type) \
\
struct spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt { \
	spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_t *crs; \
	uint64_t row_start, row_end, nnz; \
}; \
typedef struct spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt \
                   spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_t; \
\
void *spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_init_mmf( \
    char *mmf_file, uint64_t *nrows, uint64_t *ncols, uint64_t *nnz, \
    void *metadata); \
\
void spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_make_map(void *spm); \
\
uint64_t spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_size(void *spm); \
\
uint64_t spm_crs ## __idx_bits ## _ ## __elem_type ## \
         _sym_mt_map_size(void *spm); \
\
spmv_ ## __elem_type ## _sym_fn_t \
    spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_multiply;

SPM_CRS_SYM_MT_DECLARE(32, double)
SPM_CRS_SYM_MT_DECLARE(64, double)
SPM_CRS_SYM_MT_DECLARE(32, float)
SPM_CRS_SYM_MT_DECLARE(64, float)

#ifdef SPM_NUMA
#define SPM_CRS_SYM_MT_NUMA_DECLARE(__idx_bits, __elem_type) \
\
struct spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_numa { \
	spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_t *crs; \
	uint64_t row_start, row_end, nnz; \
}; \
typedef struct spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_numa \
               spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_numa_t; \
\
void * spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_numa_init_mmf( \
    char *mmf_file, uint64_t *nrows, uint64_t *ncols, uint64_t *nnz, \
    void *metadata); \
\
void spm_crs ## __idx_bits ## _ ## __elem_type ## \
         _sym_mt_numa_make_map(void *spm); \
\
uint64_t spm_crs ## __idx_bits ## _ ## __elem_type ## \
             _sym_mt_numa_size(void *spm); \
\
uint64_t spm_crs ## __idx_bits ## _ ## __elem_type ## \
             _sym_mt_numa_map_size(void *spm); \
\
spmv_ ## __elem_type ## _sym_fn_t \
    spm_crs ## __idx_bits ## _ ## __elem_type ## _sym_mt_numa_multiply;

SPM_CRS_SYM_MT_NUMA_DECLARE(32, double)
SPM_CRS_SYM_MT_NUMA_DECLARE(64, double)
SPM_CRS_SYM_MT_NUMA_DECLARE(32, float)
SPM_CRS_SYM_MT_NUMA_DECLARE(64, float)

#endif /* SPM_NUMA */

#define SPM_CRS_SYM_MT_NAME(name) CON6(spm_crs, SPM_CRS_BITS, _, ELEM_TYPE, \
                                       _sym_mt, name)
#define SPM_CRS_SYM_MT_TYPE SPM_CRS_SYM_MT_NAME(_t)
#define SPM_CRS_SYM_MT_MAP_NAME(name) CON5(map, SPM_CRS_BITS, _, ELEM_TYPE, \
                                           name)
#define SPM_CRS_SYM_MT_MAP_TYPE SPM_CRS_SYM_MT_MAP_NAME(_t)

#endif /* __SPM_CRS_SYM_MT_H__ */
