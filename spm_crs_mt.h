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
	uint64_t nnz_nr; \
}; \
typedef struct spm_crs ## __idx_bits ## _ ## __elem_type ## _mt spm_crs ## __idx_bits ## _ ## __elem_type ## _mt ## _t; \
\
spm_mt_t * \
spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_init_mmf( \
	char *mmf_file, \
	unsigned long *rows_nr, unsigned long *cols_nr, \
	unsigned long *nz_nr); \
\
/* XXX: Destroy */ \
\
spmv_ ## __elem_type ## _fn_t spm_crs ## __idx_bits ## _ ## __elem_type ## _mt_multiply;

SPM_CRS_MT_DECLARE(32, double)
SPM_CRS_MT_DECLARE(64, double)
SPM_CRS_MT_DECLARE(32, float)
SPM_CRS_MT_DECLARE(64, float)


#endif
