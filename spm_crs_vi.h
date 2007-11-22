#ifndef __SPM_CRS_H__
#define __SPM_CRS_H__

#include <inttypes.h>

#define SPM_CRS_VI_VIDX_TYPE uint8_t

typedef struct {
	double                 *values;
	uint64_t               *col_ind, *row_ptr;
	SPM_CRS_VI_VIDX_TYPE   *val_ind;
	uint64_t               nz, nrows, ncols;
} spm_crsvi64_double_t;

typedef struct {
	float                  *values;
	uint64_t               *col_ind, *row_ptr;
	SPM_CRS_VI_VIDX_TYPE   *val_ind;
	uint64_t               nz, nrows, ncols;
} spm_crsvi64_float_t;


typedef struct {
	double                 *values;
	uint32_t               *col_ind, *row_ptr;
	SPM_CRS_VI_VIDX_TYPE   *val_ind;
	uint64_t               nz, nrows, ncols;
} spm_crsvi32_double_t;

typedef struct {
	float                  *values;
	uint32_t               *col_ind, *row_ptr;
	SPM_CRS_VI_VIDX_TYPE   *val_ind;
	uint64_t               nz, nrows, ncols;
} spm_crsvi32_float_t;

spm_crsvi32_double_t *
spm_crsvi32_double_init_mmf(char *mmf_file, 
                          unsigned long *rows_nr, unsigned long *cols_nr,
                          unsigned long *nz_nr);

spm_crsvi32_float_t *
spm_crsvi32_float_init_mmf(char *mmf_file, 
                         unsigned long *rows_nr, unsigned long *cols_nr,
                         unsigned long *nz_nr);

spm_crsvi64_double_t *
spm_crsvi64_double_init_mmf(char *mmf_file, 
                          unsigned long *rows_nr, unsigned long *cols_nr,
			  unsigned long *nz_nr);

spm_crsvi64_float_t *
spm_crsvi64_float_init_mmf(char *mmf_file, 
                         unsigned long *rows_nr, unsigned long *cols_nr,
                         unsigned long *nz_nr);

void spm_crsvi32_double_destroy(spm_crsvi32_double_t *crsvi);
void spm_crsvi32_float_destroy(spm_crsvi32_float_t *crsvi);
void spm_crsvi64_double_destroy(spm_crsvi64_double_t *crsvi);
void spm_crsvi64_float_destroy(spm_crsvi64_float_t *crsvi);


#include "spmv_method.h"
spmv_double_fn_t spm_crsvi32_double_multiply;
spmv_double_fn_t spm_crsvi64_double_multiply;
spmv_float_fn_t spm_crsvi32_float_multiply;
spmv_float_fn_t spm_crsvi64_float_multiply;

#endif
