#ifndef __SPM_CRS_H__
#define __SPM_CRS_H__

#include <inttypes.h>

typedef struct {
	double         *values;
	uint64_t       *col_ind, *row_ptr; 
	uint64_t       nz, nrows, ncols;
} spm_crs64_double_t;

typedef struct {
	float          *values;
	uint64_t       *col_ind, *row_ptr; 
	uint64_t       nz, nrows, ncols;
} spm_crs64_float_t;


typedef struct {
	double         *values;
	uint32_t       *col_ind, *row_ptr; 
	uint64_t       nz, nrows, ncols;
} spm_crs32_double_t;

typedef struct {
	float         *values;
	uint32_t      *col_ind, *row_ptr; 
	uint64_t      nz, nrows, ncols;
} spm_crs32_float_t;

spm_crs32_double_t *
spm_crs32_double_init_mmf(char *mmf_file, 
                          unsigned long *rows_nr, unsigned long *cols_nr,
                          unsigned long *nz_nr);

spm_crs32_float_t *
spm_crs32_float_init_mmf(char *mmf_file, 
                         unsigned long *rows_nr, unsigned long *cols_nr,
                         unsigned long *nz_nr);

spm_crs64_double_t *
spm_crs64_double_init_mmf(char *mmf_file, 
                          unsigned long *rows_nr, unsigned long *cols_nr,
			  unsigned long *nz_nr);

spm_crs64_float_t *
spm_crs64_float_init_mmf(char *mmf_file, 
                         unsigned long *rows_nr, unsigned long *cols_nr,
                         unsigned long *nz_nr);

void spm_crs32_double_destroy(spm_crs32_double_t *crs);
void spm_crs32_float_destroy(spm_crs32_float_t *crs);
void spm_crs64_double_destroy(spm_crs64_double_t *crs);
void spm_crs64_float_destroy(spm_crs64_float_t *crs);


#include "spmv_method.h"
spmv_double_fn_t spm_crs32_double_multiply;
spmv_double_fn_t spm_crs64_double_multiply;
spmv_float_fn_t spm_crs32_float_multiply;
spmv_float_fn_t spm_crs64_float_multiply;

#if 0
void spm_crs32_double_multiply(spm_crs32_double_t *crs, vector_double_t *x, vector_double_t *y);
void spm_crs32_float_multiply(spm_crs32_float_t *crs, vector_float_t *x, vector_float_t *y);
void spm_crs64_double_multiply(spm_crs64_double_t *crs, vector_double_t *x, vector_double_t *y);
void spm_crs64_float_multiply(spm_crs64_float_t *crs, vector_float_t *x, vector_float_t *y);
#endif

#endif
