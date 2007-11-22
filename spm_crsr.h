#ifndef __SPM_CRSR__
#define __SPM_CRSR__

#include "matrix.h"
#include "spm_crs.h"

struct spm_crsr_unit {
	unsigned long   val_indx;
	unsigned long   col_indx;
	unsigned long   len;
};
typedef struct spm_crsr_unit spm_crsr_unit_t;


struct spm_crsr_row {
	spm_crsr_unit_t  *units;
	unsigned long    units_nr;
};
typedef struct spm_crsr_row spm_crsr_row_t;

struct spm_crsr_s {
	elem_t          *values;
	spm_crsr_row_t  *rows;
	unsigned long   nnz, ncols, nrows;
};
typedef struct spm_crsr_s spm_crsr_t;

struct spm_hybrid_s {
	spm_crsr_t    *crsr;
	spm_crs64_t     *crs;
};
typedef struct spm_hybrid_s spm_hybrid_t;


struct spm_crsr2_row {
	spm_crsr_unit_t  *units;
	unsigned long    units_nr;
	unsigned long    row_indx;
};
typedef struct spm_crsr2_row spm_crsr2_row_t;

struct spm_crsr2_s {
	elem_t           *values;
	spm_crsr2_row_t  *rows;
	unsigned long    nnz, ncols, nrows;
};
typedef struct spm_crsr2_s spm_crsr2_t;

struct spm_hybrid2_s {
	spm_crsr2_t *crsr2;
	spm_crs64_t *crs;
};
typedef struct spm_hybrid2_s spm_hybrid2_t;

#endif
