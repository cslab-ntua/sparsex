#ifndef __SPM_CRS_VI_H__
#define __SPM_CRS_VI_H__

#include <inttypes.h>
#include "spmv_method.h"

//#define SPM_CRS_VI_VIDX_TYPE uint8_t

#define UINT_TYPE(bits) uint ## bits ## _t
#define _NAME(val_type, ci_bits, vi_bits, name) \
	spm_crs ## ci_bits ## _vi ## vi_bits ## _ ## val_type ## name

#define DECLARE_CRS_VI(val_type, ci_bits, vi_bits)        \
	typedef struct {                                  \
		val_type            *values;              \
		UINT_TYPE(ci_bits)  *col_ind, *row_ptr;   \
		UINT_TYPE(vi_bits)  *val_ind;             \
		uint64_t            nz, nrows, ncols, nv; \
	} _NAME(val_type, ci_bits, vi_bits, _t);          \
	                                                  \
	_NAME(val_type, ci_bits, vi_bits, _t) *           \
	_NAME(val_type, ci_bits, vi_bits, _init_mmf)(     \
		char *mmf_file,                           \
		unsigned long *rows_nr,                   \
		unsigned long *cols_nr,                   \
		unsigned long *nz_nr                      \
	);                                                \
	                                                  \
	void                                              \
	_NAME(val_type, ci_bits, vi_bits, _destroy)(      \
		_NAME(val_type, ci_bits, vi_bits, _t) *m  \
	);                                                \
	                                                  \
	spmv_ ## val_type ## _fn_t                        \
	_NAME(val_type, ci_bits, vi_bits, _multiply);     \


DECLARE_CRS_VI(double, 32, 8)
DECLARE_CRS_VI(double, 32, 16)
DECLARE_CRS_VI(double, 32, 32)

DECLARE_CRS_VI(float, 32, 8)
DECLARE_CRS_VI(float, 32, 16)
DECLARE_CRS_VI(float, 32, 32)

DECLARE_CRS_VI(double, 64, 8)
DECLARE_CRS_VI(double, 64, 16)
DECLARE_CRS_VI(double, 64, 32)

DECLARE_CRS_VI(float, 64, 8)
DECLARE_CRS_VI(float, 64, 16)
DECLARE_CRS_VI(float, 64, 32)

#endif
