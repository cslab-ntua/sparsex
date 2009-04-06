#ifndef __SPM_CRS_VI_H__
#define __SPM_CRS_VI_H__

#include <inttypes.h>
#include "spmv_method.h"
#include "macros.h"

//#define SPM_CRS_VI_VIDX_TYPE uint8_t

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
	void *                                            \
	_NAME(val_type, ci_bits, vi_bits, _init_mmf)(     \
		char *mmf_file,                           \
		unsigned long *rows_nr,                   \
		unsigned long *cols_nr,                   \
		unsigned long *nz_nr                      \
	);                                                \
	                                                  \
	void                                              \
	_NAME(val_type, ci_bits, vi_bits, _destroy)(      \
		void *m                                   \
	);                                                \
	                                                  \
	uint64_t                                          \
	_NAME(val_type, ci_bits, vi_bits, _size)(         \
		void *m                                   \
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

#undef _NAME
#undef DECLARE_CRS_VI

#define _NAME(val_type, vi_bits, name) \
	spm_crsvi ## vi_bits ## _ ## val_type ## name

#define DECLARE_VI(val_type, vi_bits) \
	void _NAME(val_type, vi_bits, _vals_to_vis)( \
		val_type *vals, uint64_t vals_size, \
		UINT_TYPE(vi_bits) **uvals_idx, val_type **uvals, \
		uint64_t *uvals_size \
	);

DECLARE_VI(double, 8)
DECLARE_VI(double, 16)
DECLARE_VI(double, 32)

DECLARE_VI(float, 8)
DECLARE_VI(float, 16)
DECLARE_VI(float, 32)

#undef _NAME
#undef DECLARE_VI

#define SPM_CRSVI_CI_TYPE UINT_TYPE(SPM_CRSVI_CI_BITS)
#define SPM_CRSVI_VI_TYPE UINT_TYPE(SPM_CRSVI_VI_BITS)

#define SPM_CRSVI_NAME(name) \
	CON7(spm_crs, SPM_CRSVI_CI_BITS, _vi, SPM_CRSVI_VI_BITS, _, ELEM_TYPE, name)
#define SPM_CRSVI_TYPE SPM_CRSVI_NAME(_t)

#define SPM_VI_NAME(name) \
	CON5(spm_crsvi, SPM_CRSVI_VI_BITS, _, ELEM_TYPE, name)

#endif
