#ifndef __SPM_CRS_VI_MT_H__
#define __SPM_CRS_VI_MT_H__

#include <inttypes.h>

#define _NAME(val_type, ci_bits, vi_bits, name) \
        spm_crs ## ci_bits ## _vi ## vi_bits ## _ ## val_type ## name

#define DECLARE_CRS_VI_MT(val_type, ci_bits, vi_bits)               \
	typedef struct {                                            \
		_NAME(val_type, ci_bits, vi_bits, _t) *crsvi;       \
		uint64_t row_start, row_end, nnz;                   \
	} _NAME(val_type, ci_bits, vi_bits, _mt_t);                 \
	                                                            \
	spm_mt_t * _NAME(val_type, ci_bits, vi_bits, _mt_init_mmf)( \
		char *mmf_file,                                     \
		unsigned long *rows_nr, unsigned long *cols_nr,     \
		unsigned long *nz_nr);                              \

DECLARE_CRS_VI_MT(double, 32, 8)
DECLARE_CRS_VI_MT(double, 32, 16)
DECLARE_CRS_VI_MT(double, 32, 32)

DECLARE_CRS_VI_MT(float, 32, 8)
DECLARE_CRS_VI_MT(float, 32, 16)
DECLARE_CRS_VI_MT(float, 32, 32)

DECLARE_CRS_VI_MT(double, 64, 8)
DECLARE_CRS_VI_MT(double, 64, 16)
DECLARE_CRS_VI_MT(double, 64, 32)

DECLARE_CRS_VI_MT(float, 64, 8)
DECLARE_CRS_VI_MT(float, 64, 16)
DECLARE_CRS_VI_MT(float, 64, 32)

#undef DECLARE_CRS_VI_MT
#undef _NAME

#endif
