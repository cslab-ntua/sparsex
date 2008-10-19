#ifndef __SPM_CRS_VH_MT_H__
#define __SPM_CRS_VH_MT_H__

#include <inttypes.h>
#include "spmv_method.h"
#include "spm_mt.h"

#include "huffman.h"

//#define SPM_CRS_VI_VIDX_TYPE uint8_t

#define UINT_TYPE(bits) uint ## bits ## _t
#define _NAME(val_type, ci_bits, name) \
	spm_crs ## ci_bits ## _vh_ ## val_type ## _mt ## name

#define DECLARE_CRS_VH_MT(val_type, ci_bits)              \
	typedef struct {                                  \
		UINT_TYPE(ci_bits)  *col_ind, *row_ptr;   \
		uint64_t            nz, nrows, ncols;     \
		unsigned char       *hvals;               \
		unsigned long       row_start, row_end;   \
		unsigned long       hvals_bits;           \
		huff_node_t         *htree;               \
	} _NAME(val_type, ci_bits, _t);          \
	                                                  \
	spm_mt_t *           \
	_NAME(val_type, ci_bits, _init_mmf)(     \
		char *mmf_file,                           \
		unsigned long *rows_nr,                   \
		unsigned long *cols_nr,                   \
		unsigned long *nz_nr                      \
	);                                                \
	                                                  \
	void                                              \
	_NAME(val_type, ci_bits, _destroy)(      \
		_NAME(val_type, ci_bits, _t) *m  \
	);                                                \
	                                                  \
	unsigned long                                     \
	_NAME(val_type, ci_bits, _size)(         \
		spm_mt_t *                       \
	);                                                \

DECLARE_CRS_VH_MT(double, 32)
DECLARE_CRS_VH_MT(float,  32)
DECLARE_CRS_VH_MT(double, 64)
DECLARE_CRS_VH_MT(float,  64)

#undef UNIT_TYPE
#undef _NAME
#undef DECLARE_CRS_VI

#endif
