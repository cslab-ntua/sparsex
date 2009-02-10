#ifndef __SPM_CRS_VH_MT_H__
#define __SPM_CRS_VH_MT_H__

#include <inttypes.h>
#include "macros.h"
#include "spmv_method.h"
#include "spm_mt.h"

#include "huffman.h"

//#define SPM_CRS_VI_VIDX_TYPE uint8_t

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
		unsigned long       hsym_bits;            \
	} _NAME(val_type, ci_bits, _t);          \
	                                                  \
	void *           \
	_NAME(val_type, ci_bits, _init_mmf)(     \
		char *mmf_file,                           \
		unsigned long *rows_nr,                   \
		unsigned long *cols_nr,                   \
		unsigned long *nz_nr                      \
	);                                                \
	                                                  \
	void                                              \
	_NAME(val_type, ci_bits, _destroy)(               \
		void *m                                   \
	);                                                \
	                                                  \
	uint64_t                                          \
	_NAME(val_type, ci_bits, _size)(                  \
		void *m                                   \
	);                                                \

DECLARE_CRS_VH_MT(double, 32)
DECLARE_CRS_VH_MT(float,  32)
DECLARE_CRS_VH_MT(double, 64)
DECLARE_CRS_VH_MT(float,  64)

#undef _NAME
#undef DECLARE_CRS_VI

#define SPM_CRS_VH_MT_NAME(name) \
	CON6(spm_crs, SPM_CRSVH_CI_BITS, _vh_, ELEM_TYPE, _mt, name)
#define SPM_CRS_VH_MT_TYPE SPM_CRS_VH_MT_NAME(_t)

#endif
