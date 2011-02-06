#ifndef SPM_BCSR_H__
#define SPM_BCSR_H__

#include "spm_crs.h"
#include <inttypes.h>

#define BCSR_NAME(val_type, ci_bits, name) spm_bcsr ## ci_bits ## _ ## val_type ## name
#define BCSR_TYPE(val_type, ci_bits) BCSR_NAME(val_type, ci_bits, _t)
#define BCSR_MT_TYPE(val_type, ci_bits) BCSR_NAME(val_type, ci_bits, _mt_t)

#define BCSR_DECLARE(val_type, ci_bits)    \
typedef struct {                           \
	uint64_t nr_rows, nr_cols, nr_nzeros;  \
	uint64_t nr_blocks, nr_brows;          \
	int storage;                           \
	uint32_t br, bc;                       \
	val_type *bvalues;                     \
	UINT_TYPE(ci_bits) *brow_ptr;          \
	UINT_TYPE(ci_bits) *bcol_ind;          \
} BCSR_TYPE(val_type, ci_bits) ;           \
                                           \
typedef struct {                           \
	BCSR_TYPE(val_type, ci_bits) *bcsr;    \
	uint64_t row_start, row_end;           \
	uint64_t nnz_nr;                       \
} BCSR_MT_TYPE(val_type, ci_bits);

BCSR_DECLARE(double, 32)
BCSR_DECLARE(double, 64)
BCSR_DECLARE(float, 32)
BCSR_DECLARE(float, 64)

#undef _NAME
#undef _TYPE

#include "macros.h"

#define SPM_BCSR_NAME(name) CON5(spm_bcsr, SPM_CRS_BITS, _,ELEM_TYPE,name)
#define SPM_BCSR_TYPE SPM_BCSR_NAME(_t)

#define SPM_BCSR_MT_NAME(name) CON6(spm_bcsr, SPM_CRS_BITS, _,ELEM_TYPE,_mt,name)
#define SPM_BCSR_MT_TYPE SPM_BCSR_MT_NAME(_t)

#endif /* SPM_BCSR_H__ */
