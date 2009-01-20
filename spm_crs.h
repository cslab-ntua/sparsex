#ifndef __SPM_CRS_H__
#define __SPM_CRS_H__

#include <inttypes.h>
#include "spmv_method.h"

#define _NAME(val_type, ci_bits, name) spm_crs ## ci_bits ## _ ## val_type ## name
#define _TYPE(val_type, ci_bits) _NAME(val_type, ci_bits, _t)

#define SPM_DECLARE(val_type, ci_bits) \
typedef struct { \
	val_type            *values; \
	UINT_TYPE(ci_bits)  *col_ind, *row_ptr; \
	uint64_t            nz, nrows, ncols; \
} _TYPE(val_type, ci_bits); \
\
_TYPE(val_type, ci_bits)  *\
_NAME(val_type, ci_bits, _init_mf)(char *mmf_file, \
                                    uint64_t *rows_nr, uint64_t *cols_nr, \
				    uint64_t *nz_nr); \
void _NAME(val_type, ci_bits, _destroy)(_TYPE(val_type, ci_bits) *crs); \
spmv_  ## val_type ## _fn_t _NAME(val_type, ci_bits, _multiply);

SPM_DECLARE(double, 32)
SPM_DECLARE(double, 64)
SPM_DECLARE(float, 32)
SPM_DECLARE(float, 64)

#undef _NAME
#undef _TYPE

#include "macros.h"
#define SPM_CRS_NAME(name) CON5(spm_crs, SPM_CRS_BITS, _, ELEM_TYPE, name)
#define SPM_CRS_TYPE SPM_CRS_NAME(_t)
#define SPM_CRS_IDX_TYPE UINT_TYPE(SPM_CRS_BITS)

#endif
