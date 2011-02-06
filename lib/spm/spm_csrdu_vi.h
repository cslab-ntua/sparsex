#ifndef SPM_CSRDU_VI_H__
#define SPM_CSRDU_VI_H__

#include <inttypes.h>

#include "macros.h"
#include "spmv_method.h"

#define _NAME(val_type, vi_bits, name) spm_csrdu_vi ## vi_bits ## _ ## val_type ## name
#define _TYPE(val_type, vi_bits) _NAME(val_type, vi_bits, _t)
#define _MT_TYPE(val_type, vi_bits) _NAME(val_type, vi_bits, _mt_t)

#define SPM_CSRDU_VI_DECLARE(val_type, vi_bits)         \
	typedef struct {                                \
		uint64_t	    nnz, ncols, nrows;  \
		uint8_t             *ctl;               \
		uint64_t            ctl_size;           \
		val_type            *values;            \
		UINT_TYPE(vi_bits)  *val_ind;           \
		uint64_t            nv;                 \
	} _TYPE(val_type, vi_bits);                     \
	                                                \
	void *                                          \
	_NAME(val_type, vi_bits, _init_mmf)(            \
		char *mmf_file,                         \
		uint64_t *rows_nr,                      \
		uint64_t *cols_nr,                      \
		uint64_t *nz_nr                         \
	);                                              \
	                                                \
	void                                            \
	_NAME(val_type, vi_bits, _destroy)(             \
		void *m                                 \
	);                                              \
	                                                \
	uint64_t                                        \
	_NAME(val_type, vi_bits, _size)(                \
		void *m                                 \
	);                                              \
	                                                \
	typedef struct {                                \
		_TYPE(val_type, vi_bits) *csrdu_vi;     \
		uint64_t nnz, row_start, val_start;     \
		uint64_t ctl_start;                     \
	} _MT_TYPE(val_type, vi_bits);                  \
	                                                \
	void *                                          \
	_NAME(val_type, vi_bits, _mt_init_mmf)(         \
		char *mmf_file,                         \
		uint64_t *rows_nr,                      \
		uint64_t *cols_nr,                      \
		uint64_t *nz_nr                         \
	);                                              \
	                                                \
	void                                            \
	_NAME(val_type, vi_bits, _mt_destroy)(          \
		void *m                                 \
	);                                              \
	                                                \
	uint64_t                                        \
	_NAME(val_type, vi_bits, _mt_size)(             \
		void *m                                 \
	);                                              \



SPM_CSRDU_VI_DECLARE(float, 8)
SPM_CSRDU_VI_DECLARE(float, 16)
SPM_CSRDU_VI_DECLARE(float, 32)

SPM_CSRDU_VI_DECLARE(double, 8)
SPM_CSRDU_VI_DECLARE(double, 16)
SPM_CSRDU_VI_DECLARE(double, 32)


#ifdef SPM_NUMA
#define SPM_CSRDU_VI_NUMA_DECLARE(val_type, vi_bits)          \
	void *                                                \
	_NAME(val_type, vi_bits, _mt_numa_init_mmf)(          \
		char *mmf,                                    \
		uint64_t *nrows,                              \
		uint64_t *ncols,                              \
		uint64_t *nnz);                               \
                                                              \
	uint64_t                                              \
	_NAME(val_type, vi_bits, _mt_numa_size)(void *);      \
	                                                      \
	void                                                  \
	_NAME(val_type, vi_bits, _mt_numa_destroy)(void *);   \

SPM_CSRDU_VI_NUMA_DECLARE(float, 8)
SPM_CSRDU_VI_NUMA_DECLARE(float, 16)
SPM_CSRDU_VI_NUMA_DECLARE(float, 32)

SPM_CSRDU_VI_NUMA_DECLARE(double, 8)
SPM_CSRDU_VI_NUMA_DECLARE(double, 16)
SPM_CSRDU_VI_NUMA_DECLARE(double, 32)

#endif /* SPM_NUMA */

#undef _TYPE
#undef _MT_TYPE
#undef _NAME

#define SPM_CSRDU_VI_NAME(name) \
	CON5(spm_csrdu_vi, SPM_CRSVI_VI_BITS, _, ELEM_TYPE, name)
#define SPM_CSRDU_VI_TYPE SPM_CSRDU_VI_NAME(_t)

#define SPM_CSRDU_VI_MT_TYPE SPM_CSRDU_VI_NAME(_mt_t)

#endif /* SPM_CSRDU_VI_H__ */
