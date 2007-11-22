#ifndef __SPM_DELTA_CV_H__
#define __SPM_DELTA_CV_H__

#include <assert.h>
#include <stdint.h>

#include "spmv_method.h"
#include "spm_delta.h"


typedef struct {
	unsigned long nnz, ncols, nrows;
	unsigned long ctl_size, ci_size, vi_size, v_size;
	unsigned char *ctl;
	unsigned char *c_indices;
	unsigned char *v_indices;
	float         *values;
} spm_delta_cv_float_t;
spmv_float_fn_t spm_delta_cv_float_multiply;
spmv_float_fn_t spm_delta_cv_float_jmp_multiply;

typedef struct {
	unsigned long  nnz, ncols, nrows;
	unsigned long ctl_size, ci_size, vi_size, v_size;
	unsigned char  *ctl;
	unsigned char  *c_indices;
	unsigned char  *v_indices;
	double         *values;
} spm_delta_cv_double_t;
spmv_double_fn_t spm_delta_cv_double_multiply;
spmv_double_fn_t spm_delta_cv_double_jmp_multiply;


enum {
	SPM_DELTA_VISIZE_U8  = 0,
	SPM_DELTA_VISIZE_U16 = 1,
	SPM_DELTA_VISIZE_U32 = 2,
	SPM_DELTA_VISIZE_U64 = 3,
};


/* Additional ctl flags
 *
 * 2-3 : size of valind delta
 *   3 : u64
 *   2 : u32
 *   1 : u16
 *   0 : u8
 */
#define SPM_DELTA_FL_VISIZE_SHIFT (2)

#define SPM_DELTA_FL_CVI_SP(ci, vi) \
	(ci | (vi<<SPM_DELTA_FL_VISIZE_SHIFT) | (1<<SPM_DELTA_FL_SP_SHIFT) )

#define SPM_DELTA_FL_CVI_DE(vi) \
	(vi<<SPM_DELTA_FL_VISIZE_SHIFT)

#define SPM_DELTA_FL_CV_UMASK ( ((1<<4)-1) | (1<<SPM_DELTA_FL_SP_SHIFT) )

enum {
	SPM_DELTA_CI8_VI8   = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U8, SPM_DELTA_VISIZE_U8),
	SPM_DELTA_CI16_VI8  = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U16, SPM_DELTA_VISIZE_U8),
	SPM_DELTA_CI32_VI8  = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U32, SPM_DELTA_VISIZE_U8),
	SPM_DELTA_CI64_VI8  = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U64, SPM_DELTA_VISIZE_U8),
	SPM_DELTA_DE_VI8    = SPM_DELTA_FL_CVI_DE(SPM_DELTA_VISIZE_U8),

	SPM_DELTA_CI8_VI16  = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U8, SPM_DELTA_VISIZE_U16),
	SPM_DELTA_CI16_VI16 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U16, SPM_DELTA_VISIZE_U16),
	SPM_DELTA_CI32_VI16 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U32, SPM_DELTA_VISIZE_U16),
	SPM_DELTA_CI64_VI16 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U64, SPM_DELTA_VISIZE_U16),
	SPM_DELTA_DE_VI16   = SPM_DELTA_FL_CVI_DE(SPM_DELTA_VISIZE_U16),

	SPM_DELTA_CI8_VI32  = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U8, SPM_DELTA_VISIZE_U32),
	SPM_DELTA_CI16_VI32 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U16, SPM_DELTA_VISIZE_U32),
	SPM_DELTA_CI32_VI32 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U32, SPM_DELTA_VISIZE_U32),
	SPM_DELTA_CI64_VI32 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U64, SPM_DELTA_VISIZE_U32),
	SPM_DELTA_DE_VI32   = SPM_DELTA_FL_CVI_DE(SPM_DELTA_VISIZE_U32),

	SPM_DELTA_CI8_VI64  =  SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U8, SPM_DELTA_VISIZE_U64),
	SPM_DELTA_CI16_VI64 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U16, SPM_DELTA_VISIZE_U64),
	SPM_DELTA_CI32_VI64 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U32, SPM_DELTA_VISIZE_U64),
	SPM_DELTA_CI64_VI64 = SPM_DELTA_FL_CVI_SP(SPM_DELTA_CISIZE_U64, SPM_DELTA_VISIZE_U64),
	SPM_DELTA_DE_VI64   = SPM_DELTA_FL_CVI_DE(SPM_DELTA_VISIZE_U64)

};

static inline void spm_delta_cv_fl_set_cvi_sp(unsigned char *ctl, unsigned char ci, unsigned char vi)
{
	if (ci > SPM_DELTA_CISIZE_U64 || vi > SPM_DELTA_VISIZE_U64){
		fprintf(stderr, "error in fl_set_cvi");
		exit(1);
	}

	*ctl |= SPM_DELTA_FL_CVI_SP(ci, vi);
}

static inline void spm_delta_cv_fl_set_cvi_de(unsigned char *ctl, unsigned char vi)
{
	if (vi > SPM_DELTA_VISIZE_U64){
		fprintf(stderr, "error in fl_set_cvi");
		exit(1);
	}

	*ctl |= SPM_DELTA_FL_CVI_DE(vi);
}

#if 0
static inline void s8_set(void *ptr, int8_t val)
{
	uint8_t *u8 = ptr;
	if (val < 0){
		*u8 = (uint8_t)(1<<7);
		val = -val;
	} else {
		*u8 = (uint8_t)0;
	}

	*u8 |= (uint8_t)val;
}

#define s8_get(ptr) ({                                    \
	unsigned char _uval = u8_get(ptr);                \
	int8_t _ret;                                      \
	if ( _uval & (1<<7) ){                            \
		_ret = -((int8_t)(_uval & ~(1<<7) ));     \
	} else {                                          \
		_ret = _uval;                             \
	}                                                 \
	_ret;                                             \
})

static inline void s16_set(void *ptr, int16_t val)
{
	uint16_t *u16 = ptr;
	if (val < 0){
		*u16 = (uint16_t)(1<<15);
		val = -val;
	} else {
		*u16 = (uint16_t)0;
	}

	*u16 |= (uint16_t)val;
}

#define s16_get(ptr) ({                                  \
	uint16_t _uval = u16_get(ptr);                   \
	int16_t _ret;                                    \
	if ( _uval & (1<<15) ){                          \
		_ret = -((int16_t)(_uval & ~(1<<15) ));  \
	} else {                                         \
		_ret = _uval;                            \
	}                                                \
	_ret;                                            \
})

static inline void s32_set(void *ptr, int32_t val)
{
	uint32_t *u32 = ptr;
	if (val < 0){
		*u32 = (uint32_t)(1<<31);
		val = -val;
	} else {
		*u32 = (uint32_t)0;
	}

	*u32 |= (uint32_t)val;
}

#define s32_get(ptr) ({                                  \
	uint32_t _uval = u32_get(ptr);                   \
	int32_t _ret;                                    \
	if ( _uval & (1<<31) ){                          \
		_ret = -((int32_t)(_uval & ~(1<<31) ));  \
	} else {                                         \
		_ret = _uval;                            \
	}                                                \
	_ret;                                            \
})


static inline void si_set(void *dst, long val, int si_size)
{
	switch (si_size){
		case SPM_DELTA_CISIZE_U8:
		s8_set(dst, (int8_t)val);
		break;

		case SPM_DELTA_CISIZE_U16:
		s16_set(dst, (int16_t)val);
		break;

		case SPM_DELTA_CISIZE_U32:
		s32_set(dst, (int32_t)val);
		break;

		default:
		fprintf(stderr, "error in si_set, si_size=%d", si_size);
		exit(1);
	}
}

static inline long si_get(void *dst, int si_size)
{
	long ret;
	switch (si_size){
		case SPM_DELTA_CISIZE_U8:
		ret = s8_get(dst);
		break;

		case SPM_DELTA_CISIZE_U16:
		ret = s16_get(dst);
		break;

		case SPM_DELTA_CISIZE_U32:
		ret = s32_get(dst);
		break;

		default:
		fprintf(stderr, "error in si_set, si_size=%d", si_size);
		exit(1);
	}

	return ret;
}

#endif

#endif /* __SPM_DELTA_CV_H__ */
