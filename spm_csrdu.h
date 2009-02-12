#ifndef SPM_CSRDU_H__
#define SPM_CSRDU_H__

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "spm_mt.h"

#define _NAME(val_type, name) spm_csrdu_ ## val_type ## name
#define _TYPE(val_type) _NAME(val_type, _t)

#define SPM_CSRDU_DECLARE(val_type) \
typedef struct { \
	uint64_t        nnz, ncols, nrows, ctl_size; \
	uint8_t         *ctl; \
	val_type	*values; \
} _TYPE(val_type); \
\
void *\
_NAME(val_type,_init_mmf)(char *mmf, \
                          uint64_t *nrows, uint64_t *ncols, uint64_t *nnz); \
void * \
_NAME(val_type,_mt_init_mmf)(char *mmf, \
                             uint64_t *nrows, uint64_t *ncols, uint64_t *nnz); \
uint64_t _NAME(val_type, _size)(void *csrdu); \
uint64_t _NAME(val_type, _mt_size)(void *csrdu); \
void _NAME(val_type, _destroy)(void *csrdu); \
void _NAME(val_type, _mt_destroy)(void *csrdu);

SPM_CSRDU_DECLARE(float)
SPM_CSRDU_DECLARE(double)

#undef _NAME
#undef _TYPE

#include "macros.h"
#define SPM_CSRDU_NAME(name) CON3(spm_csrdu_, ELEM_TYPE, name)
#define SPM_CSRDU_TYPE SPM_CSRDU_NAME(_t)

/*
 * CSR-DU storage format (modified delta storage format)
 *
 * ->values: array of values
 * ->ctl   : control data
 *
 * control data is a buffer for placing control information for the spmv
 * operation. the minimum data grain in the ctl structure is a char.  The
 * general idea is that the sparse matrix can be divided to units. Each of
 * these units has different structure, which can be exploited to achieve
 * better spmv performance by applying various optimizations.
 *
 */

/*
 * CSR-DU header: each unit starts with a header of two bytes
 *  - flags:
 *    bits 0-1 : size of delta column indices when sparse (0:u8..3:u64)
 *    bit 6    : sparse flag
 *    bit 7    : new row flag
 *
 *  - size     : unit size
 */
#define SPM_CSRDU_FL_NR_SHIFT (7) // new row flag
#define SPM_CSRDU_FL_SP_SHIFT (6) // sparse flag
#define SPM_CSRDU_FL_CISIZE_MASK (3) // delta size mask

// unit selection mask
#define SPM_CSRDU_FL_UNIT_MASK (SPM_CSRDU_FL_CISIZE_MASK | (1<<SPM_CSRDU_FL_SP_SHIFT))

#define SPM_CSRDU_SIZE_MAX  ((1<<8) - 1)

// new row (NR) flag functions
static inline void spm_csrdu_fl_setnr(uint8_t *ctl)
{
	*ctl |= (1<<SPM_CSRDU_FL_NR_SHIFT);
}
static inline void spm_csrdu_fl_clearnr(uint8_t *ctl)
{
	*ctl &= ~(1<<SPM_CSRDU_FL_NR_SHIFT);
}
static inline int spm_csrdu_fl_isnr(uint8_t ctl)
{
	return (ctl & (1<<SPM_CSRDU_FL_NR_SHIFT));
}


// sparse (SP) flag functions
static inline void spm_csrdu_fl_setsp(uint8_t *ctl)
{
	*ctl |= (1<<SPM_CSRDU_FL_SP_SHIFT);
}
static inline void spm_csrdu_fl_clearsp(uint8_t *ctl)
{
	*ctl &= ~(1<<SPM_CSRDU_FL_SP_SHIFT);
}
static inline int spm_csrdu_fl_issp(uint8_t *ctl)
{
	return (*ctl & (1<<SPM_CSRDU_FL_SP_SHIFT));
}

// delta size(s) for column indices
enum {
	SPM_CSRDU_CISIZE_U8  = 0,
	SPM_CSRDU_CISIZE_U16 = 1,
	SPM_CSRDU_CISIZE_U32 = 2,
	SPM_CSRDU_CISIZE_U64 = 3,
	SPM_CSRDU_CISIZE_NR  = 4 // number of cis
};

#define SPM_CSRDU_FL_UNIT_DENSE 0
#define SPM_CSRDU_FL_UNIT_SP_U8 ((1<<SPM_CSRDU_FL_SP_SHIFT) + 0)
#define SPM_CSRDU_FL_UNIT_SP_U16 ((1<<SPM_CSRDU_FL_SP_SHIFT) + 1)
#define SPM_CSRDU_FL_UNIT_SP_U32 ((1<<SPM_CSRDU_FL_SP_SHIFT) + 2)
#define SPM_CSRDU_FL_UNIT_SP_U64 ((1<<SPM_CSRDU_FL_SP_SHIFT) + 3)

static inline void spm_csrdu_fl_setcisize(uint8_t *ctl, uint8_t size)
{
	*ctl |= (size & SPM_CSRDU_FL_CISIZE_MASK);
}
static inline uint8_t spm_csrdu_fl_getcisize(uint8_t *ctl)
{
	return (*ctl & SPM_CSRDU_FL_CISIZE_MASK);
}

// number of bytes required for specified cisize
static inline int spm_csrdu_cisize_bytes(uint8_t cisize)
{
	switch (cisize){
		case SPM_CSRDU_CISIZE_U8:
		case SPM_CSRDU_CISIZE_U16:
		case SPM_CSRDU_CISIZE_U32:
		case SPM_CSRDU_CISIZE_U64:
		return (1<<cisize);
	}

	fprintf(stderr, "error in spm_csrdu_cisize_bytes: cisize=%u \n", cisize);
	exit(1);
}

#define define_copy(fn, type)                      \
static inline void _copy##fn(void *dst,            \
                             uint64_t *src,        \
                             unsigned long items)  \
{                                                  \
	type *_dst = dst;                          \
	unsigned long i;                           \
                                                   \
	for ( i=0; i < items; i++){                \
		_dst[i] = (type)src[i];            \
	}                                          \
}

define_copy(_u8,  uint8_t)
define_copy(_u16, uint16_t)
define_copy(_u32, uint32_t)
define_copy(_u64, uint64_t)

#undef define_copy

static inline void spm_csrdu_cisize_copy(void *dst, uint64_t *src,
                                         int items, uint8_t cisize)
{
	switch (cisize){
		case SPM_CSRDU_CISIZE_U8:
			_copy_u8(dst,src,items);
			break;
		case SPM_CSRDU_CISIZE_U16:
			_copy_u16(dst,src,items);
			break;
		case SPM_CSRDU_CISIZE_U32:
			_copy_u32(dst,src,items);
			break;
		case SPM_CSRDU_CISIZE_U64:
			_copy_u64(dst,src,items);
			break;
		default:
			fprintf(stderr, "error in cisize_copy: cisize=%u\n", cisize);
			exit(1);
	}
}

// unaligned copy
#define UINT_UA(bits) ua_u ## bits ## _t
#define define_ua_copy(bits)                       \
typedef struct {  \
	UINT_TYPE(bits) val __attribute__(( packed ));  \
} UINT_UA(bits); \
\
static void _copy_ua_u##bits(void *dst,               \
                          uint64_t *src,              \
                          unsigned long items)        \
{                                                     \
	UINT_UA(bits) *_dst = dst;                    \
	unsigned long i;                              \
                                                      \
	for ( i=0; i < items; i++){                   \
		_dst->val = (UINT_TYPE(bits))src[i];  \
		_dst++;                               \
	}                                             \
}

define_ua_copy(16)
define_ua_copy(32)
define_ua_copy(64)

static inline void spm_csrdu_cisize_copy_ua(void *dst, uint64_t *src,
                                            int items, uint8_t cisize)
{
	switch (cisize){
		case SPM_CSRDU_CISIZE_U8:
			_copy_u8(dst,src,items);
			break;
		case SPM_CSRDU_CISIZE_U16:
			_copy_ua_u16(dst,src,items);
			break;
		case SPM_CSRDU_CISIZE_U32:
			_copy_ua_u32(dst,src,items);
			break;
		case SPM_CSRDU_CISIZE_U64:
			_copy_ua_u64(dst,src,items);
			break;
		default:
			fprintf(stderr, "error in cisize_copy: cisize=%u\n", cisize);
			exit(1);
	}
}
#undef define_ua_copy

#define  spm_csrdu_ucmax(ci_size) (1UL<<(8<<ci_size))
static inline uint8_t spm_csrdu_usize(uint64_t u)
{
	if ( u < spm_csrdu_ucmax(SPM_CSRDU_CISIZE_U8)){
		return SPM_CSRDU_CISIZE_U8;
	}

	if ( u < spm_csrdu_ucmax(SPM_CSRDU_CISIZE_U16)){
		return SPM_CSRDU_CISIZE_U16;
	}

	if ( u < spm_csrdu_ucmax(SPM_CSRDU_CISIZE_U32)){
		return SPM_CSRDU_CISIZE_U32;
	}

	return SPM_CSRDU_CISIZE_U64;
	fprintf(stderr, "error in spm_delta_usize: delta=%lu\n", u);
	exit(1);
}
static inline uint8_t spm_csrdu_cisize(unsigned long delta)
{
	return spm_csrdu_usize(delta);
}



#define u8_get(ptr) ({ uint8_t *_ptr = ptr;  ptr++; *_ptr;})

#define u16_get(ptr) ({                                \
	uint16_t ret = *((uint16_t *)ptr);             \
	ptr += sizeof(uint16_t);                       \
	ret;                                           \
})

#define u32_get(ptr) ({                             \
	uint32_t ret = *((uint32_t *)ptr);          \
	ptr += sizeof(uint32_t);                    \
	ret;                                        \
})

#define u64_get(ptr) ({                               \
	uint64_t ret = *((uint64_t *)ptr);            \
	ptr += sizeof(uint64_t);                      \
	ret;                                          \
})

#define uc_get_ul(ptr)                                       \
({                                                           \
	unsigned long _val;                                  \
	                                                     \
	_val = u8_get(ptr);                                  \
	if ( _val > 127 ){                                   \
		unsigned shift = 7;                          \
		unsigned long _uc;                           \
		_val -= 128;                                 \
		for (;;){                                    \
			_uc = u8_get(ptr);                   \
			if ( _uc > 127 ){                    \
				_uc -= 128;                  \
				_val += (_uc<<shift);        \
				shift += 7;                  \
			} else {                             \
				_val += (_uc<<shift);        \
				break;                       \
			}                                    \
		}                                            \
	}                                                    \
	_val;                                                \
})                                                           \

#define u8_get_ul(ptr) uc_get_ul(ptr)

#include "dynarray.h"
#define LONGUC_SHIFT (7)
static inline void da_uc_put_ul(dynarray_t *da, unsigned long val)
{
	unsigned char *uc;
	const unsigned shift = LONGUC_SHIFT;

	for (;;) {
		uc = dynarray_alloc(da);
		*uc = (val & ((1<<shift) - 1));
		if ( val < (1<<shift) ){
			break;
		}
		*uc |= (1<<shift);
		val >>= shift;
	}
}

/*
 * Multithreading Support
 */
typedef struct {
	void *csrdu;
	uint64_t nnz, ctl_start, row_start, val_start;
} spm_csrdu_mt_t;

#endif /* SPM_CSRDU_H__ */
