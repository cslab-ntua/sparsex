#ifndef __SPM_DELTA_H__
#define __SPM_DELTA_H__

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "vector.h"
#include "spmv_method.h"

typedef struct {
	unsigned long nnz, ncols, nrows;
	unsigned char *ctl;
	float         *values;
} spm_delta_float_t;
spmv_float_fn_t spm_delta_float_multiply;
spmv_float_fn_t spm_delta_float_jmp_multiply;

typedef struct {
	unsigned long  nnz, ncols, nrows;
	unsigned char  *ctl;
	double         *values;
} spm_delta_double_t;
spmv_double_fn_t spm_delta_double_multiply;
spmv_double_fn_t spm_delta_double_jmp_multiply;

#define SPM_DELTA_CTL_SIZE_MAX    ((1<<8) - 1)

/*
 * delta storage format
 * 
 * The delta storage format consists of two arrays:
 *  - delta->values array of values
 *  - delta->ctl array of control data
 *
 * control data is a (data) buffer for placing control information for the spmv
 * operation. the minimum data grain in the ctl structure is 8-bytes (unsigned
 * char). The general concept is that each spmv operation can be divided to
 * "units". Each of this units has different structure, which can be exploited
 * to achieve better performance by applying various optimizations. 
 * 
 * flags bits:
 * 7   : nr      : new row flag
 * 6   : sp      : sparse flag
 * 5   : nr_seq  : new row is sequential
 *
 * 0-1 : size of colind (sp must be set)
 *   3   : u64
 *   2   : u32
 *   1   : u16
 *   0   : u8
 */
#define SPM_DELTA_FL_NR_SHIFT     (7)
#define SPM_DELTA_FL_SP_SHIFT     (6)
#define SPM_DELTA_FL_NRSEQ_SHIFT  (5)
#define SPM_DELTA_FL_CISIZE_MASK  (3)


static inline void spm_delta_fl_setnr(unsigned char *ctl)
{
	*ctl |= (1<<SPM_DELTA_FL_NR_SHIFT);
}

static inline int spm_delta_fl_isnr(unsigned char ctl)
{
	return (ctl & (1<<SPM_DELTA_FL_NR_SHIFT) );
}

static inline void spm_delta_fl_setsp(unsigned char *ctl)
{
	*ctl |= (1<<SPM_DELTA_FL_SP_SHIFT);
}

static inline int spm_delta_fl_issp(unsigned char ctl)
{
	return (ctl & (1<<SPM_DELTA_FL_SP_SHIFT) );
}

static inline void spm_delta_fl_setnrseq(unsigned char *ctl)
{
	*ctl |= (1<<SPM_DELTA_FL_NRSEQ_SHIFT);
}

static inline int spm_delta_fl_isnrseq(unsigned char ctl)
{
	return (ctl & (1<<SPM_DELTA_FL_NRSEQ_SHIFT) );
}

enum {
	SPM_DELTA_CISIZE_U8  = 0,
	SPM_DELTA_CISIZE_U16 = 1,
	SPM_DELTA_CISIZE_U32 = 2,
	SPM_DELTA_CISIZE_U64 = 3
};


#define SPM_DELTA_FL_UNIT_MASK (~((1<<SPM_DELTA_FL_NR_SHIFT)+(1<<SPM_DELTA_FL_NRSEQ_SHIFT)))

#define SPM_DELTA_FL_UNIT_DENSE 0
#define SPM_DELTA_FL_UNIT_SP_UC ((1<<SPM_DELTA_FL_SP_SHIFT) + 0)
#define SPM_DELTA_FL_UNIT_SP_US ((1<<SPM_DELTA_FL_SP_SHIFT) + 1)
#define SPM_DELTA_FL_UNIT_SP_UI ((1<<SPM_DELTA_FL_SP_SHIFT) + 2)
#define SPM_DELTA_FL_UNIT_SP_UL ((1<<SPM_DELTA_FL_SP_SHIFT) + 3)

#define  spm_delta_u_max(ci_size) ( 1UL<<(8<<ci_size) )
static inline unsigned char spm_delta_usize(unsigned long delta)
{
	if ( delta < spm_delta_u_max(SPM_DELTA_CISIZE_U8)){
		return SPM_DELTA_CISIZE_U8;
	}

	if ( delta < spm_delta_u_max(SPM_DELTA_CISIZE_U16)){
		return SPM_DELTA_CISIZE_U16;
	}

	if ( delta < spm_delta_u_max(SPM_DELTA_CISIZE_U32)){
		return SPM_DELTA_CISIZE_U32;
	}

	return SPM_DELTA_CISIZE_U64;

	fprintf(stderr, "error in spm_delta_isize: delta=%lu\n", delta);
	exit(1);
}

#if 0
#define  spm_delta_s_max(si_size) ( 1L<<((8<<si_size)-1) )
static inline unsigned char spm_delta_ssize(long delta)
{
	if ( delta < 0){
		delta = -delta;
	}
	if (delta < spm_delta_s_max(SPM_DELTA_CISIZE_U8)){
		return SPM_DELTA_CISIZE_U8;
	}

	if (delta < spm_delta_s_max(SPM_DELTA_CISIZE_U16)){
		return SPM_DELTA_CISIZE_U16;
	}

	if (delta < spm_delta_s_max(SPM_DELTA_CISIZE_U32)){
		return SPM_DELTA_CISIZE_U32;
	}

	return SPM_DELTA_CISIZE_U64;

	fprintf(stderr, "error in spm_delta_ssize: delta=%lu\n", delta);
	exit(1);
}
#endif

static inline unsigned char spm_delta_cisize(unsigned long delta)
{
	return spm_delta_usize(delta);
}

static inline unsigned char spm_delta_visize(unsigned long delta)
{
	return spm_delta_usize(delta);
}

static inline unsigned long spm_delta_cisize_size(unsigned char cisize)
{
	switch (cisize){
		case SPM_DELTA_CISIZE_U8:
		case SPM_DELTA_CISIZE_U16:
		case SPM_DELTA_CISIZE_U32:
		case SPM_DELTA_CISIZE_U64:
		return (1UL<<cisize);
	}

	fprintf(stderr, "error in spm_delta_cisize_size: cisize=%u \n", cisize);
	exit(1);
}

static inline void spm_delta_fl_setcisize(unsigned char *ctl, unsigned char size)
{
	*ctl |= (size & SPM_DELTA_FL_CISIZE_MASK);
}

static inline unsigned char spm_delta_fl_getcisize(unsigned char ctl)
{
	return (ctl & SPM_DELTA_FL_CISIZE_MASK);
}

#define uc_get(ptr) ({ unsigned char *_ptr = ptr;  ptr++; *_ptr;})

#define u8_get(ptr) ({ uint8_t *_ptr = ptr;  ptr++; *_ptr;})

#define us_get(ptr) ({                                 \
	unsigned short ret = *((unsigned short *)ptr); \
	ptr += sizeof(unsigned short);                 \
	ret;                                           \
})

#define u16_get(ptr) ({                                \
	uint16_t ret = *((uint16_t *)ptr);             \
	ptr += sizeof(uint16_t);                       \
	ret;                                           \
})

#define ui_get(ptr) ({                             \
	unsigned int ret = *((unsigned int *)ptr); \
	ptr += sizeof(unsigned int);               \
	ret;                                       \
})

#define u32_get(ptr) ({                             \
	uint32_t ret = *((uint32_t *)ptr);          \
	ptr += sizeof(uint32_t);                    \
	ret;                                        \
})

#define ul_get(ptr) ({                               \
	unsigned long ret = *((unsigned long *)ptr); \
	ptr += sizeof(unsigned long);                \
	ret;                                         \
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
	_val = uc_get(ptr);                                  \
	if ( _val > 127 ){                                   \
		unsigned shift = 7;                          \
		unsigned long _uc;                           \
		_val -= 128;                                 \
		for (;;){                                    \
			_uc = uc_get(ptr);                   \
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


#define define_copy(fn, type)                      \
static inline void _copy##fn(void *dst,            \
                             unsigned long *src,   \
                             unsigned long items)  \
{                                                  \
	type *_dst = dst;                          \
	unsigned long i;                           \
                                                   \
	for ( i=0; i < items; i++){                \
		_dst[i] = (type)src[i];            \
	}                                          \
}


define_copy(_uc, unsigned char)
define_copy(_us, unsigned short)
define_copy(_ui, unsigned int)
define_copy(_ul, unsigned long)


static inline void spm_delta_cisize_copy(void *dst, 
                                         unsigned long *src, 
			                 unsigned long items, 
			                 unsigned char cisize)
{
	switch (cisize){
		case SPM_DELTA_CISIZE_U8:
			_copy_uc(dst,src,items);
			break;
		case SPM_DELTA_CISIZE_U16:
			_copy_us(dst,src,items);
			break;
		case SPM_DELTA_CISIZE_U32:
			_copy_ui(dst,src,items);
			break;
		case SPM_DELTA_CISIZE_U64:
			_copy_ul(dst,src,items);
			break;
		default:
			fprintf(stderr, "error in copy_cisize: cisize=%u \n", cisize);
			exit(1);
	}
}


static inline void spm_delta_cisize_set(void *dst, unsigned long val, int cisize)
{
	switch (cisize){
		case SPM_DELTA_CISIZE_U8:
			*((uint8_t *)dst) = (uint8_t)val;
			break;
		case SPM_DELTA_CISIZE_U16:
			*((uint16_t *)dst) = (uint16_t)val;
			break;
		case SPM_DELTA_CISIZE_U32:
			*((uint32_t *)dst) = (uint32_t)val;
			break;
		case SPM_DELTA_CISIZE_U64:
			*((uint64_t *)dst) = (uint64_t)val;
			break;
		default:
			fprintf(stderr, "error in copy_cisize: cisize=%u \n", cisize);
			exit(1);
	}
}

/*
spm_delta_t *spm_delta_create(unsigned long rows_nr, unsigned long cols_nr,
                               unsigned long nz_nr);
void spm_delta_ctl_print(unsigned char flags, unsigned int size);
spm_delta_t *spm_delta_init_mmf(char *mmf_file, 
                                unsigned long *rows_nr, unsigned long *cols_nr,
                                unsigned long *nz_nr);

*/

#endif /* __SPM_DELTA_H__ */
