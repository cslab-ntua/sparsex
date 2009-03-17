#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

#include "spm_csrdu.h"
#include "spmv_method.h"
#include "vector.h"

#define ALIGN(buf,a) (void *) (((unsigned long) (buf) + (a-1)) & ~(a-1))

#ifdef CSRDU_ALIGNED
#define ALIGN_PTR(ptr, align) do { ptr = ALIGN(ptr, align); } while (0)
#define _NAME(name) SPM_CSRDU_NAME(_aligned ## name)
#define _NAME_MT(name) SPM_CSRDU_NAME(_mt_aligned ## name)
#else
#define ALIGN_PTR(ptr, align) do { ; } while (0)
#define _NAME(name) SPM_CSRDU_NAME(name)
#define _NAME_MT(name) SPM_CSRDU_NAME(_mt ## name)
#endif

void static inline do_mul(uint8_t *u8, ELEM_TYPE *values, ELEM_TYPE *x, ELEM_TYPE *y, uint64_t values_size)
{
	register uint64_t y_indx = 0;
	register ELEM_TYPE yr = 0;
	register ELEM_TYPE *v = values;
	register const ELEM_TYPE *v_end = values + values_size;
	register ELEM_TYPE *myx = x;

	for (;;) {
		register uint8_t flags = u8_get(u8);
		register uint8_t size = u8_get(u8);
		//spm_delta_ctl_print(flags, size);
		if ( spm_csrdu_fl_isnr(flags) ){
			y[y_indx] = yr;
			//printf("WROTE y[%lu]=%lf\n", y_indx, y[y_indx]);
			y_indx++;
			yr = 0;
			myx = x;
		}
		//printf("y_indx: %lu x_indx: %lu\n", y_indx, (unsigned long)(myx-x));

		switch (flags & SPM_CSRDU_FL_UNIT_MASK){
			register uint64_t i;

			case SPM_CSRDU_FL_UNIT_DENSE:
			__asm__ __volatile__ (" # DE loop ");
			myx += u8_get_ul(u8);
			for (i=0; i<size; i++){
				yr += v[i] * myx[i];
				//printf("\t[DE] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v + i - values), (unsigned long)(myx + i - x), *v, myx[i], yr);
			}
			myx += (size-1);
			v += size;
			break;

			#define CSRDU_SP_SPMV(bits) \
			for (i=0; i<size; i++){ \
				myx += u ## bits ## _get(u8); \
				yr += (*v) * (*myx); \
				/* printf("\t[SP%d] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", bits, y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);  */ \
				v++;                 \
			}

			case SPM_CSRDU_FL_UNIT_SP_U8:
			CSRDU_SP_SPMV(8)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U16:
			ALIGN_PTR(u8, 2);
			CSRDU_SP_SPMV(16)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U32:
			ALIGN_PTR(u8, 4);
			CSRDU_SP_SPMV(32)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U64:
			ALIGN_PTR(u8, 8);
			CSRDU_SP_SPMV(64)
			break;

			default:
			printf("Uknown flags: %u unit=%u\n", flags, flags & SPM_CSRDU_FL_UNIT_MASK);
			assert(0);
		}

		if ( v >= v_end ){
			break;
		}
	}

	y[y_indx] = yr;
	return;
}

void static inline do_mul_jmp(uint8_t *u8, ELEM_TYPE *values, ELEM_TYPE *x, ELEM_TYPE *y, uint64_t values_size)
{
	register uint64_t y_indx = 0;
	register ELEM_TYPE yr = 0;
	register ELEM_TYPE *v = values;
	register const ELEM_TYPE *v_end = values + values_size;
	register ELEM_TYPE *myx = x;
	for (;;) {
		register uint8_t flags = u8_get(u8);
		register uint8_t size = u8_get(u8);
		//spm_delta_ctl_print(flags, size);
		if ( spm_csrdu_fl_isnr(flags) ){
			y[y_indx++] = yr;
			yr = 0;
			myx = x;
		}

		//printf("y_indx: %lu x_indx: %lu\n", y_indx, (unsigned long)(myx-x));

		myx += u8_get_ul(u8);
		switch ( flags & SPM_CSRDU_FL_UNIT_MASK ){
			register unsigned long i;

			case SPM_CSRDU_FL_UNIT_DENSE:
			__asm__ __volatile__ (" # DE loop ");
			for (i=0; i<size; i++){
				//printf("\t[DE] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += v[i] * myx[i];
			}
			myx += (size-1);
			v += size;
			break;

			#define CSRDU_JMP_SP_SPMV(bits) \
			yr += (*v++)*(*myx);            \
			for (i=1; i<size; i++){         \
				myx += u ## bits ## _get(u8);      \
				/* printf("\t[SP_U # bits] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr) */ \
				yr += (*v++) * (*myx); \
			}

			case SPM_CSRDU_FL_UNIT_SP_U8:
			CSRDU_JMP_SP_SPMV(8)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U16:
			ALIGN_PTR(u8, 2);
			CSRDU_JMP_SP_SPMV(16)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U32:
			ALIGN_PTR(u8, 4);
			CSRDU_JMP_SP_SPMV(32)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U64:
			ALIGN_PTR(u8, 8);
			CSRDU_JMP_SP_SPMV(64)
			break;

			default:
			printf("Uknown flags: %u unit=%u\n", flags, flags & SPM_CSRDU_FL_UNIT_MASK);
			assert(0);
		}

		if ( v >= v_end ){
			break;
		}
	}

	y[y_indx] = yr;
	return;
}

void _NAME(_multiply)(void *matrix, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CSRDU_TYPE *csrdu = (SPM_CSRDU_TYPE *)matrix;
	ELEM_TYPE *x = in->elements, *y = out->elements, *values = csrdu->values;
	do_mul(csrdu->ctl, values, x, y, csrdu->nnz);
}
XSPMV_METH_INIT(
 _NAME(_multiply),
 SPM_CSRDU_NAME(_init_mmf),
 SPM_CSRDU_NAME(_size),
 SPM_CSRDU_NAME(_destroy),
 sizeof(ELEM_TYPE)
)

void _NAME(_jmp_multiply)(void *matrix, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CSRDU_TYPE *csrdu = (SPM_CSRDU_TYPE *)matrix;
	ELEM_TYPE *x = in->elements, *y = out->elements, *values = csrdu->values;
	do_mul_jmp(csrdu->ctl, values, x, y, csrdu->nnz);

}
XSPMV_METH_INIT(
 _NAME(_jmp_multiply),
 SPM_CSRDU_NAME(_init_mmf),
 SPM_CSRDU_NAME(_size),
 SPM_CSRDU_NAME(_destroy),
 sizeof(ELEM_TYPE)
)

void _NAME_MT(_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	spm_csrdu_mt_t *csrdu_mt = (spm_csrdu_mt_t *)spm;
	SPM_CSRDU_TYPE *csrdu = csrdu_mt->csrdu;
	ELEM_TYPE *x = in->elements, *y = out->elements, *values = csrdu->values;

	do_mul(csrdu->ctl + csrdu_mt->ctl_start,
               values + csrdu_mt->val_start,
               x,
               y + csrdu_mt->row_start,
               csrdu_mt->nnz);
}
XSPMV_MT_METH_INIT(
 _NAME_MT(_multiply),
 SPM_CSRDU_NAME(_mt_init_mmf),
 SPM_CSRDU_NAME(_mt_size),
 SPM_CSRDU_NAME(_mt_destroy),
 sizeof(ELEM_TYPE)
)

void _NAME_MT(_jmp_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	spm_csrdu_mt_t *csrdu_mt = (spm_csrdu_mt_t *)spm;
	SPM_CSRDU_TYPE *csrdu = csrdu_mt->csrdu;
	ELEM_TYPE *x = in->elements, *y = out->elements, *values = csrdu->values;

	do_mul_jmp(csrdu->ctl + csrdu_mt->ctl_start,
                   values + csrdu_mt->val_start,
                   x,
                   y + csrdu_mt->row_start,
                   csrdu_mt->nnz);
}
XSPMV_MT_METH_INIT(
 _NAME_MT(_jmp_multiply),
 SPM_CSRDU_NAME(_mt_init_mmf),
 SPM_CSRDU_NAME(_mt_size),
 SPM_CSRDU_NAME(_mt_destroy),
 sizeof(ELEM_TYPE)
)

#ifdef SPM_NUMA

void _NAME_MT(_numa_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	spm_csrdu_mt_t *csrdu_mt = (spm_csrdu_mt_t *)spm;
	SPM_CSRDU_TYPE *csrdu = csrdu_mt->csrdu;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements + csrdu_mt->row_start;

	do_mul(csrdu->ctl, csrdu->values, x, y, csrdu_mt->nnz);
}

XSPMV_MT_METH_INIT(
 _NAME_MT(_numa_multiply),
 SPM_CSRDU_NAME(_mt_numa_init_mmf),
 SPM_CSRDU_NAME(_mt_numa_size),
 SPM_CSRDU_NAME(_mt_numa_destroy),
 sizeof(ELEM_TYPE)
)

void _NAME_MT(_jmp_numa_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	spm_csrdu_mt_t *csrdu_mt = (spm_csrdu_mt_t *)spm;
	SPM_CSRDU_TYPE *csrdu = csrdu_mt->csrdu;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements + csrdu_mt->row_start;

	do_mul_jmp(csrdu->ctl, csrdu->values, x, y, csrdu_mt->nnz);
}

XSPMV_MT_METH_INIT(
 _NAME_MT(_jmp_numa_multiply),
 SPM_CSRDU_NAME(_mt_numa_init_mmf),
 SPM_CSRDU_NAME(_mt_numa_size),
 SPM_CSRDU_NAME(_mt_numa_destroy),
 sizeof(ELEM_TYPE)
)

#endif /* SPM_NUMA */
