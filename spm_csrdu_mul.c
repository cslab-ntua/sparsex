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
#else
#define ALIGN_PTR(ptr, align) do { ; } while (0)
#define _NAME(name) SPM_CSRDU_NAME(name)
#endif

void _NAME(_multiply)(void *matrix, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CSRDU_TYPE *csrdu = (SPM_CSRDU_TYPE *)matrix;
	ELEM_TYPE *x = in->elements, *y = out->elements, *values = csrdu->values;
	register uint64_t y_indx=0;
	register ELEM_TYPE yr = 0;
	register ELEM_TYPE *v = values;
	register const ELEM_TYPE *v_end = values + csrdu->nnz;
	register ELEM_TYPE *myx = x;
	register uint8_t *u8 = csrdu->ctl;

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

		switch ( flags & SPM_CSRDU_FL_UNIT_MASK){
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

			case SPM_CSRDU_FL_UNIT_SP_U8:
			for (i=0; i<size; i++){
				//printf("\t[SP_8] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				myx += u8_get(u8);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_CSRDU_FL_UNIT_SP_U16:
			ALIGN_PTR(u8, 2);
			for (i=0; i<size; i++){
				myx += u16_get(u8);
				//printf("\t[SP_US] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_CSRDU_FL_UNIT_SP_U32:
			ALIGN_PTR(u8, 4);
			for (i=0; i<size; i++){
				myx += u32_get(u8);
				//printf("\t[SP_UI] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_CSRDU_FL_UNIT_SP_U64:
			ALIGN_PTR(u8, 8);
			for (i=0; i<size; i++){
				myx += u64_get(u8);
				//printf("\t[SP_UL] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v) * (*myx);
				v++;
			}
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
XSPMV_METH_INIT(
 _NAME(_multiply),
 SPM_CSRDU_NAME(_init_mmf),
 SPM_CSRDU_NAME(_size),
 NULL,
 sizeof(ELEM_TYPE)
)

void _NAME(_jmp_multiply)(void *matrix, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CSRDU_TYPE *csrdu = (SPM_CSRDU_TYPE *)matrix;
	ELEM_TYPE *x = in->elements, *y = out->elements, *values = csrdu->values;
	register uint64_t y_indx=0;
	register ELEM_TYPE yr = 0;
	register ELEM_TYPE *v = values;
	register const ELEM_TYPE *v_end = values + csrdu->nnz;
	register ELEM_TYPE *myx = x;
	register uint8_t *u8 = csrdu->ctl;

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

			case SPM_CSRDU_FL_UNIT_SP_U8:
			//printf("\t[SP_U8] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
			yr += (*v++)*(*myx);
			for (i=1; i<size; i++){
				myx += u8_get(u8);
				//printf("\t[SP_U8] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx);
			}
			break;

			case SPM_CSRDU_FL_UNIT_SP_U16:
			//printf("\t[SP_US] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
			ALIGN_PTR(u8, 2);
			yr += (*v++)*(*myx);
			for (i=1; i<size; i++){
				myx += u16_get(u8);
				//printf("\t[SP_US] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx);
			}
			break;

			case SPM_CSRDU_FL_UNIT_SP_U32:
			ALIGN_PTR(u8, 4);
			//printf("\t[SP_UI] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
			yr += (*v++) * (*myx);
			for (i=1; i<size; i++){
				myx += u32_get(u8);
				//printf("\t[SP_UI] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx);
			}
			break;

			case SPM_CSRDU_FL_UNIT_SP_U64:
			ALIGN_PTR(u8, 8);
			//printf("\t[SP_UL] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
			yr += (*v++) * (*myx);
			for (i=1; i<size; i++){
				myx += u32_get(u8);
				//printf("\t[SP_UL] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx);
			}
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

XSPMV_METH_INIT(
 _NAME(_jmp_multiply),
 SPM_CSRDU_NAME(_init_mmf),
 SPM_CSRDU_NAME(_size),
 NULL,
 sizeof(ELEM_TYPE)
)

