#include <stdlib.h>
#include <assert.h>
#include <pthread.h>

#include "spm_delta.h"
#include "spm_delta_mt.h"
#include "vector.h"

#define SPMDELTA_CON3_(a,b,c) a ## b ## c
#define SPMDELTA_CON3(a,b,c) SPMDELTA_CON3_(a,b,c)
#define SPMDELTA_NAME(name) SPMDELTA_CON3(spm_delta_, ELEM_TYPE, name)
#define SPMDELTA_MT_NAME(name) SPMDELTA_CON3(spm_delta_mt_, ELEM_TYPE, name)
#define SPMDELTA_TYPE SPMDELTA_NAME(_t)

#if 0
static void spm_delta_ctl_print(unsigned char flags, unsigned int size)
{
	printf("nr:%d sp:%d nr_seq:%d ci_size:%d size:%d\n",
               !!spm_delta_fl_isnr(flags),
               !!spm_delta_fl_issp(flags),
               !!spm_delta_fl_isnrseq(flags),
	       spm_delta_fl_getcisize(flags),
	       size);
}
#endif


#if 0
/* XXX: This is br0ken */
void SPMDELTA_MT_NAME(_multiply)(void *matrix, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	spm_delta_mt_t *delta_mt = matrix;
	SPMDELTA_TYPE *delta = (SPMDELTA_TYPE *)delta_mt->delta;
	ELEM_TYPE *x = in->elements, *y = out->elements, *values = delta->values;
	register unsigned long y_indx=0;
	register ELEM_TYPE yr = 0;
	register ELEM_TYPE *v = values;
	register const ELEM_TYPE *v_end = values + delta_mt->nnz;
	register ELEM_TYPE *myx = x;
	register unsigned char *uc = delta->ctl + delta_mt->ctl_start;

	for (;;) {
		register unsigned char flags = uc_get(uc);
		register unsigned char size = uc_get(uc);
		//spm_delta_ctl_print(flags, size);
		if ( spm_delta_fl_isnr(flags) ){
			y[y_indx] = yr;
			if ( spm_delta_fl_isnrseq(flags) ){
				y_indx++;
			} else {
				y_indx += uc_get_ul(uc);
			}
			//printf("--new_row: (result for %lu: %lf)\n", y_indx, yr);
			yr = 0;
			myx = x;
		}
		//printf("y_indx: %lu x_indx: %lu\n", y_indx, (unsigned long)(myx-x));

		switch ( flags & SPM_DELTA_FL_UNIT_MASK ){
			register unsigned long i;

			case SPM_DELTA_FL_UNIT_DENSE:
			__asm__ __volatile__ (" # DE loop ");
			myx += uc_get_ul(uc);
			for ( i = 0; i<size; i++){
				yr += (*v) * (*myx);
				//printf("\t[DE] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				v++;
				myx++;
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_UC:
			for ( i = 0; i<size; i++){
				//printf("\t[SP_UC] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				myx += uc_get(uc);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_US:
			for ( i = 0; i<size; i++){
				myx += us_get(uc);
				//printf("\t[SP_US] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_UI:
			for ( i = 0; i<size; i++){
				myx += ui_get(uc);
				//printf("\t[SP_UI] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_UL:
			for ( i = 0; i<size; i++){
				myx += ul_get(uc);
				//printf("\t[SP_UL] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			default:
			printf("Uknown flags: %u unit=%u\n", flags, flags & SPM_DELTA_FL_UNIT_MASK);
			assert(0);
		}

		if ( v >= v_end ){
			break;
		}
	}

	y[y_indx] = yr;
	return;
}
//METHOD_INIT(spm_delta_multiply, spm_delta_init_mmf)
#endif

void SPMDELTA_MT_NAME(_jmp_multiply)(void *matrix, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	spm_delta_mt_t *delta_mt = matrix;
	SPMDELTA_TYPE *delta = (SPMDELTA_TYPE *)delta_mt->delta;
	ELEM_TYPE *x = in->elements, *y = out->elements, *values = delta->values;
	register unsigned long y_indx = delta_mt->row_start;
	register ELEM_TYPE yr = (ELEM_TYPE)0;
	register ELEM_TYPE *v = values + delta_mt->val_start;
	register const ELEM_TYPE *v_end = v + delta_mt->nnz;
	register ELEM_TYPE *myx = x;
	register unsigned char *uc = delta->ctl + delta_mt->ctl_start;

	//printf("y_indx:%lu ctl:%p ctl_start:%lu (%lu)\n", y_indx, delta->ctl, delta_mt->ctl_start, (unsigned long)pthread_self());
	for (;;) {
		register unsigned char flags = uc_get(uc);
		register unsigned char size = uc_get(uc);
		//spm_delta_ctl_print(flags, size);
		if ( spm_delta_fl_isnr(flags) ){
			y[y_indx] = yr;
			//printf("--new_row (%lu): (result for %lu: %lf)\n", (unsigned long)pthread_self(), y_indx, yr);
			if ( spm_delta_fl_isnrseq(flags) ){
				y_indx++;
			} else {
				y_indx += uc_get_ul(uc);
			}
			yr = (ELEM_TYPE)0;
			myx = x;
		}
		//printf("y_indx: %lu x_indx: %lu\n", y_indx, (unsigned long)(myx-x));

		myx += uc_get_ul(uc);
		switch ( flags & SPM_DELTA_FL_UNIT_MASK ){
			register unsigned long i;

			case SPM_DELTA_FL_UNIT_DENSE:
			__asm__ __volatile__ (" # DE loop ");
			for (i=0; i<size; i++){
				//printf("\t[DE] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx++);
			}
			break;

			#define ALIGN(buf,a) (void *) (((unsigned long) (buf) + (a-1)) & ~(a-1))
			#define ALIGN_UC(align) (uc = ALIGN(uc, align))

			case SPM_DELTA_FL_UNIT_SP_UC:
			//printf("\t[SP_UC] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
			yr += (*v++)*(*myx);
			for (i=1; i<size; i++){
				myx += uc_get(uc);
				//printf("\t[SP_UC] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx);
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_US:
			//printf("\t[SP_US] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
			ALIGN_UC(2);
			yr += (*v++)*(*myx);
			for (i=1; i<size; i++){
				myx += us_get(uc);
				//printf("\t[SP_US] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx);
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_UI:
			ALIGN_UC(4);
			//printf("\t[SP_UI] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
			yr += (*v++) * (*myx);
			for (i=1; i<size; i++){
				myx += ui_get(uc);
				//printf("\t[SP_UI] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx);
			}
			break;

			#if 0
			case SPM_DELTA_FL_UNIT_SP_UL:
			printf("\t[SP_UL] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
			yr += (*v++) * (*myx);
			for (i=1; i<size; i++){
				myx += ul_get(uc);
				printf("\t[SP_UL] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				yr += (*v++) * (*myx);
			}
			break;
			#endif

			default:
			printf("Uknown flags: %u unit=%u\n", flags, flags & SPM_DELTA_FL_UNIT_MASK);
			assert(0);
		}

		if ( v >= v_end ){
			break;
		}
	}

	y[y_indx] = yr;
	//printf("--end (%lu): (result for %lu: %lf)\n", (unsigned long)pthread_self(), y_indx, yr);
	//printf("(%lu) : last uc: %lu\n", (unsigned long)pthread_self(), uc - delta->ctl);
	return;
}

