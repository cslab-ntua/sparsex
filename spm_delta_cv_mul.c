#include <stdlib.h>
#include <assert.h>

#include "spm_delta_cv.h"
#include "vector.h"

#define SPMDELTACV_CON3_(a,b,c) a ## b ## c
#define SPMDELTACV_CON3(a,b,c) SPMDELTACV_CON3_(a,b,c)
#define SPMDELTACV_NAME(name) SPMDELTACV_CON3(spm_delta_cv_, ELEM_TYPE, name)
#define SPMDELTACV_TYPE SPMDELTACV_NAME(_t)

#define DBG_PRINT
#undef DBG_PRINT

#define asm_comment(comment) __asm__ __volatile__ ( "# " comment  "\n")

static inline void spm_delta_cv_ctl_print(unsigned char flags, unsigned int size)
{
	#ifdef DBG_PRINT
	printf("nr:%d sp:%d nr_seq:%d ci_size:%d size:%d\n", 
               !!spm_delta_fl_isnr(flags), 
               !!spm_delta_fl_issp(flags), 
               !!spm_delta_fl_isnrseq(flags),
	       spm_delta_fl_getcisize(flags),
	       size);
	 #endif
	 return;
}

void SPMDELTACV_NAME(_jmp_multiply)(void *matrix, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPMDELTACV_TYPE *d = (SPMDELTACV_TYPE *)matrix;
	ELEM_TYPE *x = in->elements, *y = out->elements, *v = d->values;
	register unsigned long y_indx=0;
	register ELEM_TYPE yr = 0;
	register unsigned long v_idx, x_idx=0;
	register unsigned char *ctl = d->ctl, *ci = d->c_indices, *vi=d->v_indices;
	register unsigned char *ctl_end = d->ctl + d->ctl_size;

	for (;;) {
		register unsigned char flags = u8_get(ctl);
		register unsigned char size = u8_get(ctl);
		spm_delta_cv_ctl_print(flags, size);
		if ( spm_delta_fl_isnr(flags) ){
			y[y_indx] = yr;
			#ifdef DBG_PRINT
			printf("--new_row: (result for %lu: %lf)\n", y_indx, yr);
			#endif
			if ( 1 || spm_delta_fl_isnrseq(flags) ){
				y_indx++;
			} else {
				y_indx += uc_get_ul(ctl);
			}
			yr = 0;
			x_idx = 0;
		} 
		#ifdef DBG_PRINT
		printf("y_indx: %lu x_indx: %lu\n", y_indx, (unsigned long)(myx-x));
		#endif
		x_idx += uc_get_ul(ctl);
		register unsigned long i=0;
		switch ( flags & SPM_DELTA_FL_CV_UMASK ){
			/* In bp and cpp we trust .. */

			#ifdef DBG_PRINT
			#define dbg_print_sp(__ci, __vi) \
				printf("[SP ci" # __ci " vi" # __vi "]: yo=%lu x=%lf xo=%lu v=%lf vo=%ld \n", \
				y_indx, *myx, myx -x, *myv, myv -v )
			#else
			#define dbg_print_sp(__ci, __vi) do { } while(0)
			#endif

			#define CASE_SP(__ci, __vi) \
			case SPM_DELTA_CI ## __ci ## _VI ## __vi :         \
			asm_comment("Start ci=" # __ci " vi=" # __vi);     \
			for (;;){                                          \
				v_idx = u ## __vi ## _get(vi);             \
				dbg_print_sp(__ci, __vi);                  \
				yr += v[v_idx] * x[x_idx];                 \
				if ( ++i <= size ){                        \
					x_idx += u ## __ci ## _get(ci);    \
				} else {                                   \
					break;                             \
				}                                          \
			}                                                  \
			asm_comment("End ci=" # __ci " vi=" # __vi);       \
			break;

			CASE_SP(8,8)
			CASE_SP(16,8)
			CASE_SP(32,8)

			CASE_SP(8,16)
			CASE_SP(16,16)
			CASE_SP(32,16)

			CASE_SP(8,32)
			CASE_SP(16,32)
			CASE_SP(32,32)

			#ifdef DBG_PRINT
			#define dbg_print_de(__vi) \
				printf("[DE vi" # __vi "]: yo=%lu x=%lf xo=%lu v=%lf vo=%ld \n", \
				y_indx, myx[i], myx+i-x, *myv, myv -v )
			#else
			#define dbg_print_de(__vi) do { } while(0)
			#endif

			#define CASE_DE(__vi) \
			case SPM_DELTA_DE_VI ## __vi :          \
			for (;i<size;i++){                      \
				v_idx = u ## __vi ## _get(vi); \
				dbg_print_de(__vi);             \
				yr += v[v_idx] * x[x_idx];      \
			}                                       \
			x_idx += size;                           \
			break;

			CASE_DE(8)
			CASE_DE(16)
			CASE_DE(32)

			default:
			printf("Uknown flags: %u unit=%u\n", flags, flags & SPM_DELTA_FL_UNIT_MASK);
			assert(0);
		}

		if ( ctl == ctl_end ){
			break;
		}
	}
	y[y_indx] = yr;
}
