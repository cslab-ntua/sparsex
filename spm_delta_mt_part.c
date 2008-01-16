#include <stdlib.h>
#include <stdlib.h>
#include <assert.h>

#include "mt_lib.h"
#include "spm_mt.h"
#include "spm_delta.h"
#include "spm_delta_mt.h"

typedef struct {
	unsigned long nnz, ncols, nrows;
	unsigned char *ctl;
} spm_delta_idx_t;


spm_mt_t *spm_delta_jmp_mt_partition(unsigned char *ctl, unsigned long nnz, void *delta)
{
	unsigned char *uc = ctl;
	unsigned int nr_cpus, *cpus_affinity;
	spm_mt_t *spm_mt;
	spm_delta_mt_t *delta_mt;

	mt_get_options(&nr_cpus, &cpus_affinity);

	spm_mt = malloc(sizeof(spm_mt_t));
	if ( !spm_mt ){
		perror("malloc");
		exit(1);
	}

	spm_mt->nr_threads = nr_cpus;
	spm_mt->spm_threads = malloc(sizeof(spm_mt_thread_t)*nr_cpus);
	if ( !spm_mt->spm_threads ){
		perror("malloc");
		exit(1);
	}

	delta_mt = malloc(sizeof(spm_delta_mt_t)*nr_cpus);
	if ( !delta_mt ){
		perror("malloc");
		exit(1);
	}

	unsigned long elements_total=0, elements=0;
	unsigned long elements_limit = nnz / nr_cpus;
	unsigned long ctl_last_off=0;
	int delta_mt_idx=0;
	unsigned long y_indx = 0;
	for (;;) {
		register unsigned char flags = uc_get(uc);
		register unsigned char size = uc_get(uc);
		elements += size;
		//spm_delta_ctl_print(flags, size);
		if ( spm_delta_fl_isnr(flags) ){
			if ( spm_delta_fl_isnrseq(flags) ){
				y_indx++;
			} else {
				y_indx += uc_get_ul(uc);
			}
		}

		#define ALIGN(buf,a) (void *) (((unsigned long) (buf) + (a-1)) & ~(a-1))
		#define ALIGN_UC(align) (uc = ALIGN(uc, align))

		uc_get_ul(uc);
		switch ( flags & SPM_DELTA_FL_UNIT_MASK ){

			case SPM_DELTA_FL_UNIT_SP_UC:
			uc += sizeof(unsigned char)*(size-1);
			break;

			case SPM_DELTA_FL_UNIT_SP_US:
			ALIGN_UC(2);
			uc +=  sizeof(unsigned short)*(size-1);
			break;

			case SPM_DELTA_FL_UNIT_SP_UI:
			ALIGN_UC(4);
			uc +=  sizeof(unsigned int)*(size-1);
			break;

			case SPM_DELTA_FL_UNIT_SP_UL:
			ALIGN_UC(8);
			uc += sizeof(unsigned long)*(size-1);
			break;

			case SPM_DELTA_FL_UNIT_DENSE:
			break;

			default:
			printf("Uknown flags: %u unit=%u\n", flags, flags & SPM_DELTA_FL_UNIT_MASK);
			assert(0);
		}

		//printf("elements:%lu elements_total:%lu nnz:%lu\n", elements, elements_total, nnz);
		if (elements >= elements_limit) {
			spm_mt_thread_t *spm_mt_thread = spm_mt->spm_threads + delta_mt_idx;
			spm_delta_mt_t *dmt = delta_mt + delta_mt_idx;

			spm_mt_thread->cpu = cpus_affinity[delta_mt_idx];
			spm_mt_thread->spm = dmt;

			dmt->delta = delta;
			dmt->nnz = elements;
			dmt->ctl_start = ctl_last_off;

			ctl_last_off = dmt->ctl_end = (uc - ctl);
			elements_total += elements;
			elements_limit = (nnz - elements_total) / (nr_cpus - delta_mt_idx);
			elements = 0;

			if ( ++delta_mt_idx == nr_cpus){
				break;
			}
		}

		assert(elements_total < nnz);

	}

	free(cpus_affinity);
	return spm_mt;
}
