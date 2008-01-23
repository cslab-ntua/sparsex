#ifndef __SPM_DELTA_MT_H__
#define __SPM_DELTA_MT_H__

#include "spm_mt.h"
#include "spmv_method.h"

typedef struct {
	void *delta;
	unsigned long nnz, ctl_start, row_start, val_start;
} spm_delta_mt_t;

spm_mt_t *spm_delta_jmp_mt_partition(unsigned char *ctl, unsigned long nnz, void *delta);

spmv_float_fn_t spm_delta_mt_float_multiply;
spmv_float_fn_t spm_delta_mt_float_jmp_multiply;

spmv_double_fn_t spm_delta_mt_double_multiply;
spmv_double_fn_t spm_delta_mt_double_jmp_multiply;


#endif /* __SPM_DELTA_H__ */
