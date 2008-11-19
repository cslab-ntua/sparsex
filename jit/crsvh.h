#ifndef SPM_CRS_VHJIT__
#define SPM_CRS_VHJIT__

#include <inttypes.h>

#include "method.h"
#include "spm_mt.h"

typedef struct {
	uint32_t *col_ind, *row_ptr;
	uint64_t nz, nrows, ncols, row_start, row_end;
	unsigned char *hvals;
} spm_crs32_vhjit_t;

method_t *spm_crs32_vhjit_double_init_method(char *mmf_file, unsigned long *rows_nr,
                                             unsigned long *cols_nr, unsigned long *nz_nr);


spm_mt_t *spm_crs32_vhjit_mt_double_init_method(char *mmf_file, unsigned long *rows_nr,
                                                unsigned long *cols_nr, unsigned long *nz_nr);

#endif /* SPM_CRS_VHJIT__ */
