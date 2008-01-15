#include <stdlib.h>
#include <stdio.h>

#include "spm_mt.h"
#include "spm_crs_mt.h"
#include "spmv_loops_mt.h"

int main(int argc, char **argv)
{
	spm_mt_t *spm_mt;
	unsigned long rows_nr, cols_nr, nz_nr;
	
	if (argc<2){
		fprintf(stderr, "Usage: %s mmf_file\n", argv[0]);
		exit(1);
	}

	spm_mt = spm_crs32_double_mt_init_mmf(argv[1], &rows_nr, &cols_nr, &nz_nr);
	spm_crs32_double_mt_t *crs_mt;
	crs_mt = (spm_crs32_double_mt_t *)spm_mt->spm_threads->spm;

	spmv_double_check_mt_loop(crs_mt->crs,
	                          spm_mt, 
	                          spm_crs32_double_multiply, 
				  spm_crs32_double_mt_multiply,
				  1,
				  cols_nr);

	return 0;
}
