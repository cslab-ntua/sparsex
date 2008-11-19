#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <libgen.h>

#include "spm_mt.h"
#include "spm_crs.h"
#include "spmv_loops_mt.h"
#include "jit/crsvh.h"

#define LOOPS 128

int main(int argc, char **argv)
{
	if (argc < 2){
		fprintf(stderr, "Usage: %s mmf_file\n", argv[0]);
		exit(1);
	}

	spm_mt_t *spm_mt;
	unsigned long rows_nr, cols_nr, nz_nr, loops;

	loops = LOOPS;

	void *crs;
	crs = spm_crs32_double_init_mmf(argv[1], &rows_nr, &cols_nr, &nz_nr);
	spm_mt = spm_crs32_vhjit_mt_double_init_method(argv[1], &rows_nr, &cols_nr, &nz_nr);
	spmv_double_check_mt_loop(crs, spm_mt,
	                          spm_crs32_double_multiply, 1,
				  rows_nr, cols_nr, NULL);
	return 0;
}
