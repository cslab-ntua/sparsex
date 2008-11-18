#include <stdlib.h>
#include <stdio.h>

#include "mmf.h"
#include "method.h"
#include "spm_crs.h"
#include "spmv_loops.h"
#include "jit/crsvh.h"

#define LOOPS 128

int main(int argc, char **argv)
{
	if (argc < 2){
		fprintf(stderr, "Usage: %s mmf_file \n", argv[0]);
		exit(1);
	}

	void *crs;
	unsigned long rows_nr, cols_nr, nz_nr;
	method_t *m = spm_crs32_vhjit_double_init_method(argv[1], &rows_nr, &cols_nr, &nz_nr);
	crs = spm_crs32_double_init_mmf(argv[1], &rows_nr, &cols_nr, &nz_nr);
	spmv_double_check_loop(crs, m->data, spm_crs32_double_multiply, m->fn, 1, rows_nr, cols_nr, nz_nr);
	return 0;
}
