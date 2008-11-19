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

	spm_mt = spm_crs32_vhjit_mt_double_init_method(argv[1], &rows_nr, &cols_nr, &nz_nr);
	double time = spmv_double_bench_mt_loop(spm_mt, loops, rows_nr, cols_nr, NULL);
	double flops = (double)(loops*nz_nr*2)/(1000*1000*time);
	printf("%s %s %lu %lf %lf\n", basename(argv[1]), "spmv_vhjit_mt", 0L, time, flops);

	return 0;
}
