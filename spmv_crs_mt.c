#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>

#include "spm_mt.h"
#include "spm_crs_mt.h"
#include "spmv_loops_mt.h"

#define LOOPS 128

int main(int argc, char **argv)
{
	spm_mt_t *spm_mt;
	unsigned long rows_nr, cols_nr, nz_nr, loops, size;
	double time, flops;

	if (argc < 2){
		fprintf(stderr, "Usage: %s mmf_file\n", argv[0]);
		exit(1);
	}

	loops = LOOPS;

	spm_mt = spm_crs32_double_mt_init_mmf(argv[1], &rows_nr, &cols_nr, &nz_nr);
	time = spmv_double_bench_mt_loop(spm_crs32_double_mt_multiply, spm_mt,loops, cols_nr);
	flops = (double)(loops*nz_nr*2)/(1000*1000*time);
	size = nz_nr*(8 + 4) + cols_nr*4;

	printf("%s spm_crs32_double_mt %lu %lf %lf\n", basename(argv[1]), size, time, flops);

	return 0;
}
