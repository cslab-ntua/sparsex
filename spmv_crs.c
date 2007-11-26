#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>

#include "spm_crs.h"
#include "spmv_loops.h"

#define LOOPS 128

int main(int argc, char **argv)
{
	spm_crs32_double_t *crs;
	unsigned long rows_nr, cols_nr, nz_nr, loops, size;
	double time, flops;

	if (argc < 2){
		fprintf(stderr, "Usage: %s mmf_file\n", argv[0]);
		exit(1);
	}

	loops = LOOPS;

	crs = spm_crs32_double_init_mmf(argv[1], &rows_nr, &cols_nr, &nz_nr);
	time = spmv_double_bench_loop(spm_crs32_double_multiply, crs,loops, cols_nr);
	flops = (double)(loops*nz_nr*2)/(1000*1000*time);
	size = nz_nr*(8 + 4) + cols_nr*4;

	printf("spm_crs32_double %s %lu %lf %lf\n", basename(argv[1]), size, time, flops);
	
	return 0;
}
