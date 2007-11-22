#include <stdlib.h>
#include <stdio.h>

#include "spm_crs_vi.h"
#include "spmv_loops.h"

#define LOOPS 128

int main(int argc, char **argv)
{
	spm_crsvi32_double_t *crsvi;
	unsigned long rows_nr, cols_nr, nz_nr, loops;
	double time, flops;

	if (argc < 2){
		fprintf(stderr, "Usage: %s mmf_file\n", argv[0]);
		exit(1);
	}

	loops = LOOPS;

	crsvi = spm_crsvi32_double_init_mmf(argv[1], &rows_nr, &cols_nr, &nz_nr);
	time = spmv_double_bench_loop(spm_crsvi32_double_multiply, crsvi,loops, cols_nr);
	flops = (double)(loops*nz_nr*2)/(1000*1000*time);

	printf("spm_crsvi32_double: %lf %lf\n", time, flops);
	
	return 0;
}
