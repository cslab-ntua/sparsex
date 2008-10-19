#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <libgen.h>

#include "method.h"
#include "spm_mt.h"
#include "spmv_loops_mt.h"
#include "spm_crs.h"

#define LOOPS 128

int main(int argc, char **argv)
{
	char method_str[1024];
	int ci=32;
	char type[] = "double";

	if (argc < 2){
		fprintf(stderr, "Usage: %s mmf_file\n", argv[0]);
		exit(1);
	}

	snprintf(method_str, 1024, "spm_crs%d_vh_%s_mt_multiply", ci, type);
	method_t *m = method_get(method_str);
	if (!m){
		fprintf(stderr, "method: %s does not exist\n", method_str);
		method_fprint(stderr, "available methods: \n", "  ", "\n", "--\n");
		exit(1);
	}

	spmv_method_t *spmv_m = m->data;
	spmv_load_fn_t *mmf_init = spmv_m->mmf_init;
	spmv_size_fn_t *spm_size = spmv_m->size;
	spm_mt_t *spm_mt;
	unsigned long rows_nr, cols_nr, nz_nr, loops;

	loops = LOOPS;

	void *crs = spm_crs32_double_init_mmf(argv[1], &rows_nr, &cols_nr, &nz_nr);

	spm_mt = mmf_init(argv[1], &rows_nr, &cols_nr, &nz_nr);
	if ( spmv_m->elem_size != sizeof(double) ){
		fprintf(stderr, "sanity check failed (?)\n");
		exit(1);
	}

	//double time = spmv_double_bench_mt_loop(m->fn, spm_mt, loops, cols_nr);
	//double flops = (double)(loops*nz_nr*2)/(1000*1000*time);
	//printf("%s %s %lu %lf %lf\n", basename(argv[1]), method_str, spm_size(spm_mt), time, flops);

	spmv_double_check_mt_loop(crs, spm_mt,
	                          spm_crs32_double_multiply, m->fn,
				  1, cols_nr);


	return 0;
}
