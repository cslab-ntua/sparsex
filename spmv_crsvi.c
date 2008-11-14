#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <libgen.h>

#include "method.h"
#include "spm_crs_vi.h"
#include "spmv_loops.h"

#define LOOPS 128

int main(int argc, char **argv)
{
	char method_str[1024];
	int ci=32, vi=8;
	char type[] = "double";

	if (argc < 2){
		fprintf(stderr, "Usage: %s mmf_file <nr_distinct_values>\n", argv[0]);
		exit(1);
	}

	if (argc > 2){
		unsigned long distinct_vals = atol(argv[2]);
		if ( !distinct_vals ){
			fprintf(stderr, "invalid nr_distinct_values: %s\n", argv[2]);
			exit(1);
		}
		if ( distinct_vals < 1UL<<(sizeof(uint8_t)*8) ){
			vi = 8;
		} else if ( distinct_vals < 1UL<<(sizeof(uint16_t)*8) ){
			vi = 16;
		} else if ( distinct_vals < 1UL<<(sizeof(uint32_t)*8) ){
			vi = 32;
		} else {
			fprintf(stderr, "distinct_vals: %lu too big\n", distinct_vals);
			exit(1);
		}
	}

	snprintf(method_str, 1024, "spm_crs%d_vi%d_%s_multiply", ci, vi, type);
	method_t *m = method_get(method_str);
	spmv_method_t *spmv_m = m->data;
	spmv_load_fn_t *mmf_init = spmv_m->mmf_init;
	spmv_size_fn_t *spm_size = spmv_m->size;
	void *crsvi;
	unsigned long rows_nr, cols_nr, nz_nr, loops;
	double time=0, flops;

	loops = LOOPS;

	crsvi = mmf_init(argv[1], &rows_nr, &cols_nr, &nz_nr);
	if ( spmv_m->elem_size != sizeof(double) ){
		fprintf(stderr, "sanity check failed (?)\n");
		exit(1);
	}
	time = spmv_double_bench_loop(m->fn, crsvi,loops, rows_nr, cols_nr);
	flops = (double)(loops*nz_nr*2)/(1000*1000*time);

	printf("%s %s %lu %lf %lf\n", basename(argv[1]), method_str, spm_size(crsvi), time, flops);
	
	return 0;
}
