#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <libgen.h>

#include "spm_crs.h"
#include "method.h"
#include "spm_crs_vh.h"
#include "spmv_loops.h"

#define LOOPS 128

int main(int argc, char **argv)
{
	char method_str[1024];
	int ci=32;
	char type[] = "double";

	if (argc < 2){
		fprintf(stderr, "Usage: %s mmf_file \n", argv[0]);
		exit(1);
	}

	snprintf(method_str, 1024, "spm_crs%d_vh_%s_multiply", ci, type);
	method_t *m = method_get(method_str);
	printf("%s\n", method_str);
	method_print("--\n", " ", "\n", "--\n");
	spmv_method_t *spmv_m = m->data;
	spmv_load_fn_t *mmf_init = spmv_m->mmf_init;
	void *crs;
	unsigned long rows_nr, cols_nr, nz_nr, loops;

	loops = LOOPS;

	void *crsvh = mmf_init(argv[1], &rows_nr, &cols_nr, &nz_nr);
	if ( spmv_m->elem_size != sizeof(double) ){
		fprintf(stderr, "sanity check failed (?)\n");
		exit(1);
	}

	crs = spm_crs32_double_init_mmf(argv[1], &rows_nr, &cols_nr, &nz_nr);
	spmv_double_check_loop(crs, crsvh, spm_crs32_double_multiply, m->fn, 1, cols_nr, nz_nr);
	return 0;
}
