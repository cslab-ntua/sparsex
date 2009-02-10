#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>
#include <libgen.h>

#include "method.h"
#include "spmv_method.h"
#include "spmv_loops.h"

static char *progname = NULL;
extern int optind;

static void help()
{
	fprintf(stderr, "Usage: %s [-c -b -l [nr_loops]] mmf_file [method]\n", progname);
	method_fprint(stderr, "available methods:\n", "\t", "\n", "");
}

int main(int argc, char **argv)
{

	progname = argv[0];

	// parse options
	int opt_check = 0;
	int opt_bench = 0;
	int loops_nr = 128;
	int c;
	while ((c = getopt(argc, argv, "cbl:")) != -1){
		switch (c) {
			case 'c':
			opt_check = 1;
			break;

			case 'b':
			opt_bench = 1;
			break;

			default:
			fprintf(stderr, "Error parsing arguments: -%c-\n", c);
			help();
			exit(1);
		}
	}

	if (!opt_check && !opt_bench){
		opt_bench = 1;
	}

	int remargs=argc - optind; // remaining arguments
	if (remargs < 1){
		help();
		exit(1);
	}
	argv = &argv[optind-1];

	char *mmf_file = argv[1];
	char *method = (remargs > 1) ? argv[2] : "spm_crs32_double_multiply";
	method_t *meth = method_get(method);
	if (!meth) {
		fprintf(stderr, "No such method: %s\n", method);
		method_fprint(stderr, "available methods:\n", "\t", "\n", "");
		exit(1);
	}

	spmv_method_t *spmv_meth = meth->data;
	//spm_size_fn_t *size_fn = spmv_meth->size_fn;

	uint64_t nrows, ncols, nnz;
	void *m = spmv_meth->mmf_init_fn(mmf_file, &nrows, &ncols, &nnz);

	if (opt_check){
		switch (spmv_meth->elem_size) {
			method_t *meth1;
			spmv_method_t *spmv_meth1;
			uint64_t nrows1, ncols1, nnz1;
			void *m1;

			case 8:
			meth1 = method_get("spm_crs32_double_multiply");
			spmv_meth1 = (spmv_method_t *)meth1->data;
			m1 = spmv_meth1->mmf_init_fn(mmf_file, &nrows1, &ncols1, &nnz1);
			if (nrows != nrows1 || ncols != ncols1 || nnz != nnz1){
				fprintf(stderr, "sizes do not match (f:%s) (m:%s)\n", mmf_file, method);
				exit(1);
			}
			spmv_double_check_loop(m1, m, meth1->fn, meth->fn, 1, nrows, ncols, nnz);
			break;

			case 4:
			meth1 = method_get("spm_crs32_float_multiply");
			spmv_meth1 = (spmv_method_t *)meth1->data;
			m1 = spmv_meth1->mmf_init_fn(mmf_file, &nrows1, &ncols1, &nnz1);
			if (nrows != nrows1 || ncols != ncols1 || nnz != nnz1){
				fprintf(stderr, "sizes do not match (f:%s) (m:%s)\n", mmf_file, method);
				exit(1);
			}
			spmv_float_check_loop(m1, m, meth1->fn, meth->fn, 1, nrows, ncols, nnz);
			break;

			default:
			fprintf(stderr, "bogus elem_size (%d)\n", spmv_meth->elem_size);
			exit(1);
		}
		printf("check for %s in %s verified\n", method, basename(mmf_file));
	}

	if (opt_bench){
		double t;
		switch (spmv_meth->elem_size){
			case 8:
			t = spmv_double_bench_loop(meth->fn, m, loops_nr, nrows, ncols);
			break;

			case 4:
			t = spmv_float_bench_loop(meth->fn, m, loops_nr, nrows, ncols);
			break;

			default:
			fprintf(stderr, "bogus elem_size (%d)\n", spmv_meth->elem_size);
			exit(1);
		}
		printf("%s %s %lf\n", method, basename(mmf_file), t);
	}

	return 0;
}
