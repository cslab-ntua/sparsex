/*
 * spmv.c -- control program
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h> /* PRIu64 */
#include <libgen.h>
#include <assert.h>

#include "method.h"
#include "spmv_method.h"
#include "spmv_loops.h"
#include "spmv_loops_mt.h"

static char *progname = NULL;
extern int optind;

static void help()
{
	fprintf(stderr, "Usage: %s [ -h -c -b [-l nr_loops] [-L outer_loops]] mmf_file [method]\n", progname);
	method_fprint(stderr, "available methods:\n", "\t", "\n", "");
}

int main(int argc, char **argv)
{

	progname = argv[0];

	// parse options
	int opt_check = 0;
	int opt_bench = 0;
	int loops_nr = 128;
	int outer_loops = 4;
	int c;
	while ((c = getopt(argc, argv, "hcbl:L:")) != -1){
		switch (c) {
			case 'c':
			opt_check = 1;
			break;

			case 'b':
			opt_bench = 1;
			break;

			case 'h':
			method_fprint(stdout, "", " ", "\n", "");
			exit(0);

			case 'L':
			outer_loops = atol(optarg);
			break;

			case 'l':
			loops_nr = atol(optarg);
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
	int elem_size = spmv_meth->elem_size;
	assert(elem_size == 8 || elem_size == 4);

	if (opt_check){
		method_t *meth1;
		spmv_method_t *spmv_meth1;
		uint64_t nrows1, ncols1, nnz1;
		void *m1;
		char *meth1_str = (elem_size == 8) ? "spm_crs32_double_multiply" : "spm_crs32_float_multiply";
		meth1 = method_get(meth1_str);
		spmv_meth1 = (spmv_method_t *)meth1->data;
		m1 = spmv_meth1->mmf_init_fn(mmf_file, &nrows1, &ncols1, &nnz1);
		if (nrows != nrows1 || ncols != ncols1 || nnz != nnz1){
			fprintf(stderr, "sizes do not match (f:%s) (m:%s)\n", mmf_file, method);
			exit(1);
		}

		switch (elem_size + spmv_meth->mt_flag){
			case 8:
			spmv_double_check_loop(m1, m, meth1->fn, meth->fn, 1, nrows, ncols, nnz);
			break;

			case (8+1):
			spmv_double_check_mt_loop(m1, m, meth1->fn, 1, nrows, ncols, meth->fn);
			break;

			case 4:
			spmv_float_check_loop(m1, m, meth1->fn, meth->fn, 1, nrows, ncols, nnz);
			break;

			case (4+1):
			spmv_float_check_mt_loop(m1, m, meth1->fn, 1, nrows, ncols, meth->fn);
			break;

			default:
			fprintf(stderr, "woops!\n");
			exit(1);
		}
		printf("check for %s in %s verified\n", method, basename(mmf_file));
		spmv_meth1->destroy_fn(m1);
	}

	if (opt_bench){
		int count;
		double t = -666.0;
		for (count=0; count < outer_loops; count++){
			switch (elem_size + spmv_meth->mt_flag){
				case 8:
				t = spmv_double_bench_loop(meth->fn, m, loops_nr, nrows, ncols);
				break;

				case (8+1):
				t = spmv_double_bench_mt_loop(m, loops_nr, nrows, ncols, meth->fn);
				break;

				case 4:
				t = spmv_float_bench_loop(meth->fn, m, loops_nr, nrows, ncols);
				break;

				case (4+1):
				t = spmv_float_bench_mt_loop(m, loops_nr, nrows, ncols, meth->fn);
				break;

				default:
				fprintf(stderr, "woops!\n");
				exit(1);
			}
			double flops = (double)(loops_nr*nnz*2)/((double)1000*1000*t);
			printf("m:%s f:%s s:%" PRIu64 " t:%lf r:%lf\n", method, basename(mmf_file), spmv_meth->size_fn(m), t, flops);
		}
	}
	spmv_meth->destroy_fn(m);

	return 0;
}
