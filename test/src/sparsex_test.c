/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/*
 * \file sparsex_test.c
 * \brief Simple program for testing the SparseX API
 *
 * \author Athena Elafrou
 * \date 2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/sparsex.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

static const char *program_name;

static struct option long_options[] = {
    {"option",  required_argument,  0, 'o'},
    {"help",    no_argument,        0, 'h'}
};

static void print_usage()
{
    fprintf(stderr,
            "Usage: %s [-o <option=value>]... <mmf_file>\n\n"
            "\t-o, --option <option=value>\tset a preprocessing option\n"
            "\t-h, --help\t\t\tprint this help message and exit\n",
            basename(program_name));
}

static void set_option(const char *arg)
{
	// copy arg first to work with strtok
	char *temp = malloc(strlen(arg) + 1);
	temp = strncpy(temp, arg, strlen(arg) + 1);
	char *option = strtok(temp, "=");
	char *value = strtok(NULL, "=");
    spx_option_set(option, value);
	free(temp);
}

int main(int argc, char **argv)
{
    spx_value_t alpha = 0.8, beta = 0.42;
    const size_t loops = 128;

    spx_init();
    spx_log_info_console();

    char c;
    char *filename = NULL;
    char *option;
    int option_index = 0;

    program_name = argv[0];
    while ((c = getopt_long(argc, argv, "o:h", long_options,
                            &option_index)) != -1) {
        switch (c) {
        case 'o':
            option = optarg;
            set_option(option);
            break;
        case 'h':
            print_usage();
            exit(0);
        default:
            print_usage();
            exit(1);
        }
    }
    
	int remargs = argc - optind; // remaining arguments
	if (remargs < 1) {
		print_usage();
		exit(1);
	}
	argv = &argv[optind];
	filename = argv[0];

    /* Load matrix from MMF file */
    spx_input_t *input = spx_input_load_mmf(filename);

    /* Transform to CSX */
    spx_matrix_t *A = spx_mat_tune(input);

    /* Create random x and y vectors */
    spx_partition_t *parts = spx_mat_get_partition(A);
    spx_vector_t *x = spx_vec_create_random(spx_mat_get_ncols(A), parts);
    spx_vector_t *y = spx_vec_create_random(spx_mat_get_nrows(A), parts);

    /* Run 128 loops of the SpMV kernel */
    spx_timer_t t;
    double elapsed_time, flops;
    int i;

    spx_timer_clear(&t);
    spx_timer_start(&t);
    for (i = 0; i < loops; i++) {
        spx_matvec_kernel(alpha, A, x, beta, y);
    }
    spx_timer_pause(&t);
    elapsed_time = spx_timer_get_secs(&t);
    flops = (double) (2*loops*spx_mat_get_nnz(A)) /
        ((double) 1000*1000*elapsed_time);

    printf("Elapsed time: %lf secs\n", elapsed_time);
    printf("FLOPS: %lf\n", flops);

    /* Cleanup */
    spx_input_destroy(input);
    spx_mat_destroy(A);
    spx_partition_destroy(parts);
    spx_vec_destroy(x);
    spx_vec_destroy(y);

    return 0;
}
