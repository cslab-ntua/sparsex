/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file matrix_caching_example_p2.c
 * \brief Read CSX matrix from a binary file.
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/sparsex.h>
#include <stdio.h>

static char *program_name;

void print_usage()
{
    fprintf(stderr, "Usage: %s <binary_file>\n", program_name);
}

int main(int argc, char **argv)
{
    program_name = basename(argv[0]);

    /* Initialize library */
    spx_init();
    spx_log_info_console();

    if (argc < 2) {
        fprintf(stderr, "%s: too few arguments\n", program_name);
        print_usage();
        exit(1);
    }

    /* Load matrix from binary file */
    spx_matrix_t *A = spx_mat_restore(argv[1]);

    /* Change some values in the matrix */
    spx_index_t row = 1, col = 1;
    spx_value_t value = 0.42;
    spx_mat_set_entry(A, row, col, value);
    /* ... */

    /* Create random x and y vectors */
    spx_partition_t *parts = spx_mat_get_partition(A);
    spx_vector_t *x = spx_vec_create_random(spx_mat_get_ncols(A), parts);
    spx_vector_t *y = spx_vec_create_random(spx_mat_get_nrows(A), parts);

    /* Run 128 loops of the SpMV kernel */
    spx_value_t alpha = 0.8, beta = 0.42;
    const size_t nr_loops = 128;
    spx_timer_t t;
    double elapsed_time, flops;
    size_t i;

    spx_timer_clear(&t);
    spx_timer_start(&t);
    for (i = 0; i < nr_loops; i++) {
        spx_matvec_kernel(alpha, A, x, beta, y);
    }
    spx_timer_pause(&t);

    elapsed_time = spx_timer_get_secs(&t);
    flops = (double) (2 * nr_loops * spx_mat_get_nnz(A)) /
        ((double) 1000 * 1000 * elapsed_time);
    printf("SPMV time: %lf secs\n", elapsed_time);
    printf("MFLOPS: %lf\n", flops);

    /* Cleanup */
    spx_mat_destroy(A);
    spx_partition_destroy(parts);
    spx_vec_destroy(x);
    spx_vec_destroy(y);

    /* Shutdown library */
    spx_finalize();

    return 0;
}
