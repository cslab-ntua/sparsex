/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file csr_example.c
 * \brief Example 2
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/sparsex.h>

int main(int argc, char **argv)
{
    /* Define CSR data structures */
    spx_index_t nrows = 10, ncols = 10;
    spx_index_t rowptr[] = {0,5,6,10,15,18,22,24,29,33,38};
    spx_index_t colind[] = {0,1,2,3,8,7,0,1,6,9,0,1,3,5,9,0,1,9,0,1,5,9,2,3,
                            2,3,4,5,7,2,3,4,5,2,3,4,5,9};
    spx_value_t values[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                            20,21,22,23,24,25,26,26.1,26.2,27,28,29,29.1,29.2,
                            30,31,31.1,31.2,32};

    /* Define random x and y arrays */
    spx_value_t *x = (spx_value_t *) malloc(sizeof(spx_value_t) * ncols);
    spx_value_t *y = (spx_value_t *) malloc(sizeof(spx_value_t) * nrows);
    spx_value_t val = 0, max = 1, min = -1;
    size_t i;
    for (i = 0; i < nrows; i++) {
		val = ((spx_value_t) (rand()+i) / ((spx_value_t) RAND_MAX + 1));
		x[i] = min + val*(max-min);
		y[i] = max + val*(min-max);
    }

    /* Initialize library */
    spx_init();
    spx_log_info_console();

    /* Partition input matrix */
    spx_partition_t *parts = spx_partition_csr(rowptr, nrows, 2);

    /* Create x and y vector views */
    spx_vector_t *x_view = spx_vec_create_from_buff(x, ncols, parts,
                                                    SPX_VEC_SHARE);
    spx_vector_t *y_view = spx_vec_create_from_buff(y, nrows, parts,
                                                    SPX_VEC_SHARE);

    /* Declare a tuned matrix handle */
    spx_matrix_t *A = SPX_INVALID_MAT;

    /* Set tuning options */
    spx_option_set("spx.rt.nr_threads", "2");
    spx_option_set("spx.rt.cpu_affinity", "0,1");
    spx_option_set("spx.preproc.xform", "all");

    /* Transform to CSX and run the SpMV kernel */
    spx_value_t alpha = 0.8, beta = 0.42;
    const size_t nr_loops = 1500;
    spx_timer_t t;
    double elapsed_time, flops;

    spx_timer_clear(&t);
    spx_timer_start(&t);
    for (i = 0; i < nr_loops; i++) {
        spx_matvec_kernel_csr(&A, nrows, ncols, rowptr, colind, values,
                              alpha, x_view, beta, y_view);
    }
    spx_timer_pause(&t);

    elapsed_time = spx_timer_get_secs(&t);
    flops = (double) (2 * nr_loops * spx_mat_get_nnz(A)) /
        ((double) 1000 * 1000 * elapsed_time);
    printf("Elapsed time: %lf secs\n", elapsed_time);
    printf("FLOPS: %lf\n", flops);

    /* Cleanup */
    spx_mat_destroy(A);
    spx_partition_destroy(parts);
    spx_vec_destroy(x_view);
    spx_vec_destroy(y_view);

    /* Shutdown library */
    spx_finalize();

    return 0;
}
