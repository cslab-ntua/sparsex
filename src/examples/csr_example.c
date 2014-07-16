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
    /* Create CSR data structures */
    spx_index_t nrows = 10, ncols = 10;
    spx_index_t rowptr[] = {0,5,6,10,15,18,22,24,29,33,38};
    spx_index_t colind[] = {0,1,2,3,8,7,0,1,6,9,0,1,3,5,9,0,1,9,0,1,5,9,2,3,
                            2,3,4,5,7,2,3,4,5,2,3,4,5,9};
    spx_value_t values[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                            20,21,22,23,24,25,26,26.1,26.2,27,28,29,29.1,29.2,
                            30,31,31.1,31.2,32};

    /* Create random x and y vectors */
	spx_value_t *x = (spx_value_t *) malloc(sizeof(spx_value_t) * ncols);
	spx_value_t *y = (spx_value_t *) malloc(sizeof(spx_value_t) * nrows);
    spx_value_t val = 0, max = 1, min = -1;
    size_t i;
    for (i = 0; i < nrows; i++) {
		val = ((spx_value_t) (rand()+i) / ((spx_value_t) RAND_MAX + 1));
		x[i] = min + val*(max-min);
		y[i] = max + val*(min-max);
    }

    spx_init();

    /* Cretae CSR wrapper */
    spx_input_t *input = spx_input_load_csr(rowptr, colind, values, nrows,
                                            ncols, SPX_INDEX_ZERO_BASED);

    /* Transform to CSX */
    spx_option_set("spx.rt.nr_threads", "2");
    spx_option_set("spx.rt.cpu_affinity", "0,1");
    spx_option_set("spx.preproc.xform", "all");
    spx_matrix_t *A = spx_mat_tune(input);

    /* Create x and y vector views */
    spx_partition_t *parts = spx_mat_get_partition(A);
    spx_vector_t *x_view = spx_vec_create_from_buff(x, ncols, parts,
                                                    SPX_VEC_SHARE);
    spx_vector_t *y_view = spx_vec_create_from_buff(y, nrows, parts,
                                                    SPX_VEC_SHARE);

    /* Run the SpMV kernel */
    spx_value_t alpha = 0.8, beta = 0.42;
    spx_matvec_kernel(alpha, A, x_view, beta, y_view);

    /* printf("Result: "); */
    /* spx_vec_print(y_view); */

    /* Cleanup */
    spx_input_destroy(input);
    spx_mat_destroy(A);
    spx_partition_destroy(parts);
    spx_vec_destroy(x_view);
    spx_vec_destroy(y_view);

    return 0;
}
