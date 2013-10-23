/* -*- C++ -*-
 *
 * libcsx_module.h --  Implementation of the SpMV kernel with LIBCSX.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "libcsx_module.h"

/* SpMV kernel implemented with LIBCSX */
void libcsx_spmv(int *rowptr, int *colind, double *values, int nrows, int ncols,
                 int nnz, double *x, double *y)
{
    libcsx_init();
    /* 1. Matrix loading phase */
    input_t *input = libcsx_mat_create_csr(rowptr, colind, values, nrows, ncols,
                                           true);

    /* 2. Vector loading */
    vector_t *x_view = vec_create_from_buff(x, ncols);
    vector_t *y_view = vec_create_from_buff(y, nrows);

    /* 3. Tuning phase */
    t.Clear();
    t.Start();
    libcsx_set_option("libcsx.rt.cpu_affinity", "0,1");
    // libcsx_set_option("libcsx.matrix.symmetric", "true");
    libcsx_set_option("libcsx.preproc.xform", "all");
    libcsx_set_option("libcsx.preproc.sampling", "portion");
    libcsx_set_option("libcsx.preproc.sampling.portion", ".01");
    libcsx_set_option("libcsx.preproc.sampling.nr_samples", "48");
    matrix_t *A = libcsx_mat_tune(input, 0);
    t.Pause();
    double pt = t.ElapsedTime();
    cout << "pt: " << pt << endl;
    
    /* 4. SpMV benchmarking phase */
    SPMV_BENCH(libcsx_matvec_mult(A, ALPHA, x_view, BETA, y_view));
    double mt = t.ElapsedTime() / OUTER_LOOPS;
    cout << "mt: " << mt << endl;
//    vec_print(y_view);

    /* 5. Cleanup */
    libcsx_mat_destroy_input(input);
    libcsx_mat_destroy_tuned(A);
    vec_destroy(x_view);
    vec_destroy(y_view);
}
