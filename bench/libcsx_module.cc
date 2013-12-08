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
#include <algorithm>
#include <vector>

/* SpMV kernel implemented with LIBCSX */
void libcsx_spmv(int *rowptr, int *colind, double *values, int nrows, int ncols,
                 int nnz, double *x, double *y)
{
    libcsx_init();
    /* 1. Matrix loading phase */
    input_t *input = libcsx_mat_create_csr(rowptr, colind, values, nrows, ncols,
                                           INDEXING_ZERO_BASED);

    /* 2. Tuning phase */
    libcsx_set_options_from_env();
    // libcsx_set_option("libcsx.matrix.symmetric", "true");
    t.Clear();
    t.Start();
    matrix_t *A = libcsx_mat_tune(input);
    t.Pause();
    double pt = t.ElapsedTime();

    /* 3. Vector loading */
    vector_t *x_view = libcsx_vec_create_from_buff(x, ncols, A);
    vector_t *y_view = libcsx_vec_create_from_buff(y, nrows, A);

    /* 4. SpMV benchmarking phase */
    std::vector<double> mt(OUTER_LOOPS);
    SPMV_BENCH(libcsx_matvec_mult(A, ALPHA, x_view, BETA, y_view));
    sort(mt.begin(), mt.end());
    double mt_median = 
        (OUTER_LOOPS % 2) ? mt[((OUTER_LOOPS+1)/2)-1]
        : ((mt[OUTER_LOOPS/2] + mt[OUTER_LOOPS/2+1])/2);  
    double flops = (double)(LOOPS*nnz*2)/((double)1000*1000*mt_median);
    cout << "m: " << MATRIX
         << " pt: " << pt
         << " mt(median): " << mt_median
         << " flops: " << flops << endl;
    // for (size_t i=0;i<OUTER_LOOPS;i++) {
    //     cout << "m: " << MATRIX
    //          << " pt: " << pt
    //          << " mt: " << mt[i]
    //          << " flops: " << (double)(LOOPS*nnz*2)/((double)1000*1000*mt[i]) << endl;
    // }

    /* 5. Cleanup */
    libcsx_mat_destroy_input(input);
    libcsx_mat_destroy_tuned(A);
    libcsx_vec_destroy(x_view);
    libcsx_vec_destroy(y_view);
}
