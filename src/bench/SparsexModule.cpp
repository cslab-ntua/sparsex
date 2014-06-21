/*
 * SparsexModule.cpp --  Implementation of the SpMV kernel with SparseX.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "SparsexModule.hpp"

#include <algorithm>
#include <vector>

/* SpMV kernel implemented with SparseX */
void sparsex_spmv(int *rowptr, int *colind, double *values, int nrows, int ncols,
                  int nnz, double *x, double *y)
{
    spx_init();
    /* 1. Matrix loading phase */
    spx_input_t *input = spx_input_load_csr(rowptr, colind, values, nrows, ncols,
                                            INDEXING_ZERO_BASED);

    /* 2. Tuning phase */
    // spx_options_set_from_env();
    spx_option_set("spx.rt.nr_threads", "2");
    spx_option_set("spx.rt.cpu_affinity", "0,1");
    spx_option_set("spx.preproc.xform", "all");
    spx_option_set("spx.preproc.sampling", "portion");
    spx_option_set("spx.preproc.sampling.nr_samples", "48");
    spx_option_set("spx.preproc.sampling.portion", "0.01");
    // spx_option_set("spx.matrix.symmetric", "true");

    t.Clear();
    t.Start();
    spx_matrix_t *A = spx_mat_tune(input);
    t.Pause();
    double pt = t.ElapsedTime();

    /* 3. Vector loading */
    spx_partition_t *parts = spx_mat_get_partition(A);
    spx_vector_t *x_view = spx_vec_create_from_buff(x, ncols, parts, OP_SHARE);
    spx_vector_t *y_view = spx_vec_create_from_buff(y, nrows, parts, OP_SHARE);

    /* 4. SpMV benchmarking phase */
    vector<double> mt(OUTER_LOOPS);
    for (unsigned int i = 0; i < OUTER_LOOPS; i++) {
        t.Clear();
        t.Start();
        for (unsigned long int j = 0; j < LOOPS; j++) {
            spx_matvec_kernel(ALPHA, A, x_view, BETA, y_view);
        }
        t.Pause();
        mt[i] = t.ElapsedTime();
    }
    sort(mt.begin(), mt.end());
    double mt_median = 
        (OUTER_LOOPS % 2) ? mt[((OUTER_LOOPS+1)/2)-1]
        : ((mt[OUTER_LOOPS/2] + mt[OUTER_LOOPS/2+1])/2);  
    double flops = (double)(LOOPS*nnz*2)/((double)1000*1000*mt_median);
    cout << "m: " << MATRIX
         << " pt: " << pt
         << " mt(median): " << mt_median
         << " flops: " << flops << endl;
    
    /* 5. Cleanup */
    spx_input_destroy(input);
    spx_mat_destroy(A);
    spx_partition_destroy(parts);
    spx_vec_destroy(x_view);
    spx_vec_destroy(y_view);
}
