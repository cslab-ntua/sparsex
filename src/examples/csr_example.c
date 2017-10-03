/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file csr_example.c
 * \brief Read matrix from the CSR format and tune it into the CSX format. This
 *  example uses the spx_mat_vec_kernel_csr() routine that hides the
 *  preprocessing phase of CSX.
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

  /* Initialize library */
  spx_init();
  spx_log_info_console();

  /* Partition input matrix */
  spx_partition_t *parts = spx_partition_csr(rowptr, nrows, 2);

  /* Create random x and y vectors */
  spx_vector_t *x = spx_vec_create_random(nrows, parts);
  spx_vector_t *y = spx_vec_create_random(nrows, parts);

  /* Declare a tuned matrix handle */
  spx_matrix_t *A = SPX_INVALID_MAT;

  /* Set tuning options */
  spx_option_set("spx.rt.nr_threads", "2");
  spx_option_set("spx.rt.cpu_affinity", "0,1");

  /* Transform to CSX and run the SpMV kernel */
  spx_value_t alpha = 0.8, beta = 0.42;
  const size_t nr_loops = 1500;
  size_t i;
  spx_timer_t t;
  double elapsed_time, flops;

  spx_timer_clear(&t);
  spx_timer_start(&t);
  for (i = 0; i < nr_loops; i++) {
    spx_matvec_kernel_csr(&A, nrows, ncols, rowptr, colind, values,
			  alpha, x, beta, y);
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
