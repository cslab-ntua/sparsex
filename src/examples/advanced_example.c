/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file advanced_example.c
 * \brief This example shows usage of the vector tuning feature.
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
  fprintf(stderr, "Usage: %s <mmf_file>\n", program_name);
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

  /* Load matrix from MMF file */
  spx_input_t *input = spx_input_load_mmf(argv[1]);

  /* Set tuning options */
  spx_option_set("spx.rt.nr_threads", "2");
  spx_option_set("spx.rt.cpu_affinity", "0,1");

  /* Transform to CSX */
  spx_matrix_t *A = spx_mat_tune(input);

  /* Define random x and y arrays */
  spx_index_t nrows = spx_mat_get_nrows(A);
  spx_index_t ncols = spx_mat_get_ncols(A);
  spx_value_t *x = (spx_value_t *) malloc(sizeof(spx_value_t) * ncols);
  spx_value_t *y = (spx_value_t *) malloc(sizeof(spx_value_t) * nrows);
  spx_value_t val = 0, max = 1, min = -1;
  size_t i;
  for (i = 0; i < nrows; i++) {
    val = ((spx_value_t) (rand()+i) / ((spx_value_t) RAND_MAX + 1));
    x[i] = min + val*(max-min);
    y[i] = max + val*(min-max);
  }

  /* Create x and y vectors from the corresponding buffers */
  spx_value_t *x_tuned, *y_tuned;
  spx_partition_t *parts = spx_mat_get_partition(A);
  spx_vector_t *x_view = spx_vec_create_from_buff(
						  x, &x_tuned, ncols, parts, SPX_VEC_TUNE);
  spx_vector_t *y_view = spx_vec_create_from_buff(
						  y, &y_tuned, nrows, parts, SPX_VEC_TUNE);

  /* Run 128 loops of the SpMV kernel */
  spx_value_t alpha = 0.8, beta = 0.42;
  const size_t nr_loops = 128;
  spx_timer_t t;
  double elapsed_time, flops;

  spx_timer_clear(&t);
  spx_timer_start(&t);
  for (i = 0; i < nr_loops; i++) {
    spx_matvec_kernel(alpha, A, x_view, beta, y_view);
  }
  spx_timer_pause(&t);
  elapsed_time = spx_timer_get_secs(&t);
  flops = (double) (2 * nr_loops * spx_mat_get_nnz(A)) /
    ((double) 1000 * 1000 * elapsed_time);
  printf("SPMV time: %lf secs\n", elapsed_time);
  printf("MFLOPS: %lf\n", flops);

  /* From this point on the user can use the tuned buffers */
  if (x_tuned != x) {
    x = x_tuned;
  }

  if (y_tuned != y) {
    y = y_tuned;
  }

  /* ... */
  /* ... */
  /* ... */

  /* Cleanup */
  spx_input_destroy(input);
  spx_mat_destroy(A);
  spx_partition_destroy(parts);
  spx_vec_destroy(x_view);
  spx_vec_destroy(y_view);

  /* The user can apply the free() function on the tuned buffers */
  free(x);
  free(y);

  /* Shutdown library */
  spx_finalize();

  return 0;
}
