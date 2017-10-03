/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file sparsex_test.c
 * \brief Simple program for testing the SparseX API
 *
 * \author Athena Elafrou
 * \date 2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/sparsex.h>
#include "CsxCheck.hpp"
#include <stdio.h>
#include <string.h>
#include <getopt.h>

static const char *program_name;

static struct option long_options[] = {
  {"option",              required_argument,  0, 'o'},
  {"enable-timing",       no_argument,        0, 't'},
  {"enable-reordering",   no_argument,        0, 'r'},
  {"help",                no_argument,        0, 'h'}
};

static void print_usage()
{
  fprintf(
	  stderr,
	  "Usage: %s [-o <option=value>]... [-t] <mmf_file>\n\n"
	  "\t-o, --option <option=value>\tset a preprocessing option\n"
	  "\t-r, --enable-timing \t\tenable timing of the SpMV kernel\n"
	  "\t-t, --enable-reordering \tenable reordering of the input matrix\n"   
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
  const size_t loops = 128;

  spx_init();
  spx_log_verbose_console();

  char c;
  char *matrix_file = NULL;
  char *option;
  int option_index = 0;
  int enable_timing = 0;
  int enable_reordering = 0;
  const spx_value_t alpha = 0.5;

  program_name = argv[0];
  while ((c = getopt_long(argc, argv, "o:rth", long_options,
			  &option_index)) != -1) {
    switch (c) {
    case 'o':
      option = optarg;
      set_option(option);
      break;
    case 'r':
      enable_reordering = 1;
      break;
    case 't':
      enable_timing = 1;
      break;
    case 'h':
      print_usage();
      exit(0);
    default:
      print_usage();
      exit(1);
    }
  }
    
  int remargs = argc - optind;
  if (remargs < 1) {
    print_usage();
    exit(1);
  }

  argv = &argv[optind];
  matrix_file = argv[0];

  /* Load matrix from MMF file */
  spx_input_t *input = spx_input_load_mmf(matrix_file);
  if (input == SPX_INVALID_INPUT) {
    SETERROR_0(SPX_ERR_INPUT_MAT);
    exit(1);
  }

  /* Transform to CSX */
  spx_matrix_t *A = SPX_INVALID_MAT;
  if (enable_reordering) {
    A = spx_mat_tune(input, SPX_MAT_REORDER);
  } else {
    A = spx_mat_tune(input);
  }

  if (A == SPX_INVALID_MAT) {
    SETERROR_0(SPX_ERR_TUNED_MAT);
    exit(1);
  }

  spx_partition_t *parts = spx_mat_get_partition(A);
  if (parts == SPX_INVALID_PART) {
    SETERROR_0(SPX_ERR_PART);
    exit(1);
  }

  /* Create random x and y vectors */
  spx_vector_t *x = spx_vec_create_random(spx_mat_get_ncols(A), parts);
  if (x == SPX_INVALID_VEC) {
    SETERROR_0(SPX_ERR_VEC);
    exit(1);
  }

  spx_vector_t *y = spx_vec_create(spx_mat_get_nrows(A), parts);
  if (y == SPX_INVALID_VEC) {
    SETERROR_0(SPX_ERR_VEC);
    exit(1);
  }

  /* Reorder vectors */
  spx_perm_t *p = SPX_INVALID_PERM;
  if (enable_reordering) {
    p = spx_mat_get_perm(A);
    if (p == SPX_INVALID_PERM) {
      SETERROR_0(SPX_ERR_PERM);
      exit(1);
    }
    spx_vec_reorder(x, p);
  }

  /* Run a matrix-vector multiplication: y <-- alpha*A*x */
  spx_timer_t t;
  if (enable_timing) {
    spx_timer_clear(&t);
    spx_timer_start(&t);
  }
    
  for (size_t i = 0; i < loops; i++) {
    spx_matvec_mult(alpha, A, x, y);
  }

  if (enable_timing) {
    spx_timer_pause(&t);
    double elapsed_time = spx_timer_get_secs(&t);
    double flops = (double) (2*loops*spx_mat_get_nnz(A) + 
			     spx_mat_get_nrows(A)) /
      ((double) 1000*1000*elapsed_time);
    printf("SPMV time: %lf secs\n", elapsed_time);
    printf("MFLOPS: %lf\n", flops);
  }

  /* Restore original ordering of resulting vector */
  if (enable_reordering) {
    spx_vec_inv_reorder(y, p);
    spx_vec_inv_reorder(x, p);
  }

  /* Check the result */
  check_result(y, alpha, x, matrix_file);

  /* Cleanup */
  spx_input_destroy(input);
  spx_mat_destroy(A);
  spx_partition_destroy(parts);
  spx_vec_destroy(x);
  spx_vec_destroy(y);
  return 0;
}
