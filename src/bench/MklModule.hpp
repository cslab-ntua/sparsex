/*
 * MklModule.hpp --  The SpMV kernel with Intel MKL.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef MKL_MODULE_HPP
#define MKL_MODULE_HPP

#include "Timer.hpp"
#include <sparsex/types.h>
#include <mkl.h>
#include <iostream>

using namespace std;

extern string MATRIX; 
extern unsigned int OUTER_LOOPS;
extern unsigned long LOOPS;
extern unsigned int NR_THREADS;
extern spx_value_t ALPHA, BETA;
extern Timer t;

void mkl_spmv(spx_index_t *rowptr, spx_index_t *colind,
              spx_value_t *values, spx_index_t nrows, spx_index_t ncols,
              spx_index_t nnz, spx_value_t *x, spx_value_t *y);

#endif  // MKL_MODULE_HPP

