/* -*- C++ -*-
 *
 * mkl_module.h --  The SpMV kernel with Intel MKL.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef MKL_MODULE_H__
#define MKL_MODULE_H__

#include <iostream>
#include "mkl.h"
#include "timer.h"

using namespace std;

extern std::string MATRIX; 
extern unsigned int OUTER_LOOPS;
extern unsigned long LOOPS;
extern unsigned int NR_THREADS;
extern double ALPHA, BETA;
extern Timer t;

void mkl_spmv(int *rowptr, int *colind, double *values, int nrows, int ncols,
              int nnz, double *x, double *y);

#endif  // MKL_MODULE_H__

