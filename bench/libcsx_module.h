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
#ifndef LIBCSX_MODULE_H__
#define LIBCSX_MODULE_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "libcsx.h"

#ifdef __cplusplus
}
#endif

#include <iostream>
#include "macros.h"

using namespace std;

extern std::string MATRIX; 
extern unsigned int OUTER_LOOPS;
extern unsigned long LOOPS;
extern unsigned int NR_THREADS;
extern double ALPHA, BETA;
extern Timer t;

void libcsx_spmv(int *rowptr, int *colind, double *values, int nrows, int ncols,
                 int nnz, double *x, double *y);

#endif  // LIBCSX_MODULE_H__
