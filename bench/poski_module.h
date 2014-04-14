/* -*- C++ -*-
 *
 * poski_module.h --  The SpMV kernel with pOSKI.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef POSKI_MODULE_H__
#define POSKI_MODULE_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <poski.h>

#ifdef __cplusplus
}
#endif

#include <iostream>
#include "timer.h"

using namespace std;

extern std::string MATRIX; 
extern unsigned int OUTER_LOOPS;
extern unsigned long LOOPS;
extern unsigned int NR_THREADS;
extern double ALPHA, BETA;
extern Timer t;

void poski_spmv(int *Aptr, int *Aind, double *Aval, int nrows, int ncols,
                int nnz, double *x, double *y);

#endif  // POSKI_MODULE_H__

