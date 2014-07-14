/*
 * PoskiModule.hpp --  The SpMV kernel with pOSKI.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef POSKI_MODULE_HPP
#define POSKI_MODULE_HPP

#include <sparsex/internals/cdecl.h>

SPX_BEGIN_C_DECLS__
#include <poski/poski.h>
SPX_END_C_DECLS__

#include "Timer.hpp"
#include <sparsex/types.h>
#include <iostream>

using namespace std;

extern string MATRIX; 
extern unsigned int OUTER_LOOPS;
extern unsigned long LOOPS;
extern unsigned int NR_THREADS;
extern spx_value_t ALPHA, BETA;
extern Timer t;

void poski_spmv(spx_index_t *Aptr, spx_index_t *Aind, spx_value_t *Aval,
                spx_index_t nrows, spx_index_t ncols,
                spx_index_t nnz, spx_value_t *x, spx_value_t *y);

#endif  // POSKI_MODULE_HPP

