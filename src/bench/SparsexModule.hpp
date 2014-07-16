/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file SparsexModule.hpp
 * \brief The SpMV kernel implemented with SparseX
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_MODULE_HPP
#define SPARSEX_MODULE_HPP

#ifdef __cplusplus
extern "C" {
#endif

#include <sparsex/sparsex.h>

#ifdef __cplusplus
}
#endif

#include "Timer.hpp"
#include <iostream>

using namespace std;

extern string MATRIX; 
extern unsigned int OUTER_LOOPS;
extern unsigned long LOOPS;
extern unsigned int NR_THREADS;
extern spx_value_t ALPHA, BETA;
extern Timer t;

void sparsex_spmv(spx_index_t *rowptr, spx_index_t *colind,
                  spx_value_t *values, spx_index_t nrows,
                  spx_index_t ncols, spx_index_t nnz,
                  spx_value_t *x, spx_value_t *y);

#endif  // SPARSEX_MODULE_HPP
