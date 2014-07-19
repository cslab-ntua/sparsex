/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file MklModule.cpp
 * \brief The SpMV kernel implemented with Intel MKL
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include "MklModule.hpp"

#include <algorithm>
#include <vector>
#include <sched.h>
#include <stdio.h>

using namespace std;

extern string MATRIX; 
extern unsigned int OUTER_LOOPS;
extern unsigned long LOOPS;
extern Timer t;

/* SpMV kernel implemented with Intel MKL */
void mkl_spmv(spx_index_t *rowptr, spx_index_t *colind, spx_value_t *values,
              spx_index_t nrows, spx_index_t ncols, spx_index_t nnz,
              spx_value_t *x, spx_value_t *y,
              spx_value_t ALPHA, spx_value_t BETA)
{
    /* 1. Matrix loading phase */
    MKL_INT *pointerB, *pointerE;
    char    transa;
    char    matdescra[6];

    transa = 'n';
    matdescra[0] = 'g';
    matdescra[1] = '-';
    matdescra[2] = '-';
    matdescra[3] = 'c';

	pointerB = (MKL_INT *) malloc(sizeof(MKL_INT) * nrows);
	pointerE = (MKL_INT *) malloc(sizeof(MKL_INT) * nrows);
    for (int i = 0; i < nrows; i++) {
        pointerB[i] = rowptr[i];
        pointerE[i] = rowptr[i+1];
    }
    
    /* 2. SpMV benchmarking phase */
    vector<double> mt(OUTER_LOOPS);
    for (size_t i = 0; i < OUTER_LOOPS; i++) {
        t.Clear();
        t.Start();
        for (size_t j = 0; j < LOOPS; j++) {
            mkl_dcsrmv(&transa, &nrows, &ncols, &ALPHA, matdescra, values,
                       colind, pointerB, pointerE, x, &BETA, y);            
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
         << " mt(median): " << mt_median
         << " flops: " << flops << endl;

    /* 3. Cleanup */
    free(pointerB);
    free(pointerE);
}
