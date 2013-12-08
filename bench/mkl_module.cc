/* -*- C++ -*-
 *
 * mkl_module.cc --  The SpMV kernel with Intel MKL.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "mkl_module.h"
#include <algorithm>
#include <vector>
#include <sched.h>
#include <stdio.h>
#include <omp.h>

/* SpMV kernel implemented with Intel MKL */
void mkl_spmv(int *rowptr, int *colind, double *values, int nrows, int ncols,
              int nnz, double *x, double *y)
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
    mkl_set_num_threads(NR_THREADS);                                                                                                                                                                                              
    /* 2. SpMV benchmarking phase */
    std::vector<double> mt(OUTER_LOOPS);
    /*#pragma omp parallel num_threads(NR_THREADS) default(shared)
    {
        cpu_set_t new_mask;
        int tid = omp_get_thread_num();
        CPU_ZERO(&new_mask);
        CPU_SET(tid%4, &new_mask);
        if (sched_setaffinity(0, sizeof(new_mask), &new_mask) == -1) {
            printf("Error: sched_setaffinity(%d, sizeof(new_mask), &new_mask)\n", tid);
        }
        printf("tid=%d new_mask=%08X\n", tid, *(unsigned int*)(&new_mask));
    }*/

    SPMV_BENCH(mkl_dcsrmv(&transa, &nrows, &ncols, &ALPHA, matdescra, values,
                           colind, pointerB, pointerE, x, &BETA, y));            
    sort(mt.begin(), mt.end());
    double mt_median = 
        (OUTER_LOOPS % 2) ? mt[((OUTER_LOOPS+1)/2)-1]
        : ((mt[OUTER_LOOPS/2] + mt[OUTER_LOOPS/2+1])/2);  
    double flops = (double)(LOOPS*nnz*2)/((double)1000*1000*mt_median);
    cout << "m: " << MATRIX
         << " mt(median): " << mt_median
         << " flops: " << flops << endl;
    // for (size_t i=0;i<OUTER_LOOPS;i++) {
    //     cout << "m: " << MATRIX
    //          << " mt: " << mt[i]
    //          << " flops: " << (double)(LOOPS*nnz*2)/((double)1000*1000*mt[i]) << endl;
    // }
    std::cout << y[0] << " " << y[nrows-64]<<std::endl;

    /* 3. Cleanup */
    free(pointerB);
    free(pointerE);
}
