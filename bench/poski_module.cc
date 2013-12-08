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

#include "poski_module.h"
#include <algorithm>
#include <vector>

/* SpMV kernel implemented with pOSKI */
void poski_spmv(int *Aptr, int *Aind, double *Aval, int nrows, int ncols,
                int nnz, double *x, double *y)
{
    poski_Init();

    /* 1. Matrix loading phase */
    poski_threadarg_t *poski_thread = poski_InitThreads();
    poski_ThreadHints(poski_thread, NULL, POSKI_THREADPOOL, NR_THREADS);
    poski_partitionarg_t *partitionMat = poski_partitionMatHints(OneD, NR_THREADS,
                                                                 KERNEL_MatMult,
                                                                 OP_NORMAL);
    poski_mat_t A_tunable = 
        poski_CreateMatCSR(Aptr, Aind, Aval, nrows, ncols, nnz, SHARE_INPUTMAT,
                           poski_thread, partitionMat, 2, INDEX_ZERO_BASED, 
                           MAT_GENERAL);

    /* 2. Vector loading */
    poski_partitionvec_t *partitionVecX = 
        poski_PartitionVecHints(A_tunable, KERNEL_MatMult, OP_NORMAL, INPUTVEC);
    poski_vec_t x_view = poski_CreateVec(x, ncols, STRIDE_UNIT, partitionVecX);
    poski_partitionvec_t *partitionVecY = 
        poski_PartitionVecHints(A_tunable, KERNEL_MatMult, OP_NORMAL, OUTPUTVEC);
    poski_vec_t y_view = poski_CreateVec(y, nrows, STRIDE_UNIT, partitionVecY);

    /* 3. Tuning phase */
    t.Clear();
    t.Start();
    poski_TuneHint_MatMult(A_tunable, OP_NORMAL, ALPHA, x_view, BETA, 
                           y_view, ALWAYS_TUNE_AGGRESSIVELY);
    poski_TuneMat(A_tunable);
    t.Pause();
    double pt = t.ElapsedTime();

    /* 4. SpMV benchmarking phase */
    std::vector<double> mt(OUTER_LOOPS);
    SPMV_BENCH(poski_MatMult(A_tunable, OP_NORMAL, ALPHA, x_view, BETA, 
                             y_view));
    sort(mt.begin(), mt.end());
    // for (size_t i=0;i<OUTER_LOOPS;i++) {
    //     cout << "m: " << MATRIX
    //          << " mt: " << mt[i]
    //          << " flops: " << (double)(LOOPS*nnz*2)/((double)1000*1000*mt[i]) << endl;
    // }
    double mt_median = 
        (OUTER_LOOPS % 2) ? mt[((OUTER_LOOPS+1)/2)-1]
        : ((mt[OUTER_LOOPS/2] + mt[OUTER_LOOPS/2+1])/2);  
    double flops = (double)(LOOPS*nnz*2)/((double)1000*1000*mt_median);
    cout << "m: " << MATRIX 
         << " pt: " << pt 
         << " mt(median): " << mt_median
         << " flops: " << flops << endl;

    /* 4. Cleanup */
    poski_DestroyPartitionVec(partitionVecX);
    poski_DestroyPartitionVec(partitionVecY);
    poski_DestroyVec(x_view);
    poski_DestroyVec(y_view);
    poski_DestroyPartitionMat(partitionMat);
    poski_DestroyMat(A_tunable);
    poski_DestroyThreads(poski_thread);
    poski_Close();
}
