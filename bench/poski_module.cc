/* -*- C++ -*-
 *
 * poski_module.h --  Implementation of the SpMV kernel with pOSKI.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "poski_module.h"

/* SpMV kernel implemented with pOSKI */
void poski_spmv(int *Aptr, int *Aind, double *Aval, int nrows, int ncols,
                int nnz, double *x, double *y)
{
    poski_Init();

    /* 1. Matrix loading phase */
    poski_threadarg_t *poski_thread = poski_InitThreads();
    poski_mat_t A_tunable = 
        poski_CreateMatCSR(Aptr, Aind, Aval, nrows, ncols, nnz, SHARE_INPUTMAT,
                           poski_thread, NULL, 2, INDEX_ZERO_BASED, MAT_GENERAL);

    /* 2. Vector loading */
    poski_vec_t x_view = poski_CreateVec(x, ncols, STRIDE_UNIT, NULL);
    poski_vec_t y_view = poski_CreateVec(y, nrows, STRIDE_UNIT, NULL);

    /* 3. Tuning phase */
    poski_MatMult(A_tunable, OP_NORMAL, ALPHA, x_view, BETA, y_view);
    t.Clear();
    t.Start();
    poski_TuneHint_MatMult(A_tunable, OP_NORMAL, ALPHA, x_view, BETA, 
                           y_view, LOOPS);
    poski_TuneMat(A_tunable);
    t.Pause();
    double pt = t.ElapsedTime();
    cout << "pt: " << pt << endl;

    /* 4. SpMV benchmarking phase */
    SPMV_BENCH(poski_MatMult(A_tunable, OP_NORMAL, ALPHA, x_view, BETA, 
                             y_view));
    double mt = t.ElapsedTime() / OUTER_LOOPS;
    cout << "mt: " << mt << endl;
    // for (int i = 0; i < nrows; i++) {
    //     std::cout << y[i] << " ";
    // }
    // cout << endl;

    /* 4. Cleanup */
    poski_DestroyMat(A_tunable);
    poski_DestroyVec(x_view);
    poski_DestroyVec(y_view);
    poski_DestroyThreads(poski_thread);
    poski_Close();
}
