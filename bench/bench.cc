/* -*- C++ -*-
 *
 * bench.cc --  Benchmarking interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "bench.h"
#include "mmf.h"
#include <cassert>
#include <libgen.h>

using namespace std;
using namespace bench;

/* Global variables */
unsigned int OUTER_LOOPS;    /**< Number of SpMV iterations */
unsigned long LOOPS;         /**< Number of repeats */
unsigned int NR_THREADS = 2; /**< Number of threads for a multithreaded
                                execution */
double ALPHA = 2, BETA = 1;  /**< Scalar parameters of the SpMV kernel
                                (y->APLHA*A*x + BETA*y) */
Timer t;                     /**< Timer for benchmarking */

static SpmvFn GetSpmvFn(library type);
static void MMFtoCSR(const char *filename, int **rowptr, int **colind,
                     double **values, int *nrows, int *ncols, int *nnz);


void Bench_Directory(const char *directory, const char *library,
                     const char *stats_file, unsigned long loops,
                     unsigned int outer_loops)
{
    DirectoryIterator it(directory);
 
    for (; it; ++it) {
        const char *filename = it.filename().c_str();
        Bench_Matrix(filename, library, stats_file, loops, outer_loops);
    }
}

void Bench_Matrix(const char *filename, const char *library,
                  const char *stats_file, unsigned long loops,
                  unsigned int outer_loops)
{
    cout << "\n==========================================================\n";
    cout << " COMPUTING SPMV PRODUCT WITH MATRIX: " 
         << basename(const_cast<char *>(filename)) << endl;

    LOOPS = loops;
    OUTER_LOOPS = outer_loops;
    SpmvFn  fn;
    int     *rowptr, *colind;
    double  *values;
    int     nrows, ncols, nnz;

    MMFtoCSR(filename, &rowptr, &colind, &values, &nrows, &ncols, &nnz);

	double *x = (double *) malloc(sizeof(double) * ncols);
	double *y = (double *) malloc(sizeof(double) * nrows);
    for (int i = 0; i < nrows; i++) {
        x[i] = 1;
        y[i] = 2;
    }

    if (!library) {
        cout << "Using library Intel MKL..." << endl;
        mkl_spmv(rowptr, colind, values, nrows, ncols, nnz, x, y);            
        cout << "Using library pOSKI..." << endl;
        poski_spmv(rowptr, colind, values, nrows, ncols, nnz, x, y);            
        cout << "Using library LIBCSX..." << endl;
        libcsx_spmv(rowptr, colind, values, nrows, ncols, nnz, x, y); 
    } else {
        cout << "Using library " << library << "...\n";
        if (strcmp(library, "MKL") == 0) {
            fn = GetSpmvFn(MKL);
        } else if (strcmp(library, "pOSKI") == 0) {
            fn = GetSpmvFn(pOSKI);
        } else if (strcmp(library, "LIBCSX") == 0) {
            fn = GetSpmvFn(LIBCSX);
        } else {
            cerr << "Library doesn't exist." << endl;
        }
        fn(rowptr, colind, values, nrows, ncols, nnz, x, y);            
    }

    /* Cleanup */
    free(rowptr);
    free(colind);
    free(values);
    free(x);
    free(y);
}

void Bench_Matrix(const char *mmf_file, SpmvFn fn, const char *stats_file,
                  unsigned long loops, unsigned int outer_loops)
{
    cout << "\n==========================================================\n";
    cout << " COMPUTING SPMV PRODUCT WITH MATRIX: " 
         << basename(const_cast<char *>(mmf_file)) << endl;

    LOOPS = loops;
    OUTER_LOOPS = outer_loops;
    int     *rowptr, *colind;
    double  *values;
    int     nrows, ncols, nnz;

    MMFtoCSR(mmf_file, &rowptr, &colind, &values, &nrows, &ncols, &nnz);

	double *x = (double *) malloc(sizeof(double) * ncols);
	double *y = (double *) malloc(sizeof(double) * nrows);
    for (int i = 0; i < nrows; i++) {
        x[i] = 1;
        y[i] = 2;
    }
    
    fn(rowptr, colind, values, nrows, ncols, nnz, x, y);            

    /* Cleanup */
    free(rowptr);
    free(colind);
    free(values);
    free(x);
    free(y);
}

static SpmvFn GetSpmvFn(library type)
{
    SpmvFn ret;

    switch(type) {
    case LIBCSX:
        ret = libcsx_spmv;
        break;
    case MKL:
        ret = mkl_spmv;
        break;
    case pOSKI:
        ret = poski_spmv;
        break;
    default:
        cerr << "Unknown library" << endl;
        assert(false);
    }

    return ret;
}

static void MMFtoCSR(const char *filename, int **rowptr, int **colind,
                     double **values, int *nrows, int *ncols, int *nnz)
{
    MMF<int, double> mmf(filename);
    *nrows = mmf.GetNrRows(); *ncols = mmf.GetNrCols();
    *nnz = mmf.GetNrNonzeros();
	*values = (double *) malloc(sizeof(double) * mmf.GetNrNonzeros());
	*colind = (int *) malloc(sizeof(int) * mmf.GetNrNonzeros());
	*rowptr = (int *) malloc(sizeof(int) * (mmf.GetNrRows() + 1));

    MMF<int, double>::iterator iter = mmf.begin();
    MMF<int, double>::iterator iter_end = mmf.end();   
    int row_i = 0, val_i = 0, row_prev = 0;
    int row, col;
    double val;

	(*rowptr)[row_i++] = val_i;
    for (;iter != iter_end; ++iter) {
        row = (*iter).row - 1;
        col = (*iter).col - 1;
        val = (*iter).val;
		assert(row >= row_prev);
		if (row != row_prev) {
			for (int i = 0; i < row - row_prev; i++) {
				(*rowptr)[row_i++] = val_i;
            }
			row_prev = row;
		}
		(*values)[val_i] = val;
		(*colind)[val_i] = col;
		val_i++;
    }
	(*rowptr)[row_i++] = val_i;
}
