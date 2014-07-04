/*
 * Bench.cpp --  Benchmarking interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "Bench.hpp"
#include "Mmf.hpp"

#include <cassert>
#include <libgen.h>

using namespace std;
using namespace bench;

/* Global variables */
string MATRIX;
unsigned int OUTER_LOOPS = 5;   /**< Number of SpMV iterations */
unsigned long LOOPS = 128;      /**< Number of repeats */
unsigned int NR_THREADS = 1;    /**< Number of threads for a multithreaded
                                   execution */
spx_value_t ALPHA = 1.32, BETA = 0.48;     /**< Scalar parameters of the SpMV kernel
                                        (y->APLHA*A*x + BETA*y) */
Timer t;                        /**< Timer for benchmarking */

static SpmvFn GetSpmvFn(library type);
static void MMFtoCSR(const char *filename,
                     spx_index_t **rowptr, spx_index_t **colind,
                     spx_value_t **values,
                     spx_index_t *nrows, spx_index_t *ncols, spx_index_t *nnz);
static inline int elems_neq(spx_value_t a, spx_value_t b);
static int vec_compare(const spx_value_t *v1, const spx_value_t *v2, size_t size)
    __attribute__ ((unused));

void Bench_Directory(const char *directory, const char *library,
                     const char *stats_file)
{
    DirectoryIterator it(directory);
 
    for (; it; ++it) {
        const char *filename = it.filename().c_str();
        Bench_Matrix(filename, library, stats_file);
    }
}

void Bench_Matrix(const char *filename, const char *library,
                  const char *stats_file)
{
    MATRIX = basename(const_cast<char *>(filename));
    cout << "\n==========================================================\n";
    cout << " COMPUTING SPMV PRODUCT WITH MATRIX: " << MATRIX << endl;

    SpmvFn  fn;
    spx_index_t *rowptr, *colind;
    spx_value_t *values;
    spx_index_t nrows, ncols, nnz;

    MMFtoCSR(filename, &rowptr, &colind, &values, &nrows, &ncols, &nnz);

	spx_value_t *x = (spx_value_t *) malloc(sizeof(spx_value_t) * ncols);
	spx_value_t *y = (spx_value_t *) malloc(sizeof(spx_value_t) * nrows);
#if defined(MKL) || defined(POSKI)
	spx_value_t *y_cmp = (spx_value_t *) malloc(sizeof(spx_value_t) * nrows);
#endif
    spx_value_t val = 0, max = 1, min = -1;
 
    for (spx_index_t i = 0; i < nrows; i++) {
		val = ((spx_value_t) (rand()+i) / ((spx_value_t) RAND_MAX + 1));
		x[i] = min + val*(max-min);
		y[i] = max + val*(min-max);
#if defined(MKL) || defined(POSKI)
        y_cmp[i] = y[i];
#endif
    }

    ALPHA = min + val*(max-min);
    BETA = max + val*(min-max);

    if (!library) {
#ifdef MKL
        cout << "Using library Intel MKL..." << endl;
        mkl_spmv(rowptr, colind, values, nrows, ncols, nnz, x, y_cmp);            
#endif
#ifdef POSKI
        cout << "Using library pOSKI..." << endl;
        poski_spmv(rowptr, colind, values, nrows, ncols, nnz, x, y_cmp);
#endif
        cout << "Using library SparseX..." << endl;
        sparsex_spmv(rowptr, colind, values, nrows, ncols, nnz, x, y);

#ifdef COMPARE_RESULT
        if (vec_compare(y, y_cmp, nrows) < 0)
            cout << "Error in resulting vector!" << endl;
        else
            cout << "Checked passed!" << endl;
#endif
    } else {
        cout << "Using library " << library << "...\n";
        if (strcmp(library, "MKL") == 0) {
            fn = GetSpmvFn(MKL);
        } else if (strcmp(library, "pOSKI") == 0) {
            fn = GetSpmvFn(pOSKI);
        } else if (strcmp(library, "SparseX") == 0) {
            fn = GetSpmvFn(SparseX);
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
#if defined(MKL) || defined(POSKI)
    free(y_cmp);
#endif
}

void Bench_Matrix(const char *mmf_file, SpmvFn fn, const char *stats_file)
{
    cout << "\n==========================================================\n";
    cout << " COMPUTING SPMV PRODUCT WITH MATRIX: " 
         << basename(const_cast<char *>(mmf_file)) << endl;

    spx_index_t *rowptr, *colind;
    spx_value_t *values;
    spx_index_t nrows, ncols, nnz;

    MMFtoCSR(mmf_file, &rowptr, &colind, &values, &nrows, &ncols, &nnz);

	spx_value_t *x = (spx_value_t *) malloc(sizeof(spx_value_t) * ncols);
	spx_value_t *y = (spx_value_t *) malloc(sizeof(spx_value_t) * nrows);
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
    case SparseX:
        ret = sparsex_spmv;
        break;
#ifdef MKL
    case MKL:
        ret = mkl_spmv;
        break;
#endif
#ifdef POSKI
    case pOSKI:
        ret = poski_spmv;
        break;
#endif
    default:
        cerr << "Unknown library" << endl;
        assert(false);
    }

    return ret;
}

static void MMFtoCSR(const char *filename,
                     spx_index_t **rowptr, spx_index_t **colind,
                     spx_value_t **values,
                     spx_index_t *nrows, spx_index_t *ncols, spx_index_t *nnz)
{
    MMF<spx_index_t, spx_value_t> mmf(filename);
    *nrows = mmf.GetNrRows(); *ncols = mmf.GetNrCols();
    *nnz = mmf.GetNrNonzeros();
	*values = (spx_value_t *) malloc(sizeof(spx_value_t) * mmf.GetNrNonzeros());
	*colind = (spx_index_t *) malloc(sizeof(**colind) * mmf.GetNrNonzeros());
	*rowptr = (spx_index_t *) malloc(sizeof(**rowptr) * (mmf.GetNrRows() + 1));

    MMF<spx_index_t, spx_value_t>::iterator iter = mmf.begin();
    MMF<spx_index_t, spx_value_t>::iterator iter_end = mmf.end();   
    spx_index_t row_i = 0, val_i = 0, row_prev = 0;
    spx_index_t row, col;
    spx_value_t val;

	(*rowptr)[row_i++] = val_i;
    for (;iter != iter_end; ++iter) {
        row = (*iter).row - 1;
        col = (*iter).col - 1;
        val = (*iter).val;
		assert(row >= row_prev);
		if (row != row_prev) {
			for (spx_index_t i = 0; i < row - row_prev; i++) {
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

static inline int elems_neq(spx_value_t a, spx_value_t b)
{
	if (fabs((spx_value_t) (a - b) / (spx_value_t) a)  > 1.e-7)
		return 1;
	return 0;
}

static int vec_compare(const spx_value_t *v1,
                       const spx_value_t *v2, size_t size)
{
	for (size_t i = 0; i < size; i++) {
		if (elems_neq(v1[i], v2[i])) {
            cout << "element " << i << " differs: "<< v1[i] << " != "
                 << v2[i] << endl;
			return -1;
		}
	}

	return 0;
}
