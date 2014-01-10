/*
 * elmerif.cc -- Interface for Elmer integration
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <iostream>
#include <fstream>
#include "spm.h"
#include "spmv.h"

extern "C" {
#include "elmerif.h"
#include "mt_lib.h"
#include "spmv_matvec_mt.h"
#include "string.h"
}

template<typename IndexType, typename ValueType>
void write_mmf(std::ostream &os, IndexType *rowptr, IndexType *colind,
               ValueType *values, IndexType nr_rows)
{
    IndexType nr_nzeros = rowptr[nr_rows] - 1;
    os << nr_rows << " " << nr_rows << " " << nr_nzeros << std::endl;

    for (IndexType i = 0; i < nr_rows; ++i) {
        for (IndexType j = rowptr[i] - 1; j < rowptr[i+1] - 1; ++j)
            os << i+1 << " " << colind[j] << " " << values[j] << std::endl;
    }
}

template<typename MatrixType>
void compare(MatrixType mat, std::vector<csx::CooElem> &elems)
{
    typename MatrixType::iterator iter = mat.begin();
    typename MatrixType::iterator iter_end = mat.end();
    std::vector<csx::CooElem>::iterator elem_iter = elems.begin();

    assert(mat.GetNrNonzeros() == elems.size());
    for (; iter != iter_end; ++iter, ++elem_iter) {
        csx::CooElem mat_elem = *iter;
        csx::CooElem target_elem = *elem_iter;
        assert(!CooCmp(mat_elem, target_elem) &&
               mat_elem.val == target_elem.val);
    }
}

template<typename T>
static void print_array(T *array, size_t n)
{
    std::cout << "[ ";
    for (size_t i = 0; i < n; i++)
        std::cout << array[i] << " ";
    std::cout << "]" << std::endl;
}

static void print_prod(const elmer_index_t *rowptr,
                       const elmer_index_t *colind,
                       const elmer_value_t *values,
                       const elmer_value_t *x,
                       elmer_index_t row_no)
{
    elmer_index_t i;
    register elmer_value_t y_ = 0;
    for (i = rowptr[row_no] - 1; i < rowptr[row_no+1] - 1; ++i) {
        std::cout << "i = " << i << ", "
                  << "x = " << colind[i] << ","
                  << "v = " << values[i] << "|";
        y_ += values[i]*x[colind[i]-1];
    }

    std::cout << "y = " << y_ << "|" << std::endl;
}

static void matvec_csr(const elmer_index_t *rowptr,
                       const elmer_index_t *colind,
                       const elmer_value_t *values,
                       elmer_index_t nr_rows,
                       const elmer_value_t *x,
                       elmer_value_t *y)
{
    elmer_index_t i, j;
    for (i = 0; i < nr_rows; ++i) {
        register elmer_value_t yi_ = 0;
        for (j = rowptr[i] - 1; j < rowptr[i+1] - 1; ++j) {
            yi_ += values[j]*x[colind[j] - 1];
        }

        y[i] = yi_;
    }
}

struct csr_struct {
    elmer_index_t   *rowptr;
    elmer_index_t   *colind;
    elmer_value_t   *values;
    elmer_index_t   nr_rows;
    elmer_index_t   nr_nzeros;
};

static unsigned long nr_calls = 0;
static struct csr_struct csr_orig;

void elmer_matvec_(void **tuned, void *n, void *rowptr, void *colind,
                   void *values, void *u, void *v, void *reinit)
{
    unsigned int nr_threads;
    unsigned int *cpus __attribute__ ((unused));
    csx::SPM *spms;
    spm_mt_t *spm_mt;
    elmer_index_t n_ = *((elmer_index_t *) n);
    elmer_index_t *rowptr_ = (elmer_index_t *) rowptr;
    elmer_index_t *colind_ = (elmer_index_t *) colind;
    elmer_value_t *values_ = (elmer_value_t *) values;
    bool reinit_ = *((bool *) reinit);
    elmer_value_t *x_ = (elmer_value_t *) u;
    elmer_value_t *y_ = (elmer_value_t *) v;

    elmer_index_t nr_nzeros = rowptr_[n_] - 1;
    ++nr_calls;
    // char mmf_name[32];
    // snprintf(mmf_name, 32, "test_%ld.mtx", nr_calls);
    // std::ofstream mmf(mmf_name);
    // write_mmf(mmf, rowptr_, colind_, values_, n_);
    // mmf.close();
    if (!*tuned || reinit_) {
        std::cout << "Hello!" << std::endl;
        mt_get_options(&nr_threads, &cpus);
        spms = csx::SPM::LoadCSR_mt<elmer_index_t, elmer_value_t>
            (rowptr_, colind_, values_, n_, n_, false, nr_threads);
        *tuned = spm_mt = GetSpmMt(NULL, spms);
        // std::ofstream mmf_file("test.mtx");
        // write_mmf(mmf_file, rowptr_, colind_, values_, n_);
        // mmf_file.close();
        // *tuned = spm_mt = GetSpmMt("test.mtx");
        csr_orig.rowptr = (elmer_index_t *) malloc(n_*sizeof(elmer_index_t));
        csr_orig.colind =
            (elmer_index_t *) malloc(nr_nzeros*sizeof(elmer_index_t));
        csr_orig.values =
            (elmer_value_t *) malloc(nr_nzeros*sizeof(elmer_value_t));
        memcpy(csr_orig.rowptr, rowptr_, n_*sizeof(elmer_index_t));
        memcpy(csr_orig.colind, colind_, nr_nzeros*sizeof(elmer_index_t));
        memcpy(csr_orig.values, values_, nr_nzeros*sizeof(elmer_value_t));
        csr_orig.nr_rows = n_;
        csr_orig.nr_nzeros = nr_nzeros;
    } else {
        spm_mt = (spm_mt_t *) *tuned;
    }

    assert(csr_orig.nr_rows == n_);
    assert(csr_orig.nr_nzeros == nr_nzeros);

    for (elmer_index_t i = 0; i < n_; ++i) {
        assert(csr_orig.rowptr[i] == rowptr_[i]);
    }

    for (elmer_index_t i = 0; i < nr_nzeros; ++i) {
        assert(csr_orig.colind[i] == colind_[i]);
        assert(csr_orig.values[i] == values_[i]);
    }

    // FIXME: Although Elmer uses doubles by default, this is not portable
    vector_double_t *vec_x = vector_double_create_from_buff(x_, n_);
    vector_double_t *vec_y = vector_double_create_from_buff(y_, n_);
    vector_double_t *y_check = vector_double_create(n_);
    vector_double_init(vec_y, 0);
    vector_double_init(y_check, 0);
    std::cout << "Call No. " << nr_calls << std::endl;
    // std::cout << "Y_before = ";
    // print_array(y_, n_);
    // printf("y_before[22] = %lf\n", y_[22]);
    spmv_double_matvec_mt(spm_mt, vec_x, vec_y);
    matvec_csr(rowptr_, colind_, values_, n_, vec_x->elements,
               y_check->elements);
    // print_prod(rowptr_, colind_, values_, vec_x->elements, 22);
    // std::cout << "Y_after = ";
    // print_array(y_, n_);
    // printf("y_after[22] = %lf\n", y_[22]);
    vector_double_compare(vec_y, y_check);
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
