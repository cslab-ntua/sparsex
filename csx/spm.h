/* -*- C++ -*-
 *
 * spm.h --  Internal representation of sparse matrices.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_SPM_H__
#define CSX_SPM_H__

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/function.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include "spm_bits.h"
#include "mmf.h"
#include "csr.h"
#include "dynarray.h"

namespace csx {

class MMF;

/**
 *  Internal representation of a sparse matrix. This class can also represent a
 *  sub-matrix of a larger one, starting at a specific row of the original
 *  matrix.
 */
class SPM
{
private:
    uint64_t nr_rows_, nr_cols_, nr_nzeros_;
    SpmIterOrder type_;
    SpmRowElem *elems_;
    uint64_t elems_size_;
    uint64_t *rowptr_;
    uint64_t rowptr_size_;
    uint64_t row_start_;    /* Row of the original matrix, where this
                               sub-matrix starts. */
    bool elems_mapped_;
    
    // Maximum possible rowptr_size after transformations.
    uint64_t max_rowptr_size_;

    // These are mainly SPM construction classes, so make them friends
    friend class Builder;
    friend class VirtualBuilder;
    friend class CsxManager;
    friend class DRLE_Manager;

public:
    SPM() : type_(NONE), elems_(NULL), rowptr_(NULL), elems_mapped_(false) {}

    ~SPM()
    {
        if (!elems_mapped_ && elems_)
            free(elems_);

        if (rowptr_)
            free(rowptr_);
    };

    uint64_t GetNrRows()
    {
        return rowptr_size_ - 1;
    }

    uint64_t GetNrCols()
    {
        return nr_cols_;
    }

    uint64_t GetNrNonzeros()
    {
        return nr_nzeros_;
    }

    SpmIterOrder GetType()
    {
        return type_;
    }
    
    uint64_t GetElemsSize()
    {
        return elems_size_;
    }
    
    SpmRowElem *GetElems()
    {
        return elems_;
    }
    
    uint64_t *GetRowPtr()
    {
        return rowptr_;
    }
    
    uint64_t GetRowStart()
    {
        return row_start_;
    }
    
    uint64_t GetMaxRowPtrSize()
    {
        return max_rowptr_size_;
    }

    void SetNrRows(uint64_t nr_rows)
    {
        nr_rows_ = nr_rows;
    }

    void SetNrCols(uint64_t nr_cols)
    {
        nr_cols_ = nr_cols;
    }

    void SetNrNonzeros(uint64_t nr_nzeros)
    {
        nr_nzeros_ = nr_nzeros;
    }

    void SetType(SpmIterOrder type)
    {
        type_ = type;
    }
    
    void SetRowStart(uint64_t row_start)
    {
        row_start_ = row_start;
    }
    
    void SetRowPtrSize(uint64_t rowptr_size)
    {
        rowptr_size_ = rowptr_size;
    }
    
    void SetMaxRowPtrSize(uint64_t max_rowptr_size)
    {
        max_rowptr_size_ = max_rowptr_size;
    }
    
    SpmRowElem *RowBegin(uint64_t ridx = 0);
    SpmRowElem *RowEnd(uint64_t ridx = 0);

    class Builder;
    class VirtualBuilder;

    /**
     *  Function for filling the matrix using point iterators
     *
     *  @param pnts_start point iterators start.
     *  @param pnts_end   point iterators end.
     *  @param first_row  first row of the matrix.
     *  @param limit      limit for the number of elements.
     *  @param nr_elems   size for elems allocation.
     *  @param nr_rows    size for rows allocation.
     *  @param SpmBld     point to a builder responsible for allocations.
     *  @return           number of elements allocated.
     */
    template <typename IterT>
    uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end,
                      uint64_t first_row, unsigned long limit = 0,
                      uint64_t nr_elems = 0, uint64_t nr_rows = 0);
    template <typename IterT>
    uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end,
                      uint64_t first_row, unsigned long limit,
                      uint64_t nr_elems, uint64_t nr_rows,
                      SPM::Builder *SpmBld);

    /**
     *  Loads a matrix from file.
     *
     *  @param filename name of the file that matrix is kept.
     *  @param in       buffer from which the matrix is taken.
     *  @return         spm class object with the characteristics of the matrix.
     */
    static SPM *LoadMMF(const char *filename);
    static SPM *LoadMMF(std::istream &in=std::cin);

    /**
     *  Loads matrix from a file specifying the number of threads.
     *
     *  @param filename name of the file that matrix is kept.
     *  @param in       buffer from which the matrix is taken.
     *  @param mmf      handler of MMF class.
     *  @param nr       number of threads to be used.
     *  @return         spm class object with the characteristics of the matrix.
     */
    static SPM *LoadMMF_mt(const char *mmf_file, const long nr);
    static SPM *LoadMMF_mt(std::istream &in, const long nr);
    static SPM *LoadMMF_mt(MMF &mmf, const long nr);
    
    /**
     *  Loads matrix from CSR format.
     *
     *  @param rowptr     array "rowptr" of CSR format.
     *  @param colind     array "colind" of CSR format.
     *  @param values     array "values" of CSR format.
     *  @param nr_rows    number of rows.
     *  @param nr_cols    number of columns.
     *  @param zero_based """bkk"""
     *  @param nr         number of threads to be used.
     *  @return           spm class object with the characteristics of the 
     *                    matrix.
     */
    template<typename IndexType, typename ValueType>
    static SPM *LoadCSR_mt(IndexType *rowptr, IndexType *colind,
                           ValueType *values, IndexType nr_rows,
                           IndexType nr_cols, bool zero_based, long nr)
    {
        CSR<IndexType, ValueType> csr(rowptr, colind, values, nr_rows,
                                      nr_cols, zero_based);
        return DoLoadMatrix<CSR<IndexType, ValueType> >(csr, nr);
    };

    /**
     *  Printing methods.
     */
    void Print(std::ostream &out = std::cout);
    void PrintElems(std::ostream &out = std::cout);
    void PrintRows(std::ostream &out = std::cout);
    void PrintStats(std::ostream &out = std::cout);

    class PntIter;

    PntIter PointsBegin(uint64_t ridx = 0);
    PntIter PointsEnd(uint64_t ridx = 0);

    /**
     *  Transforms matrix to specific type.
     *
     *  @param type type of pattern to transform to.
     *  @param rs   line to start the transformation. If it is zero, it is the
     *              the beggining of matrix.
     *  @param re   line to end the transformation. If it is zero, it is the
     *              the end of matrix.
     */
    void Transform(SpmIterOrder type, uint64_t rs = 0, uint64_t re = 0);

    /**
     *  Finds the reverse transform function. Transforms from <tt>type</tt> to
     *  <tt>HORIZONTAL</tt>.
     *
     *  @param type type of pattern to transform from.
     *  @return     reverse transform function.
     */
    TransformFn GetRevXformFn(SpmIterOrder type);

    /**
     *  Finds the transform function. Transforms from HORIZONTAL to
     *  <tt>type</tt>.
     *
     *  @param type type of pattern to transform to.
     *  @return     transform function.
     */
    TransformFn GetXformFn(SpmIterOrder type);

    /**
     *  Finds the complete transform function. Tranforms from <tt>from</tt> to
     *  <tt>to</tt>.
     *
     *  @param from type of pattern to transform from.
     *  @param to   type of pattern to transform to.
     *  @return     complete transform function.
     */
    TransformFn GetTransformFn(SpmIterOrder from, SpmIterOrder to);

    /**
     *  Extracts a window (sub-matrix) from this matrix. This methods copies the
     *  nonzero elements of the original matrix to the window.
     *
     *  @param rs     the row in this matrix, from where the window will be
     *                extracted.
     *  @param length number of rows to include in the window.
     *  @return       the selected window.
     *  @see GetWindow()
     */
    SPM *ExtractWindow(uint64_t rs, uint64_t length);

    /**
     *  Extracts a window (sub-matrix) from this matrix using a backing
     *  store. This method uses the nonzero elements of this matrix without
     *  performing any copy.
     *
     *  @param rs     the row in this matrix, from where the window will be
     *                extracted.
     *  @param length number of rows to include in the window.
     *  @return       the selected window.
     *  @see ExtractWindow()
     */
    SPM *GetWindow(uint64_t rs, uint64_t length);

    /**
     *  Puts a (possibly processed) window of the matrix back to this matrix.
     *
     *  @param window the (possibly processed) window.
     */
    void PutWindow(const SPM *window);

/**
 *  """bkk"""
 */
private:
    template<typename MatrixType>
    static SPM *DoLoadMatrix(MatrixType &mat, long nr)
    {
        SPM *ret, *spm;
        long limit, cnt, row_start;
        typename MatrixType::iterator iter = mat.begin();
        typename MatrixType::iterator iter_end = mat.end();

        ret = new SPM[nr];
        row_start = limit = cnt = 0;
        for (long i = 0; i < nr; ++i) {
            spm = ret + i;
            limit = (mat.GetNrNonzeros() - cnt) / (nr - i);
            spm->nr_nzeros_ = spm->SetElems(iter, iter_end, row_start + 1,
                                            limit,
                                            limit + mat.GetNrRows() + 1,
                                            mat.GetNrRows() + 1);
            spm->nr_rows_ = spm->rowptr_size_ - 1;
            spm->nr_cols_ = mat.GetNrCols();
            spm->row_start_ = row_start;
            row_start += spm->nr_rows_;
            spm->type_ = HORIZONTAL;
            spm->max_rowptr_size_ = spm->nr_rows_ + spm->nr_cols_ + 1;
            cnt += spm->nr_nzeros_;
        }

        assert((uint64_t) cnt == mat.GetNrNonzeros());
        return ret;
    }

    void print_prod(uint64_t row_no)
    {
        for (uint64_t i = rowptr_[row_no]; i < rowptr_[row_no+1]; ++i) {
            std::cout << "i = " << i << ", "
                      << "elem = " << elems_[i] << "|";
        }
        std::cout << std::endl;
    }
};


void TestMMF(SPM *spm, const char *mmf_file);

/**
 *  Utility class for building an SPM.
 */
class SPM::Builder
{
public:
    /**
     *  Constructs a Builder object.
     *
     *  @param spm      the matrix to be built.
     *  @param nr_elems an initial guess of the nonzero elements.
     *  @param nr_rows  an initial guess of the number of rows.
     */
    Builder(SPM *spm, uint64_t nr_elems = 0, uint64_t nr_rows = 0);
    virtual ~Builder();

    /**
     *  Allocates a new element in the new matrix.
     *
     *  @return pointer to element allocated.
     */
    virtual SpmRowElem *AllocElem();

    /**
     *  Allocates a number of elements in the new matrix.
     *
     *  @param  nr number of elements to be allocated.
     *  @return    pointer to the first allocated element.
     */
    virtual SpmRowElem *AllocElems(uint64_t nr);

    /**
     *  Counts elements already allocated in the matrix.
     *
     *  @return number of elements allocated.
     */
    virtual uint64_t GetElemsCnt();

    /**
     *  Allocate elements for the next row.
     *
     *  @param rdiff number of empty rows between the current and the new row.
     */
    virtual void NewRow(uint64_t rdiff = 1);

    /**
     *  Finalizes the construction of the new matrix.
     */
    virtual void Finalize();

private:
    SPM *spm_;
    dynarray_t *da_elems_;
    dynarray_t *da_rowptr_;
};

/**
 *  Iterators for moving between elements inside the matrix.
 */
class SPM::PntIter : public std::iterator<std::forward_iterator_tag, CooElem>
{
public:
    PntIter(): spm_(NULL), row_idx_(0), elem_idx_(0) { }
    PntIter(SPM *spm, uint64_t row_idx);

    bool operator==(const PntIter &pi);
    bool operator!=(const PntIter &pi);
    void operator++();
    SpmCooElem operator*();
    std::ostream &operator<<(std::ostream &out);

private:
    SPM *spm_;
    uint64_t row_idx_;
    uint64_t elem_idx_;
};


/**
 *  Internal representation of a sparse symmetric matrix. This class can also 
 *  represent a sub-matrix of a larger one.
 */
class SPMSym
{
private:
    SPM         *lower_matrix_;   // Representation of lower matrix.
    double      *diagonal_;       // Values of diagonal elements.
    uint64_t    diagonal_size_;   // Size of the diagonal.
    SPM         *m1_;             /* Matrix that contains the elems of SpmSym
                                     for which column of the element is smaller
                                     than the row start. */
    SPM         *m2_;             /* Matrix that contains the rest of the
                                     elements. */
    
public:
    SPMSym() : lower_matrix_(NULL), diagonal_(NULL), m1_(NULL), m2_(NULL)
    {
        lower_matrix_ = new SPM;
        m1_ = new SPM;
        m2_ = new SPM;
    }

    ~SPMSym()
    {
        if (lower_matrix_)
            delete lower_matrix_;
        
        if (m1_)
            delete m1_;
        
        if (m2_)
            delete m2_;
        
        if (diagonal_)
            free(diagonal_);
    };
    
    SPM *GetLowerMatrix()
    {
        return lower_matrix_;
    }
    
    SPM *GetFirstMatrix()
    {
        return m1_;
    }
    
    SPM *GetSecondMatrix()
    {
        return m2_;
    }
    
    double *GetDiagonal()
    {
        return diagonal_;
    }
    
    uint64_t GetDiagonalSize()
    {
        return diagonal_size_;
    }
    
    void SetDiagonalSize(uint64_t size)
    {
        diagonal_size_ = size;
    }
    
    void SetDiagonal(double *diagonal)
    {
        diagonal_ = diagonal;
    }
    
    /**
     *  Loads matrix from a file specifying the number of threads.
     *
     *  @param filename name of the file that matrix is kept.
     *  @param in       buffer from which the matrix is taken.
     *  @param mmf      handler of MMF class.
     *  @param nr       number of threads to be used.
     *  @return         spmsym class object with the characteristics of the 
     *                  matrix.
     */
    static SPMSym *LoadMMF_mt(const char *mmf_file, const long nr);
    static SPMSym *LoadMMF_mt(std::istream &in, const long nr);
    static SPMSym *LoadMMF_mt(MMF &mmf, const long nr);
    
    /**
     *  Divide the matrix (m) to two submatrices (m1 and m2).
     */
    void DivideMatrix();

    /**
     *  Merge the two submatrices (m1 and m2) into one matrix (m).
     */
    void MergeMatrix();
    
    class Builder;
    
    /**
     *  Function for filling the matrix using point iterators
     *
     *  @param pnts_start   point iterators start.
     *  @param pnts_end     point iterators end.
     *  @param first_row    first row of the matrix.
     *  @param limit        limit for the number of elements.
     *  @param nr_elems     size for elems allocation.
     *  @param nr_rows      size for rows allocation.
     *  @param spm_sym_bld  point to a builder responsible for allocations.
     *  @return             number of elements allocated.
     */
    template <typename IterT>
    uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end,
                      uint64_t first_row, unsigned long limit = 0,
                      uint64_t nr_elems = 0, uint64_t nr_rows = 0);
    template <typename IterT>
    uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end,
                      uint64_t first_row, unsigned long limit,
                      uint64_t nr_elems, uint64_t nr_rows,
                      SPMSym::Builder *spm_sym_bld);
                      
    void PrintDiagElems(std::ostream &out = std::cout);
};

class SPMSym::Builder
{
public:
    /**
     *  Constructs a Builder object.
     *
     *  @param spm_sym  the matrix to be built.
     *  @param nr_elems an initial guess of the nonzero elements.
     *  @param nr_rows  an initial guess of the number of rows.
     */
    Builder(SPMSym *spm_sym, uint64_t nr_elems = 0, uint64_t nr_rows = 0);
    virtual ~Builder();
    
    SPM::Builder *GetLowerBuilder()
    {
        return spm_bld_;
    }
    
    double *AllocDiagElem();
    uint64_t GetDiagElemsCnt();
    uint64_t GetElemsCnt();
    void Finalize();
    
private:
    SPMSym *spm_sym_;
    SPM::Builder *spm_bld_;
    dynarray_t *da_dvalues_;
};

/*
 *  SetElems functions
 */ 
template<typename IterT>
uint64_t SPM::SetElems(IterT &pi, const IterT &pnts_end, uint64_t first_row,
                       unsigned long limit, uint64_t nr_elems, uint64_t nrows,
                       SPM::Builder *SpmBld)
{
    SpmRowElem *elem;
    uint64_t row_prev, row;

    row_prev = first_row;
    for (; pi != pnts_end; ++pi) {
        row = (*pi).y;
        if (row != row_prev) {
            assert(row > row_prev);
            if (limit && SpmBld->GetElemsCnt() >= limit)
                break;
            SpmBld->NewRow(row - row_prev);
            row_prev = row;
        }

        elem = SpmBld->AllocElem();
        MakeRowElem(*pi, elem);
    }

    return elems_size_;
}

template<typename IterT>
uint64_t SPM::SetElems(IterT &pi, const IterT &pnts_end, uint64_t first_row,
                       unsigned long limit, uint64_t nr_elems, uint64_t nrows)
{
    SPM::Builder *SpmBld = new SPM::Builder(this, nr_elems, nrows);

    SetElems(pi, pnts_end, first_row, limit, nr_elems, nrows, SpmBld);
    SpmBld->Finalize();
    delete SpmBld;
    return elems_size_;
}

template<typename IterT>
uint64_t SPMSym::SetElems(IterT &pi, const IterT &pnts_end, uint64_t first_row,
                          unsigned long limit, uint64_t nr_elems,
                          uint64_t nrows, SPMSym::Builder *spm_sym_bld)
{
    double *value;
    SpmRowElem *elem;
    uint64_t row_prev, row, col;

    row_prev = first_row;
    for (; pi != pnts_end; ++pi) {
        row = (*pi).y;
        col = (*pi).x;
        
        if (row > col) {
            if (row != row_prev) {
                assert(row > row_prev);
                if (limit && spm_sym_bld->GetElemsCnt() >= limit && 
                    row_prev == row - 1)
                    break;
                spm_sym_bld->GetLowerBuilder()->NewRow(row - row_prev);
                row_prev = row;
            }
            elem = spm_sym_bld->GetLowerBuilder()->AllocElem();
            MakeRowElem(*pi, elem);
        } else if (row == col) {
            value = spm_sym_bld->AllocDiagElem();
            *value = (*pi).val;
        }
    }
    return lower_matrix_->GetElemsSize() + GetDiagonalSize();
}

template<typename IterT>
uint64_t SPMSym::SetElems(IterT &pi, const IterT &pnts_end, uint64_t first_row,
                          unsigned long limit, uint64_t nr_elems,
                          uint64_t nrows)
{
    SPMSym::Builder *spm_sym_bld = new SPMSym::Builder(this, nr_elems, nrows);

    SetElems(pi, pnts_end, first_row, limit, nr_elems, nrows, spm_sym_bld);
    spm_sym_bld->Finalize();
    delete spm_sym_bld;
    return lower_matrix_->GetElemsSize() + GetDiagonalSize();
}

}   // end of csx namespace

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
