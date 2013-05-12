/* -*- C++ -*-
 *
 * mmf.h -- Matrix Market Format routines
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011,      Vasileios Karakasis
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef MMF_H__
#define MMF_H__

#include "spm_bits.h"

#include <algorithm>
#include <boost/bind/bind.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

namespace csx {

// To be removed
void ReadMmfSizeLine(const char *mmf_file, uint64_t &nr_rows, uint64_t &nr_cols,
                     uint64_t &nr_nzeros);

bool DoRead(istream &in, vector<string> &arguments);
template <typename IndexType, typename ValueType>
void ParseElement(vector<string> &arguments, IndexType &y,IndexType &x,
                  ValueType &v);

class MMF
{
public:

    MMF(istream &in);
    
    size_t GetNrRows() const 
    {
        return nr_rows_;
    }

    size_t GetNrCols() const
    {
        return nr_cols_;
    }

    size_t GetNrNonzeros() const
    {
        if (symmetric_ || col_wise_)
            return matrix_.size();
        else 
            return nr_nzeros_;
    }

    bool IsSymmetric() const
    {
        return symmetric_;
    }

    bool IsColWise() const
    {
        return col_wise_;
    }
    
    bool IsZeroBased() const
    {
        return zero_based_;
    }

    void GetCoordinates(size_t index, uint64_t &row, uint64_t &col);
    void SetCoordinates(size_t index, uint64_t row, uint64_t col);
    void Sort();
    void Print();

    void InitMatrix(size_t size)
    {
        matrix_.reserve(size);
    }
    
    void InsertElement(CooElem elem)
    {
        matrix_.push_back(elem);
    }

    void SetReordered()
    {
        reordered_ = true;
    }

    class iterator;
    iterator begin();
    iterator end();

private:
    size_t nr_rows_, nr_cols_, nr_nzeros_;
    istream &in_;
    bool symmetric_, col_wise_, zero_based_, reordered_;
    int file_mode_;     // 0 for MMF files, 1 for regular files
    vector<CooElem> matrix_;

    enum MmfInfo {
        Banner,
        Matrix,
        Coordinate,
        Real,
        Double,
        Integer,
        General,
        Symmetric,
        Indexing0,
        Indexing1,
        ColumnWise,
        RowWise
    };

    static boost::unordered_map<MmfInfo, const string> names_;

    void ParseMmfHeaderLine(vector<string> &arguments);
    void ParseMmfSizeLine(vector<string> &arguments); 
    void DoLoadMmfMatrix();
    bool GetNext(uint64_t &y, uint64_t &x, double &val);
};

class MMF::iterator : public std::iterator<forward_iterator_tag, CooElem>
{
public:
    iterator() {}
    
    iterator(MMF *mmf, uint64_t cnt)
      :
        mmf_(mmf),
        cnt_(cnt)
    {
        if (mmf_->symmetric_ || mmf_->col_wise_ || mmf_->reordered_)
            return;

        // this is the initializer
        if (cnt_ == 0) {
            this->DoSet();
        }
    }
    
    bool operator==(const iterator &i)
    {
        //cout << "me: " << mmf_ << " " << cnt_
        //     << " i: " << i.mmf << " " << i.cnt << "\n";
        return (mmf_ == i.mmf_) && (cnt_ == i.cnt_);
    }

    bool operator!=(const iterator &i)
    {
        return !(*this == i);
    }

    void operator++()
    {
        ++cnt_;
        if (!mmf_->symmetric_ || !mmf_->col_wise_ || !mmf_->reordered_) {
            this->DoSet();
        }
    }

    CooElem operator*()
    {
        if (mmf_->symmetric_ || mmf_->col_wise_ || mmf_->reordered_) {
            return mmf_->matrix_[cnt_];
        } else {
            if (!valid_) {
                cerr << "Requesting dereference, but mmf ended\n"
                          << "cnt: " << cnt_ << endl;
                exit(1);
            }
            assert(valid_);
            return elem_;
        }
    }

private:
    void DoSet()
    {
        valid_ = mmf_->GetNext(elem_.row, elem_.col, elem_.val);
    }

    MMF *mmf_;
    uint64_t cnt_;
    CooElem elem_;
    bool valid_;
};

MMF::iterator MMF::begin()
{
    return MMF::iterator(this, 0);
}

MMF::iterator MMF::end()
{
    if (this->symmetric_ || this->col_wise_ || this->reordered_) {
        return MMF::iterator(this, this->matrix_.size());
    } else {
        return MMF::iterator(this, this->nr_nzeros_);
    }
}

}  //csx namespace end

#endif  // MMF_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
