/*
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011,      Vasileios Karakasis
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Mmf.hpp
 * \brief Matrix Market Format utilities
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_MMF_HPP
#define SPARSEX_INTERNALS_MMF_HPP

#include <sparsex/internals/Element.hpp>
#include <sparsex/internals/logger/Logger.hpp>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/bind/bind.hpp>
#include <boost/container/vector.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace sparsex::csx;

namespace sparsex {
namespace io {

// For testing purposes; to be removed
void ReadMmfSizeLine(const char *mmf_file, size_t& nr_rows, size_t& nr_cols,
                     size_t& nr_nzeros);

// Returns false at EOF
bool DoRead(ifstream &in, vector<string> &arguments);

template<typename IndexType, typename ValueType>
void ParseElement(vector<string> &arguments, IndexType &y, IndexType &x, 
                  ValueType &v);

template<typename IndexType, typename ValueType>
class MMF
{
public:
    typedef IndexType idx_t;
    typedef ValueType val_t;

    MMF(const char *filename);
    MMF(IndexType *rowptr, IndexType *colind, ValueType *values,
        size_t nr_rows, size_t nr_cols, bool zero_based) {} //Dummy
    ~MMF()
    {
        if (in.is_open())
            in.close();
    }

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

    bool IsReordered() const
    {
        return reordered_;
    }

    pair<IndexType, IndexType> GetCoordinates(IndexType idx) const
    {
        return matrix_[idx].GetCoordinates();
    }

    void SetCoordinates(IndexType idx, IndexType row, IndexType col)
    {
        matrix_[idx] = Element<IndexType, ValueType>(row, col,
                                                     matrix_[idx].GetValue());
    }

    void Sort()
    {
        sort(matrix_.begin(), matrix_.end());
    }

    /*void Print(std::ostream &os) const
    {
        os << "Elements of Matrix" << endl;
        os << "------------------" << endl;
        if (symmetric_ || col_wise_) {
            for (size_t i = 0; i < matrix_.size(); i++) {
                os << matrix_[i].row << " " << matrix_[i].col << " "
                   << matrix_[i].val << endl;
            }
        }//  else {
        //     iterator iter = begin();
        //     iterator iter_end = end();
        //     for (;iter != iter_end; ++iter) {
        //         cout << (*iter).row << " " << (*iter).col << " "
        //              << (*iter).val << endl;
        //     }
        // }
        os << endl;
        }*/

    void InitMatrix(size_t size)
    {
        matrix_.reserve(size);
    }

    void ResetStream();
    
    void InsertElement(const Element<IndexType, ValueType> &elem)
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
    ifstream in;
    bool symmetric_, col_wise_, zero_based_, reordered_;
    int file_mode_;     // 0 for MMF files, 1 for regular files
    boost::container::vector<Element<IndexType, ValueType> > matrix_;

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
    bool GetNext(IndexType &r, IndexType &c, ValueType &v);
};

template<typename IndexType, typename ValueType>
class MMF<IndexType, ValueType>::iterator 
    : public std::iterator<forward_iterator_tag, Element<IndexType, ValueType> >
{
public:
    // In case subsequent calls to DoLoadMatrix() are made
    ~iterator() {
        mmf_->ResetStream();
    }
    
    iterator(MMF *mmf, size_t cnt)
      : mmf_(mmf),
        cnt_(cnt),
        elem_(Element<IndexType, ValueType>(0, 0, 0)),
        valid_(false)
    {
        if (mmf_->symmetric_ || mmf_->col_wise_ || mmf_->reordered_)
            return;

        // this is the initializer
        if (cnt_ == 0)
            DoSet();
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
        if (mmf_->symmetric_ || mmf_->col_wise_ || mmf_->reordered_) {
            return;
        }

        DoSet();
    }

    Element<IndexType, ValueType> &operator*()
    {
        if (mmf_->symmetric_ || mmf_->col_wise_ || mmf_->reordered_) {
            return mmf_->matrix_[cnt_];
        } else {
            if (!valid_) {
                LOG_ERROR << "Requesting dereference, but mmf ended (cnt: " 
                          << cnt_ << "/" << mmf_->nr_nzeros_ << ").\n";
                exit(1);
            }

            return elem_;
        }
    }

private:
    void DoSet()
    {
        IndexType row, col;
        ValueType val;
        valid_ = mmf_->GetNext(row, col, val);
        elem_ = Element<IndexType, ValueType>(row, col, val);
    }

    MMF *mmf_;
    size_t cnt_;
    Element<IndexType, ValueType> elem_;
    bool valid_;
};

template<typename IndexType, typename ValueType>
typename MMF<IndexType, ValueType>::iterator MMF<IndexType, ValueType>::begin()
{
    return iterator(this, 0);
}

template<typename IndexType, typename ValueType>
typename MMF<IndexType, ValueType>::iterator MMF<IndexType, ValueType>::end()
{
    if (this->symmetric_ || this->col_wise_ || this->reordered_) {
        return iterator(this, matrix_.size());
    } else {
        return iterator(this, nr_nzeros_);
    }
}

/*
 * Implementation of class MMF
 */
template<typename IndexType, typename ValueType>
boost::unordered_map<typename MMF<IndexType, ValueType>::MmfInfo, const string>
MMF<IndexType, ValueType>::names_ =
                     boost::assign::map_list_of 
                     (Banner, "%%MatrixMarket")
                     (Matrix, "matrix")
                     (Coordinate, "coordinate")
                     (Real, "real")
                     (Double, "double")
                     (Integer, "integer") 
                     (General, "general")
                     (Symmetric, "symmetric")
                     (Indexing0, "0-base")
                     (Indexing1, "1-base")
                     (ColumnWise, "column")
                     (RowWise, "row");

template<typename IndexType, typename ValueType>
MMF<IndexType, ValueType>::MMF(const char* filename)
  : 
    nr_rows_(0),
    nr_cols_(0),
    nr_nzeros_(0),
    symmetric_(false), 
    col_wise_(false),
    zero_based_(false),
    reordered_(false),
    file_mode_(0)
{
    try {
        in.open(filename);
        if (!in.is_open()) {
            throw ios_base::failure("");
        }
    } catch (ios_base::failure &e) {
        LOG_ERROR << "MMF file error\n";
        exit(1);
    }
    vector<string> arguments;

    DoRead(in, arguments);
    ParseMmfHeaderLine(arguments);
    ParseMmfSizeLine(arguments);

    if (symmetric_ || col_wise_) {
        DoLoadMmfMatrix();
    }
}

template<typename IndexType, typename ValueType>
void MMF<IndexType, ValueType>::ParseMmfHeaderLine(vector<string> &arguments)
{
    // Check if header line exists
    if (arguments[0] != names_[Banner]) {
        if (arguments[0].length() > 2 && arguments[0][0] == '%'&&
            arguments[0][1] == '%') {
            // Header exists but is erroneous so exit
            LOG_ERROR << "invalid header line in MMF file\n";
            exit(1);
        } else {
            // Parse as size line
            file_mode_ = 1;
            col_wise_ = false;
            return;
        }        
    }

    size_t length;
    if ((length = arguments.size()) < 5) {
        LOG_ERROR << "less arguments in header line of MMF file\n";
        exit(1);
    }

    // Convert to lowercase just in case
    BOOST_FOREACH(string &t, arguments) {
        boost::algorithm::to_lower(t);
    }

    if (arguments[1] != names_[Matrix]) {
        LOG_ERROR << "unsupported object in header line of MMF file\n";
        exit(1);
    }

    if (arguments[2] != names_[Coordinate]) {
        LOG_ERROR << "unsupported matrix format in header line of MMF file\n";
        exit(1);
    }

    if (arguments[4] == names_[General]) {
        symmetric_ = false;
    } else if (arguments[4] == names_[Symmetric]) {
        symmetric_ = true;
    } else {
        LOG_ERROR << "unsupported symmetry in header line of MMF file\n";
        exit(1);
    }
    
    if (length > 5) {
        for (size_t i = 5; i < length; i++) {
            if (arguments[i] == names_[Indexing0]) zero_based_ = true;
            else if (arguments[i] == names_[Indexing1]) zero_based_ = false;
            else if (arguments[i] == names_[ColumnWise]) col_wise_ = true;
            else if (arguments[i] == names_[RowWise]) col_wise_ = false;
        }
    }
}

template<typename IndexType, typename ValueType>
void MMF<IndexType, ValueType>::ParseMmfSizeLine(vector<string> &arguments)
{
    bool ignore_comments = false;

    if (file_mode_ && arguments[0][0] == '%') {
        ignore_comments = true;
    }

    if (!file_mode_ || ignore_comments) {
        while (in.peek() == '%') {
            in.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        if (!DoRead(in, arguments)) {
            LOG_ERROR << "size line error in MMF file\n";
            exit(1);
        }
    }
    ParseElement(arguments, nr_rows_, nr_cols_, nr_nzeros_);
}

template<typename IndexType, typename ValueType>
void MMF<IndexType, ValueType>::DoLoadMmfMatrix()
{
    IndexType row, col;
    ValueType val;

    if (symmetric_) {
        matrix_.reserve(nr_nzeros_ << 1);
        for (size_t i = 0; i < nr_nzeros_; i++) {
            if (!MMF::GetNext(row, col, val)) {
                LOG_ERROR << "Requesting dereference, but mmf ended.\n";
                exit(1);
            }

            matrix_.emplace_back(row, col, val);
            if (row != col)
                // Emplace also the symmetric element
                matrix_.emplace_back(col, row, val);
        }
    } else {
        matrix_.reserve(nr_nzeros_);
        for (size_t i = 0; i < nr_nzeros_; i++) {
            if (!MMF::GetNext(row, col, val)) {
                LOG_ERROR << "Requesting dereference, but mmf ended.\n";
                exit(1);
            }

            matrix_.emplace_back(row, col, val);
        }
    }

    sort(matrix_.begin(), matrix_.end());
}

template<typename IndexType, typename ValueType>
bool MMF<IndexType, ValueType>::GetNext(IndexType &r, IndexType &c,
                                        ValueType &v)
{
    vector<string> arguments;

    if (!DoRead(in, arguments)) {
        return false;
    }

    ParseElement(arguments, r, c, v);
    if (zero_based_) {
        r++;
        c++;
    }

    return true;
}

template<typename IndexType, typename ValueType>
void MMF<IndexType, ValueType>::ResetStream()
{
    in.clear();
    in.seekg(0, ios::beg);
    while (in.peek() == '%') {
        in.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    in.ignore(numeric_limits<streamsize>::max(), '\n');
}

template<typename IndexType, typename ValueType>
void ParseElement(vector<string> &arguments, IndexType &y, IndexType &x, 
                  ValueType &v)
{
    if (arguments.size() == 3) {
        y = boost::lexical_cast<IndexType,string>(arguments[0]);
        x = boost::lexical_cast<IndexType,string>(arguments[1]);
        v = boost::lexical_cast<ValueType,string>(arguments[2]);
    } else {
        LOG_ERROR << "bad input, less arguments in line of MMF file\n";
        exit(1);
    }
}

} // end of namespace io
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_MMF_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
