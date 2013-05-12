/*
 * mmf.cc -- Matrix Market Format routines
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011,      Vasileios Karakasis
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "mmf.h"

#include <iostream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace csx;
using namespace std;

namespace csx {

// To be removed
void ReadMmfSizeLine(const char *mmf_file, uint64_t &nr_rows, uint64_t &nr_cols,
                     uint64_t &nr_nzeros)
{
    ifstream in;

    in.open(mmf_file);
    in.seekg(0, ios::beg);

    // Ignore comments
    while (in.peek() == '%') {
        in.ignore(numeric_limits<std::streamsize>::max(), '\n');
    }

    vector<string> arguments;
    if (!(DoRead(in, arguments))) {
        cerr << "Size line error" << endl;
        exit(1);
    }

    ParseElement(arguments, nr_rows, nr_cols, nr_nzeros);
    in.close();
}

// Returns false at EOF
bool DoRead(istream &in, vector<string> &arguments)
{
    string buff;

    if (getline(in, buff).eof()) {
        return false;
    }

    boost::split(arguments, buff, boost::algorithm::is_any_of(" "),
        boost::algorithm::token_compress_on);
    
    return true;
}

template <typename IndexType, typename ValueType>
void ParseElement(vector<string> &arguments, IndexType &y, IndexType &x, 
                  ValueType &v)
{
    if (arguments.size() == 3) {
        y = boost::lexical_cast<IndexType,string>(arguments[0]);
        x = boost::lexical_cast<IndexType,string>(arguments[1]);
        v = boost::lexical_cast<ValueType,string>(arguments[2]);
    } else {
        cerr << "Bad input: less arguments in line" << endl;
        exit(1);
    }
}

} //end csx namespace

boost::unordered_map<csx::MMF::MmfInfo, const string> csx::MMF::names_ =
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

MMF::MMF(istream &in)
  : 
    nr_rows_(0),
    nr_cols_(0),
    nr_nzeros_(0),
    in_(in),
    symmetric_(false), 
    col_wise_(false),
    zero_based_(false),
    reordered_(false),
    file_mode_(0)
{
    vector<string> arguments;

    DoRead(in, arguments);
    ParseMmfHeaderLine(arguments);
    ParseMmfSizeLine(arguments);

    if (symmetric_ || col_wise_) {
       DoLoadMmfMatrix();
    }
}

void MMF::ParseMmfHeaderLine(vector<string> &arguments)
{
    // Check if header line exists
    if (arguments[0] != names_[Banner]) {
        if (arguments[0].length() > 2 && arguments[0][0] == '%'&&
            arguments[0][1] == '%') {
            // Header exists but is erroneous so exit
            cerr << "Header line error" << endl;
            exit(1);
        } else {
            // Parse as size line
            file_mode_ = 1;
            return;
        }        
    }

    size_t length;
    if ((length = arguments.size()) < 5) {
        cerr << "Header line error: less arguments" << endl;
        exit(1);
    }

    // Convert to lowercase just in case
    BOOST_FOREACH(string &t, arguments) {
        boost::algorithm::to_lower(t);
    }

    if (arguments[1] != names_[Matrix]) {
        cerr << "Unsupported object" << endl;
        exit(1);
    }

    if (arguments[2] != names_[Coordinate]) {
        cerr << "Unsupported matrix format" << endl;
        exit(1);
    }

    /*if (arguments[3] == names_[Real]) {
        //set(REAL);
    } else if (arguments[3] == names_[Double]) {
        //set(DOUBLE);
    } else if (arguments[3] == names_[Integer]) {
       //set(INTEGER);
    } else {
        cerr << "Unsupported value format" << endl;
        exit(1);
    }*/

    if (arguments[4] == names_[General]) {
        symmetric_ = false;
    } else if (arguments[4] == names_[Symmetric]) {
        symmetric_ = true;
    } else {
        cerr << "Unsupported symmetry" << endl;
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

void MMF::ParseMmfSizeLine(vector<string> &arguments)
{
    bool ignore_comments = false;

    if (file_mode_ && arguments[0][0] == '%') {
        ignore_comments = true;
    }

    if (!file_mode_ || ignore_comments) {
        while (in_.peek() == '%') {
            in_.ignore(numeric_limits<std::streamsize>::max(), '\n');
        }
        if (!DoRead(in_, arguments)) {
            cerr << "Size line error" << endl;
            exit(1);
        }
    }

    ParseElement(arguments, nr_rows_, nr_cols_, nr_nzeros_);
}

void MMF::DoLoadMmfMatrix()
{
    CooElem elem;
    uint64_t tmp;

    if (symmetric_) {
        matrix_.reserve(nr_nzeros_ << 1);
        for (size_t i = 0; i < nr_nzeros_; i++) {
            if (!MMF::GetNext(elem.row, elem.col, elem.val)) {
                cerr << "Requesting dereference, but mmf ended" << endl;
                exit(1);
            }
            matrix_.push_back(elem);
            if (elem.row != elem.col) {
                tmp = elem.row;
                elem.row = elem.col;
                elem.col = tmp;
                matrix_.push_back(elem);          
            }
        }
        nr_nzeros_ = matrix_.size();
    } else {
        matrix_.reserve(nr_nzeros_);
        for (size_t i = 0; i < nr_nzeros_; i++) {
            if (!MMF::GetNext(elem.row, elem.col, elem.val)) {
                cerr << "Requesting dereference, but mmf ended" << endl;
                exit(1);
            }
            matrix_.push_back(elem);
        }
    }

    sort(matrix_.begin(), matrix_.end(), CooElemSorter());
}

bool MMF::GetNext(uint64_t &y, uint64_t &x, double &v)
{
    vector<string> arguments;

    if (!DoRead(in_, arguments)) {
        return false;
    }

    ParseElement(arguments, y, x, v);

    if (zero_based_) {
        y++;
        x++;
    }

    return true;
}

void MMF::GetCoordinates(size_t index, uint64_t &row, uint64_t &col)
{
    row = matrix_[index].row;
    col = matrix_[index].col;
}

void MMF::SetCoordinates(size_t index, uint64_t row, uint64_t col)
{
    matrix_[index].row = row;
    matrix_[index].col = col;
}

void MMF::Sort()
{
    sort(matrix_.begin(), matrix_.end(), CooElemSorter());
}

void MMF::Print()
{
    cout << "Elements of Matrix" << endl;
    cout << "------------------" << endl;
    for (size_t i = 0; i < nr_nzeros_; i++) {
        cout << matrix_[i].row << " " << matrix_[i].col << " " << matrix_[i].val << endl;
    }
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
