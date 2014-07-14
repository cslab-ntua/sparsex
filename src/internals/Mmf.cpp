/*
 * \file Mmf.cpp
 *
 * \brief Matrix Market Format utilities
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/internals/Mmf.hpp>

namespace sparsex {
namespace io {

// For testing purposes
void ReadMmfSizeLine(const char *mmf_file, size_t &nr_rows, size_t &nr_cols,
                     size_t &nr_nzeros)
{
    ifstream in;

    try {
        in.open(mmf_file);
        if (!in.is_open()) {
            throw ios_base::failure("");
        }
        in.seekg(0, ios::beg);
        // Ignore comments
        while (in.peek() == '%') {
            in.ignore(numeric_limits<streamsize>::max(), '\n');
        }
    } catch (ios_base::failure &e) {
        LOG_ERROR << "MMF file error\n";
        exit(1);
    }

    vector<string> arguments;
    if (!(DoRead(in, arguments))) {
        LOG_ERROR << "size line error in MMF file\n";
        exit(1);
    }

    ParseElement(arguments, nr_rows, nr_cols, nr_nzeros);
    in.close();
}

// Returns false at EOF
bool DoRead(ifstream &in, vector<string> &arguments)
{
    string buff;

    try {
        if (getline(in, buff).eof()) {
            return false;
        }
    } catch (ios_base::failure &e) {
        LOG_ERROR << "error reading from MMF file: " << (string) e.what();
        exit(1);
    }

    buff = boost::trim_left_copy(buff);
    buff = boost::trim_right_copy(buff);
    boost::split(arguments, buff, boost::algorithm::is_any_of(" "),
        boost::algorithm::token_compress_on);
    
    return true;
}

} // end of namespace io
} // end of namespace sparsex
