/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file BenchUtil.hpp
 * \brief Benchmarking utilities
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef BENCH_UTIL_HPP
#define BENCH_UTIL_HPP

#if SPX_BENCH_POSKI
#   include "PoskiModule.hpp"
#endif
#if SPX_BENCH_MKL
#   include "MklModule.hpp"
#endif
#include "SparsexModule.hpp"

#include <cstdlib>
#include <iostream>
#include <dirent.h>
#include <boost/function.hpp>

using namespace std;

namespace bench {

/* The available libraries */
enum library {
    SparseX,
    MKL,
    pOSKI
};

typedef boost::function<void (spx_index_t*, spx_index_t*, spx_value_t*,
                              spx_index_t, spx_index_t, spx_index_t,
                              spx_value_t*, spx_value_t*,
                              spx_value_t, spx_value_t)> SpmvFn;

class DirectoryIterator
{
public:
    DirectoryIterator(const char *directory) 
        : directory_(directory),
          directory_id_(0),
          file_(),
          is_valid_(0)
    {
        directory_id_ = opendir(directory);
        if (!directory_id_) {
            cerr << "[BENCH]: failed to open directory." << endl;
            exit(1);
        } else {
            is_valid_ = true;
        }
        ++(*this);
    }

    ~DirectoryIterator()
    {
        if (directory_id_)
            closedir(directory_id_);
    }

    void operator++()
    {
        struct dirent *cur_entry = 0;
        while (((cur_entry = readdir(directory_id_)) != 0)
               && dot_or_dot_dot(cur_entry->d_name));
        if (cur_entry != 0) {
            if (cur_entry->d_type != DT_DIR) {
                string matrix_name_ = static_cast<string>(cur_entry->d_name);
                file_ = directory_ + "/" + matrix_name_;
                is_valid_ = true;
            } else {
                is_valid_ = false;
            }
        } else {
            is_valid_ = false;
        }
    }

    operator bool() const 
    {
        return is_valid_;
    }

    string &filename()
    {
        return file_;
    }

private:
    const string directory_;
    DIR *directory_id_;
    string file_;
    bool is_valid_;

    DirectoryIterator() {}

    inline bool dot_or_dot_dot(const char *name)
    {
        return name[0] == '.' 
            && (name[1] == '\0' || (name[1] == '.' && name[2] == '\0'));
    }
};

} // end of namespace bench

#endif  // BENCH_UTIL_HPP
