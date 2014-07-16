/*
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file JitUtil.hpp
 * \brief Jit-related utility functions
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_JIT_UTIL_HPP
#define SPARSEX_INTERNALS_JIT_UTIL_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>

using namespace std;

namespace sparsex {
namespace jit {

string SourceFromFile(const char *filename);
string SourceFromFile(const string &filename);
void SourceToFile(const char *filename, const string &source);
void SourceToFile(const string &filename, const string &source);
const char *UniqueFilename(string &tmpl);
void RemoveFile(const char *path);
string Tabify(int num);

template<typename T>
string Stringify(const T& t)
{
    stringstream ss;
    ss << t;
    return ss.str();
}

} // end of namespace jit
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_JIT_UTIL_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
