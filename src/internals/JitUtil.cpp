/*
 * JitUtil.cpp -- Implementation of jit-related utility functions.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "sparsex/internals/JitUtil.hpp"

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>

using namespace std;

string SourceFromFile(const char *filename)
{
    ifstream ifs(filename);
    string ret;
    for(;;) {
        char c = ifs.get();
        if (!ifs.good())
            break;
        ret.push_back(c);
    }

    ifs.close();
    return ret;
}

string SourceFromFile(const string &filename)
{
    return SourceFromFile(filename.c_str());
}

void SourceToFile(const char *filename, const string &source)
{
    ofstream ofs(filename);
    ofs << source;
    ofs.close();
}

void SourceToFile(const string &filename, const string &source)
{
    SourceToFile(filename.c_str(), source);
}

void RemoveFile(const char *path)
{
    if (remove(path) < 0) {
        perror("RemoveFile()");
        exit(1);
    }
}

const char *UniqueFilename(string &tmpl)
{
    // Cannot play with C++ strings here, since mkstemp() alters is char *
    // argument
    char *local_tmpl = new char[tmpl.size()+1];
    int fd;

    tmpl.copy(local_tmpl, tmpl.size());
    local_tmpl[tmpl.size()] = '\0';

    fd = mkstemp(local_tmpl);
    if (fd < 0) {
        perror("mkstemp() failed");
        exit(1);
    }

    // close file descriptor and let the caller open the file in the C++ way
    // FIXME: an implementation with boost iostreams would be more portable
    close(fd);
    tmpl.assign(local_tmpl);
    delete[] local_tmpl;
    return tmpl.c_str();
}

string Tabify(int num)
{
    stringstream ss;
    
    for (int i = 0; i < num; i++)
        ss << "\t";
    
    return ss.str();
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
