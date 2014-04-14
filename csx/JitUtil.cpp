/*
 * JitUtil.cpp -- Implementation of jit-related utility functions.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "JitUtil.hpp"

#include <fstream>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

std::string SourceFromFile(const char *filename)
{
    std::ifstream ifs(filename);
    std::string ret;
    for(;;) {
        char c = ifs.get();
        if (!ifs.good())
            break;
        ret.push_back(c);
    }

    ifs.close();
    return ret;
}

void SourceToFile(const char *filename, const std::string &source)
{
    std::ofstream ofs(filename);
    ofs << source;
    ofs.close();
}

const char *UniqueFilename(std::string &tmpl)
{
    char local_tmpl[tmpl.size()+1];
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
    return tmpl.c_str();
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
