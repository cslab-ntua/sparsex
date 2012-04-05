/*
 * mmf.cc -- Matrix Market Format routines
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011,      Vasileios Karakasis
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include "mmf.h"

#include <fstream>
#include <cstdlib>

extern "C" {
#include <inttypes.h> /* PRIu64 */
}

using namespace csx;

namespace csx {

void getMmfHeader(const char *mmf_file, uint64_t &nrows,
                  uint64_t &ncols, uint64_t &nnz)
{
    std::ifstream in;

    in.open(mmf_file);
    getMmfHeader(in, nrows, ncols, nnz);
    in.close();
}

void getMmfHeader(std::istream &in, uint64_t &nrows, uint64_t &ncols,
                  uint64_t &nnz)
{
    char buff[512];
    int ret;

    do {
        in.getline(buff, sizeof(buff));
    } while (buff[0] == '#');

    ret = sscanf(buff, "%" PRIu64 " %" PRIu64 " %" PRIu64,
                 &nrows, &ncols, &nnz);
    if (ret != 3) {
        std::cerr << "mmf header error: sscanf" << std::endl;
        exit(1);
    }
}

} // end csx namespace

MMF::MMF(std::istream &in) : in_(in)
{
    getMmfHeader(in_, this->nrows, this->ncols, this->nnz);
}

bool MMF::next(uint64_t &y, uint64_t &x, double &v)
{
    char buff[512];
    int ret;

    if (in_.getline(buff, sizeof(buff)).eof()) {
        return false;
    }

    ret = sscanf(buff, "%" PRIu64 " %" PRIu64 " %lf", &y, &x, &v);
    assert(ret == 3);
    return true;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
