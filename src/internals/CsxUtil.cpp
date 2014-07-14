/*
 * \file CsxUtil.cpp
 *
 * \breif CSX-related utilities
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/internals/CsxUtil.hpp>

namespace sparsex {
namespace csx {

const unsigned long DeltaIdOffset = 0;
const unsigned long PatternIdOffset = 10000;

unsigned long GetPatternId(size_t delta_size)
{
    return (delta_size << 3) + DeltaIdOffset;
}

} // end of namespace csx
} // end of namespace sparsex
