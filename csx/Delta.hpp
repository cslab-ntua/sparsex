/* -*- C++ -*-
 *
 * Delta.hpp -- Utilities for delta units.
 *
 * Copyright (C) 2009-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2013-2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef DELTA_HPP
#define DELTA_HPP

#include <inttypes.h>

namespace csx {

#define PID_DELTA_BASE  0UL

typedef union {
    uint8_t d8;
    uint16_t d16;
    uint32_t d32;
    uint64_t d64;
} DeltaType;

size_t GetDeltaSize(size_t val)
{
    DeltaType du;
    if ( (du.d8 = val) == val)
        return sizeof(du.d8);

    if ( (du.d16 = val) == val)
        return sizeof(du.d16);

    if ( (du.d32 = val) == val)
        return sizeof(du.d32);

    return sizeof(du.d64);
}

unsigned long GetDeltaPatternId(size_t delta_size)
{
    return (delta_size << 3) + PID_DELTA_BASE;
}

} // end csx namespace

#endif // DELTA_HPP
