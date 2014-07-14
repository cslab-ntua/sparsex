/*
 * \file TimerCollection.cpp
 *
 * \brief A collection of arbitrary timers
 *
 * Copyright (C) 2009-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013       Athena Elafrou
 * Copyright (C) 2014       Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/internals/TimerCollection.hpp>

namespace sparsex {
namespace timing {

void TimerCollection::PrintAllTimers(ostream &os) const
{
    boost::unordered_map<const char *, Timer>::const_iterator it;

    for (it = timers_.cbegin(); it != timers_.cend(); ++it) {
        Timer tm = it->second;
        double t = tm.ElapsedTime();
        os << tm.GetDescription() << ": " << t << "\n";
    }
}

} // end of namespace timing
} // end of namespace sparsex
