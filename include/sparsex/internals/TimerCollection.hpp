/*
 * \file TimerCollection.hpp
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

#ifndef SPARSEX_INTERNALS_TIMING_FRAMEWORK_HPP
#define SPARSEX_INTERNALS_TIMING_FRAMEWORK_HPP

#include <sparsex/internals/Timer.hpp>
#include <boost/unordered_map.hpp>
#include <iostream>

using namespace std;

namespace timing {

class TimerCollection
{
public:
    TimerCollection() {}

    void CreateTimer(const char *name, const char *description)
    {
        timers_.emplace(name, Timer(description));
    }

    void StartTimer(const char *name)
    {
        timers_[name].Start();
    }

    void PauseTimer(const char *name)
    {
        timers_[name].Pause();
    }

    void PrintAllTimers(ostream &os) const;

private:
    boost::unordered_map<const char *, Timer> timers_;
    TimerCollection(const TimerCollection &); // Do not implement
    TimerCollection &operator=(const TimerCollection &); // Do not implement
};

} // end of namespace timing

#endif  // SPARSEX_INTERNALS_TIMING_FRAMEWORK_HPP
