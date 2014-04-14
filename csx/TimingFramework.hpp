/*
 * TimingFramework.hpp --  Timing utilities
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013       Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef TIMING_FRAMEWORK_HPP
#define TIMING_FRAMEWORK_HPP

#include "Timer.hpp"
#include <boost/unordered_map.hpp>

using namespace std;

namespace timing {

class TimingFramework
{
public:
    TimingFramework() {}

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

    void PrintAllStats(ostream &os) const
    {
        boost::unordered_map<const char *, Timer>::const_iterator it;

        os << "=== TIMING STATISTICS ===\n";
        for (it = timers_.cbegin(); it != timers_.cend(); ++it) {
            Timer tm = it->second;
            double t = tm.ElapsedTime();
            os << tm.GetDescription() << ": " << t << "\n";
        }
    }

private:
    boost::unordered_map<const char *, Timer> timers_;
    TimingFramework(const TimingFramework &); // Do not implement
    TimingFramework &operator=(const TimingFramework &); // Do not implement
};

} // end of namespace timing

#endif // TIMING_FRAMEWORK_HPP
