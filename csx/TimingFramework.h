/* -*- C++ -*-
 *
 * TimingFramework.h --  Timing utilities
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013       Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef TIMINGFRAMEWORK_H__
#define TIMINGFRAMEWORK_H__

#include "timer.h"
#include <boost/unordered_map.hpp>

using namespace std;

namespace csx {

    enum AvailableTimers {
        // Preprocessing related timers
        PRE_TIMER_TOTAL = 0,
        PRE_TIMER_STATS,
        PRE_TIMER_ENCODE,
        PRE_TIMER_ALLOC,
        // SpmV related timers
        TIMER_END
    };

}

using namespace csx;

class TimingFramework
{
public:

    TimingFramework() 
    {
        descriptions_[PRE_TIMER_TOTAL] = "Total Time: ";
        descriptions_[PRE_TIMER_STATS] = "Stats Time: ";
        descriptions_[PRE_TIMER_ENCODE] = "Encode Time: ";
        descriptions_[PRE_TIMER_ALLOC] = "Alloc Time: ";
    }

    void EnableTimer(AvailableTimers type)
    {
        csx::Timer t(descriptions_[type]);
        timers_.emplace(type, t);
    }

    void DisableTimer(AvailableTimers type)
    {
        if (timers_.find(type) != timers_.end())
            timers_.erase(type);
    }

    bool Enabled(AvailableTimers type)
    {
        if (timers_.find(type) != timers_.end())
            return true;
        else
            return false;
    }

    void PrintStats(AvailableTimers type, ostream &os)
    {
        if (timers_.find(type) != timers_.end())
            os << descriptions_[type] << timers_[type].ElapsedTime() << endl;
        else
            os << "Timer doesn't exist!" << endl;
    }

    void PrintAllStats(ostream &os) const
    {
        double t;
        boost::unordered_map<AvailableTimers, csx::Timer>::const_iterator it;
        csx::Timer timer;

        os << "Timing Statistics" << endl;
        os << "-----------------" << endl;
        for (it = timers_.cbegin(); it != timers_.cend(); ++it) {
            timer = (*it).second;
            t = timer.ElapsedTime();
            os << descriptions_[(*it).first] << t << endl;
        }
    }

    boost::unordered_map<AvailableTimers, csx::Timer> timers_;

private:
    string descriptions_[TIMER_END];
};

#endif // TIMINGFRAMEWORK_H__
