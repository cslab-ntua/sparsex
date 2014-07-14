/*
 * \file Timer.cpp
 *
 * \brief Timing utilities
 *
 * Copyright (C) 2010-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2010-2012, Vasileios Karakasis
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/internals/Timer.hpp>
#include <stdlib.h>
#include <iostream>

using namespace std;

namespace sparsex {
namespace timing {

void Timer::Start()
{
    if (gettimeofday(&timestamp_, NULL) < 0) {
        cerr << "timer error: gettimeofday()" << endl;
        exit(1);
    }
}

void Timer::Pause()
{
    struct timeval t_stop;
    struct timeval t_interval;

    if (gettimeofday(&t_stop, NULL) < 0) {
        cerr << "timer error: gettimeofday()" << endl;
        exit(1);
    }

    timersub(&t_stop, &timestamp_, &t_interval);
    timeradd(&elapsed_time_, &t_interval, &elapsed_time_);
}

void Timer::Stop()
{
    Timer::Pause();
    Timer::Clear();
}

void Timer::Clear()
{
    timerclear(&elapsed_time_);
    timerclear(&timestamp_);
}

double Timer::ElapsedTime()
{
    return (elapsed_time_.tv_sec + 
            elapsed_time_.tv_usec / (double) 1000000);//USEC_PER_SEC
}

} // end of namespace timing
} // end of namespace sparsex
