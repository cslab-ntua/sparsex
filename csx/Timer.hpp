/* -*- C++ -*-
 * 
 * Timer.hpp -- Timing framework
 *
 * Copyright (C) 2012-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2010-2012, Vasileios Karakasis
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef TIMER_HPP
#define TIMER_HPP

#include <sys/time.h>
#include <string>

namespace timing {

class Timer
{
public:

    Timer()
        : description_(),
          elapsed_time_(),
          timestamp_()
    {
        Clear();
    }

    Timer(std::string desc)
        : description_(desc), 
          elapsed_time_(), 
          timestamp_()  
    {
        Clear();
    }

    Timer(const char *desc, const char *desc2)
        : description_((std::string) desc), 
          elapsed_time_(), 
          timestamp_()
    {
        Clear();
    }

    ~Timer() {}

    void Start();
    void Pause();
    void Stop();
    void Clear();
    double ElapsedTime();

    void SetDescription(const char *desc)
    {
        description_ = (std::string) desc;
    }

    void SetDescription(std::string desc)
    {
        description_ = desc;
    }

    std::string GetDescription()
    {
        return description_;
    }

private:
    std::string description_;
    struct timeval elapsed_time_;
    struct timeval timestamp_;
};

} // end of namespace timing

#endif //TIMER_HPP
