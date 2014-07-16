/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * Copyright (C) 2014,      Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file timing.h
 * \brief Timing routines.
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_TIMING_H
#define SPARSEX_TIMING_H

#include <sys/time.h>

struct timer {
    struct timeval elapsed_time;
    struct timeval timestamp;
};
typedef struct timer spx_timer_t;

/**
 *  Initializes a timer object.
 *
 *  @param[in] t       timer object to be initialized.
 */
static inline void spx_timer_clear(spx_timer_t *t)
{
    timerclear(&(t->elapsed_time));
    timerclear(&(t->timestamp));
}

/**
 *  Starts a timer object.
 *
 *  @param[in] t       timer object to be launched.
 */
static inline void spx_timer_start(spx_timer_t *t)
{
    if (gettimeofday(&(t->timestamp), NULL) < 0) {
        /* cerr << "timer error: gettimeofday()" << endl; */
        exit(1);
    }
}

/**
 *  Pauses a timer object.
 *
 *  @param[in] t       timer object to be paused.
 */
static inline void spx_timer_pause(spx_timer_t *t)
{
    struct timeval t_stop;
    struct timeval t_interval;

    if (gettimeofday(&t_stop, NULL) < 0) {
        /* cerr << "timer error: gettimeofday()" << endl; */
        exit(1);
    }

    timersub(&t_stop, &(t->timestamp), &t_interval);
    timeradd(&(t->elapsed_time), &t_interval, &(t->elapsed_time));
}

/**
 *  Returns the elapsed time in seconds.
 *
 *  @param[in] t        a timer object.
 *  @return             elapsed seconds.
 */
static inline double spx_timer_get_secs(spx_timer_t *t)
{
    return t->elapsed_time.tv_sec +
        t->elapsed_time.tv_usec / (double) 1000000;
}

#endif /* SPARSEX_TIMING_H */
