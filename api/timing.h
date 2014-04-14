/**
 * \file SparseX/timing.h -- \brief Timing routines.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sys/time.h>

struct timer {
    struct timeval elapsed_time;
    struct timeval timestamp;
};
typedef struct timer spx_timer_t;

static inline void spx_timer_clear(spx_timer_t *t)
{
    timerclear(&(t->elapsed_time));
    timerclear(&(t->timestamp));
}

static inline void spx_timer_start(spx_timer_t *t)
{
    if (gettimeofday(&(t->timestamp), NULL) < 0) {
        /* cerr << "timer error: gettimeofday()" << endl; */
        exit(1);
    }
}

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

static inline double spx_timer_get_secs(spx_timer_t *t)
{
    return (t->elapsed_time.tv_sec + t->elapsed_time.tv_usec / (double) 1000000);
}
