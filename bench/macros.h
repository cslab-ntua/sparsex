/* -*- C++ -*-
 *
 * macros.h --  Useful macros.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef MACROS_H__
#define MACROS_H__

#include "timer.h"

/* Usefull macros */
#define BENCH(timer, loops, outer_loops, code)          \
    timer.Clear();                                      \
    timer.Start();                                      \
    for (unsigned int i = 0; i < outer_loops; i++) {    \
        for (unsigned long int j = 0; j < loops; j++) { \
            code;                                       \
        }                                               \
    }                                                   \
    timer.Pause();

#define SPMV_BENCH(code) BENCH(t, LOOPS, OUTER_LOOPS, code);

#endif  // MACROS_H__
