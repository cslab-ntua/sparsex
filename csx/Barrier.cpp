/*
 * Barrier.cpp --  A centralized barrier implementation.
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "Barrier.hpp"

using namespace std;

static atomic<bool> sense(true);

void central_barrier(bool& local_sense, size_t nr_threads)
{
    local_sense = !local_sense; // each processor toggles its own sense
    if (atomic_fetch_sub(&cnt, 1) == 1) {
        cnt = nr_threads;   // atomic store?
        sense.store(local_sense); // last processor toggles global sense 
    } else {
        while (local_sense != sense.load()) {}
    }
}
