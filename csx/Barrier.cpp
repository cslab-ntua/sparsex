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

static atomic<bool> global_sense(true);
atomic<int> barrier_cnt;

void central_barrier(bool *local_sense, size_t nr_threads)
{
    *local_sense = !(*local_sense); // each processor toggles its own sense
    if (atomic_fetch_sub(&barrier_cnt, 1) == 1) {
        barrier_cnt = nr_threads;   // atomic store?
        global_sense.store(*local_sense); // last processor toggles global sense 
    } else {
        while ((*local_sense) != global_sense.load()) {}
    }
}
