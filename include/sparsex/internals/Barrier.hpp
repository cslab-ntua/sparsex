/*
 * Barrier.hpp --  A centralized barrier with timeout.
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 *
 * In a centralized barrier, every thread spins on a global “sense” flag. A
 * global count variable holds the number of threads. The global sense flag is
 * initially set to an initial state to indicate that the barrier hasn’t been 
 * reached. When a thread reaches the barrier, it decrements the global count to
 * indicate that it has reached the barrier. The last thread to reach the barrier
 * resets count to the number of threads, and toggles the global “sense”, which
 * indicates to the other threads spinning (or sleeping) that all threads have
 * reached the barrier. At this point, all threads stop spinning and continue
 * execution.
 */

#ifndef SPARSEX_INTERNALS_BARRIER_HPP
#define SPARSEX_INTERNALS_BARRIER_HPP

#include <atomic>
#include <cstdlib>

void centralized_barrier(int *local_sense, size_t nr_threads);

#endif  // SPARSEX_INTERNALS_BARRIER_HPP
