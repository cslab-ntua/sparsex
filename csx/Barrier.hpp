/*
 * Barrier.hpp --  A centralized barrier implementation.
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef BARRIER_HPP
#define BARRIER_HPP

#include <atomic>
#include <cstdlib>

std::atomic<int> cnt;
void central_barrier(bool& local_sense, size_t nr_threads);

#endif  // BARRIER_HPP
