/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Barrier.hpp
 * \brief A centralized barrier with timeout
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_BARRIER_HPP
#define SPARSEX_INTERNALS_BARRIER_HPP

#include <atomic>
#include <cstdlib>

namespace sparsex {
namespace runtime {

/*
 * In a centralized barrier, every thread spins on a global “sense” flag. A
 * global count variable holds the number of threads. The global sense flag is
 * initially set to an initial state to indicate that the barrier hasn’t been 
 * reached. When a thread reaches the barrier, it decrements the global count to
 * indicate that it has reached the barrier. The last thread to reach the
 * barrier resets count to the number of threads, and toggles the global
 * “sense”, which indicates to the other threads spinning (or sleeping) that
 * all threads have reached the barrier. At this point, all threads stop
 * spinning and continue execution.
 */
void centralized_barrier(int *local_sense, size_t nr_threads);

} // end of namespace runtime
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_BARRIER_HPP
