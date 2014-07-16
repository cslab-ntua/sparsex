/*
 * Copyright (C) 2007-2014, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Affinity.hpp
 * \brief Routine that sets the cpu affinity
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_AFFINITY_HPP
#define SPARSEX_INTERNALS_AFFINITY_HPP

#include <sparsex/internals/logger/Logger.hpp>
#include <stdlib.h>
#include <sched.h>

namespace sparsex {
namespace runtime {

void setaffinity_oncpu(unsigned int cpu);
int get_cpu_count();

} // end of namespace runtime
} // end of namespace sparsex

#endif // SPARSEX_INTERNALS_AFFINITY_HPP
