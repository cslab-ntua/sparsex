/* -*- C++ -*-
 *
 * affinity.h -- Routine that sets the cpu affinity.
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef AFFINITY_H__
#define AFFINITY_H__

#include "logger.hpp"

#include <stdlib.h>
#include <sched.h>

void setaffinity_oncpu(unsigned int cpu);

#endif /* AFFINITY_H__ */
