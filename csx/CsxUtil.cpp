/* -*- C++ -*-
 *
 * CsxUtil.cc -- CSX-related routines.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "CsxUtil.hpp"

uint64_t MapSize(void *spm)
{
    unsigned int i;
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    spm_mt_thread_t *spm_thread;
    uint64_t size = 0;
    
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        size += spm_thread->map->length * sizeof(uint32_t);
        size += spm_thread->map->length * sizeof(uint32_t);
        size += spm_thread->map->length * sizeof(double);
    }
    
    return size;
}
