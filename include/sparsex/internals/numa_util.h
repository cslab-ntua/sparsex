/*
 * numa_util.h -- NUMA utilitiy functions
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_NUMA_UTIL_H
#define SPARSEX_INTERNALS_NUMA_UTIL_H

#include "sparsex/internals/cdecl.h"
#include <stddef.h>

SPX_BEGIN_C_DECLS__
/*
 *  Custom interleaved allocation.
 *
 *  @size size in bytes of the full allocation (will be rounded to page size)
 *  @parts sizes in bytes of individual partitions (will be rounded to
 *      page size)
 *  @nr_parts  number of partitions
 *  @nodes memory nodes to bind each partition
 *
 *  @return pointer to the newly allocated area. Free with free_interleaved().
 */
void *alloc_interleaved(size_t size, size_t *parts, size_t nr_parts,
                        int *nodes);
void *alloc_onnode(size_t size, int node);
void free_interleaved(void *addr, size_t length);
int check_interleaved(void *addr, const size_t *parts,
                      size_t nr_parts, const int *nodes);
int check_region(void *addr, size_t size, int node);
void print_alloc_status(const char *data_descr, int err);

SPX_END_C_DECLS__

#endif  /* SPARSEX_INTERNALS_NUMA_UTIL_H */
