/*
 * numa_util.h -- NUMA utilitiy functions
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef __NUMA_UTIL_H__
#define __NUMA_UTIL_H__

#include <stddef.h>

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
                        const int *nodes);
void free_interleaved(void *addr, size_t length);
int check_interleaved(void *addr, size_t size, size_t *parts, size_t nr_parts,
                       const int *nodes);
int check_onnode(void *addr, size_t size, const int node);

#endif  /* __NUMA_UTIL_H__ */
